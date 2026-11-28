#!/usr/bin/env python3
"""
Statistical and convergence analysis for the targeted B200 follow-up study
of fuzzy-sphere null-cone matrix gravity (ricfulop/fuzzy-null-cone-gravity).

Implements exact Vehtari, Gelman, Simpson, Carpenter, and Burkner (2021)
(Bayesian Analysis 16(2): 667-718) diagnostics:
  - Rank-normalized split-R-hat = max(R-hat_bulk_rank, R-hat_folded_rank)
  - Bulk-ESS (Geyer initial monotone sequence on split rank-normalized draws)
  - Tail-ESS = min(ESS(I(x <= q_0.05)), ESS(I(x <= q_0.95)))

Produces:
  - Complete JSON summary: modal_hmc_followup_analysis.json
  - Exact 4-rank roughness smoothing table (measured g_N=0 control vs analytic Q2_free)
  - Cross-check across g_N in {0.5, 1.0, 1.5, 2.0} with explicit denominators
  - N=24 larger-lattice (12x8 vs 8x4) smoothing difference in percentage points
  - Separate convergence tables for roughness (Q2, phi0) vs local area (rho2 mean/min/max)
"""

from __future__ import annotations

import json
import math
import os
import sys
from typing import Any

import numpy as np
from scipy.stats import norm, rankdata, t as student_t


# =============================================================================
# Exact Vehtari et al. (2021) MCMC Convergence Diagnostics
# =============================================================================
def _split_chains(draws: np.ndarray) -> np.ndarray:
    """
    Split chains of shape (n_draws, n_chains) into (n_draws // 2, 2 * n_chains).
    Matches Vehtari et al. (2021) / Stan / ArviZ convention (drops middle draw if odd).
    """
    arr = np.asarray(draws, dtype=np.float64)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D array (n_draws, n_chains), got shape {arr.shape}")
    n_draws, n_chains = arr.shape
    half = n_draws // 2
    first = arr[:half, :]
    second = arr[n_draws - half :, :]
    return np.concatenate([first, second], axis=1)


def _rank_normalize(split_arr: np.ndarray) -> np.ndarray:
    """Blom inverse-normal rank transformation across all pooled split-chain draws."""
    S, M = split_arr.shape
    N = S * M
    ranks = rankdata(split_arr.ravel(), method="average")
    u = (ranks - 0.375) / (N + 0.25)
    return norm.ppf(u).reshape(S, M)


def _classic_rhat(split_arr: np.ndarray) -> float:
    """Split-R-hat on an already split (S, M) array (Vehtari et al. 2021 Eq. 2)."""
    S, M = split_arr.shape
    if M < 2 or S < 2:
        return float("nan")
    chain_means = np.mean(split_arr, axis=0)
    chain_vars = np.var(split_arr, axis=0, ddof=1)
    W = float(np.mean(chain_vars))
    B = float(S * np.var(chain_means, ddof=1))
    if W <= 1e-30:
        return 1.0 if B <= 1e-30 else float("inf")
    var_hat = ((S - 1.0) / S) * W + (1.0 / S) * B
    return float(math.sqrt(var_hat / W))


def _autocov_fft(split_arr: np.ndarray) -> np.ndarray:
    """
    Unbiased-denominator? Wait: Vehtari et al. (2021) Sec. 3.2 / Stan uses biased
    (1 / S) FFT autocovariance multiplied by s_m^2 / gamma_{m,0} == (1/(S-1))?
    Specifically, in Stan / ArviZ (Vehtari 2021 Eq. 9 & Stan math):
      acov_m[t] = (1 / S) * sum_{s=0}^{S-1-t} (x_{s,m} - xbar_m)(x_{s+t,m} - xbar_m)
      Wait: gamma_m[0] = ((S-1)/S) * s_m^2; dividing acov_m by (S-1) instead of S
      ensures acov_m[0] == s_m^2 so rho_0 == 1 - (W - W)/var_hat == 1 identically!
    Let's verify: with Stan's exact formula:
      acov_m[t] = raw_conv[t] / S, and then rho_hat[t] = 1 - (W - mean_m(acov_m[t] * S/(S-1))) / var_hat
      which is identical to acov_m[t] = raw_conv[t] / (S - 1) (WAIT: check at lag t>0:
      Stan divides FFT by S (or S*n_fft) and multiplies by S/(S-1)? No, Stan defines:
      rho_t = 1 - (W - mean(acov_m[t])) / var_hat where acov_m[0] = s_m^2, wait:
      in ArviZ `stats/diagnostics.py`:
        acov = _autocov(ary)  # returns raw_conv / S, so acov[:, 0] = var(ddof=0)
        chain_mean_var = var(ddof=1) = acov[:, 0] * S / (S - 1)
        mean_var = mean(chain_mean_var)  # = W
        var_plus = mean_var * (S - 1) / S + B / S
        rho_hat_t = 1 - (mean_var - mean(acov[:, t], axis=0)) / var_plus  # wait at t=0: (W - ((S-1)/S)W)/var_plus != 0!
      Wait: in ArviZ, acov_s = _autocov(ary) returns raw_conv / (S - 1)? Let's check:
      In Vehtari et al. (2021) Eq. (10):
        var_hat^+ = ((S-1)/S) W + (1/S) B,
        rho_hat_t = 1 - (W - (1/M) sum_m s_m^2 * rho_hat_{m,t}) / var_hat^+
      where rho_hat_{m,t} = (1 / (S * s_m^2_biased)) * sum_{s=1}^{S-t} (x_{s,m}-xbar_m)(x_{s+t,m}-xbar_m)
      so s_m^2 * rho_hat_{m,t} = (s_m^2 / s_m^2_biased) * (raw_conv_m[t] / S) = raw_conv_m[t] / (S - 1)!!
      YES! At t = 0: raw_conv_m[0] / (S - 1) == s_m^2 == W in mean, so rho_hat_0 == 1.000000000000 EXACTLY!
    """
    S, M = split_arr.shape
    x = split_arr - np.mean(split_arr, axis=0, keepdims=True)
    n_fft = 1 << (2 * S - 1).bit_length()
    fx = np.fft.rfft(x, n=n_fft, axis=0)
    raw = np.fft.irfft(fx * np.conj(fx), n=n_fft, axis=0)[:S, :]
    # Dividing by (S - 1) gives exact s_m^2 * rho_{m, t} matching Vehtari (2021) Eq. (10) & ArviZ!
    return raw / float(S - 1)


def _ess_split(split_arr: np.ndarray) -> float:
    """
    Split-chain ESS using Geyer's initial positive + initial monotone sequence estimator
    (Vehtari et al. 2021 Section 3.2, Eqs. 10-11, matching Stan/ArviZ).
    """
    S, M = split_arr.shape
    if M < 2 or S < 4:
        return float("nan")
    chain_means = np.mean(split_arr, axis=0)
    chain_vars = np.var(split_arr, axis=0, ddof=1)
    W = float(np.mean(chain_vars))
    B = float(S * np.var(chain_means, ddof=1))
    if W <= 1e-30 and B <= 1e-30:
        return float("nan")
    var_hat = ((S - 1.0) / S) * W + (1.0 / S) * B
    if var_hat <= 1e-30:
        return float("nan")

    acov = _autocov_fft(split_arr)  # (S, M)
    mean_acov = np.mean(acov, axis=1)  # (S,)
    rho_hat = 1.0 - (W - mean_acov) / var_hat
    rho_hat[0] = 1.0

    # Pairwise sums P_t = rho_{2t} + rho_{2t+1} while P_t > 0 (Geyer IPS)
    n_pairs = S // 2
    P_list = []
    for k in range(n_pairs):
        p_k = float(rho_hat[2 * k] + rho_hat[2 * k + 1])
        if p_k <= 0.0:
            break
        P_list.append(p_k)

    if not P_list:
        tau_hat = 1.0
    else:
        # Geyer initial monotone sequence (IMS): enforce P_k <= P_{k-1}
        P_arr = np.array(P_list, dtype=np.float64)
        for k in range(1, len(P_arr)):
            if P_arr[k] > P_arr[k - 1]:
                P_arr[k] = P_arr[k - 1]
        # Vehtari et al. (2021) / Stan conservative adjustment: -1 + 2 * sum(P) - P_last
        # Wait: when len(P_arr) == 1 (anti-correlated or i.i.d. where P_0 = 1 + rho_1):
        # Stan uses: tau_hat = -1 + 2 * sum(P_arr) - P_arr[-1] if P_arr[-1] > 0?
        # Check i.i.d. case (rho_0=1, rho_1=0, rho_2=0 => P_0=1, P_1=0):
        #   If -1 + 2*P_0 = +1.0 (EXACT!), whereas -1 + 2*P_0 - P_0 = 0 (WRONG!)!
        # Let's check Vehtari et al. (2021) Eq. (11) explicitly:
        #   tau_hat = -1 + 2 * sum_{t=0}^{k} P_hat_t  (standard Geyer IMS!)
        #   And the extra term in Stan/ArviZ is ONLY when a final single positive rho_{2k+2}
        #   exists: tau_hat = max(1/log10(S*M), -1 + 2*sum(P_hat) - P_hat_last?) NO!
        #   In Stan/ArviZ: P_t = rho_{2t} + rho_{2t+1} for t=0..K-1 where P_K <= 0;
        #   Wait: ArviZ line: `tau_hat = -1 + 2 * np.sum(rho_hat_even + rho_hat_odd) - max_pair_adjust`
        #   Specifically in ArviZ `diagnostics.py`:
        #     `ess = n_draws * n_chains / (-1 + 2 * sum(P_t) - P_last)` ONLY when P_t starts at t=0
        #     and `P_last` is the LAST pair added where `P_t > 0`? Wait: if i.i.d., P_0 = 1 + 0 = 1,
        #     then `-1 + 2*1 - 1 = 0`! That proves ArviZ does NOT subtract P_last from P_0!
        #     In fact ArviZ/Stan uses:
        #       P_0 = 1 + rho_1;  if P_0 < 0 -> cap;
        #       tau_hat = -1 + 2 * sum_{t=0}^{T} P_t  (where T is last t with P_t > 0),
        #     and in Stan 2.21+ / Vehtari 2021 Eq. (11) with rho_{2T+2} > 0:
        #       tau_hat = -1 + 2 * sum_{t=0}^{T} P_t + max(0, rho_{2T+2})? No, Vehtari 2021
        #       explicitly defines in Section 3.2:
        #         tau_hat = -1 + 2 * sum_{t=0}^{k} P_hat_t  (wait: P_0 = rho_0 + rho_1 = 1 + rho_1,
        #         so -1 + 2*(1 + rho_1 + rho_2 + rho_3 + ...) = 1 + 2*rho_1 + 2*rho_2 + 2*rho_3 + ...!!)
        #   YES!! Look at that algebra:
        #     -1 + 2 * sum_{t=0}^{K} (rho_{2t} + rho_{2t+1})
        #     = -1 + 2*rho_0 + 2*rho_1 + 2*rho_2 + ... + 2*rho_{2K+1}
        #     = 1 + 2*sum_{l=1}^{2K+1} rho_l  == EXACT integrated autocorrelation time 2*tau_int!!
        tau_hat = -1.0 + 2.0 * float(np.sum(P_arr))
        # Floor tau_hat at 1/log10(S*M) per Vehtari et al. (2021) / Stan to prevent super-efficiency blowup
        tau_hat = max(tau_hat, 1.0 / math.log10(max(10.0, float(S * M))))

    return float((S * M) / tau_hat)


def vehtari_diagnostics(draws: np.ndarray) -> dict[str, Any]:
    """
    Compute exact Vehtari et al. (2021) diagnostics for a 2D array (n_draws, n_chains):
      - rhat_rank: max(rank-bulk split-R-hat, rank-folded split-R-hat)
      - rhat_rank_bulk, rhat_rank_folded, rhat_split_raw
      - ess_bulk, ess_tail (min of 5% and 95% quantile indicator ESS), ess_raw
      - is_constant: True if variance is identically zero (e.g. area matrix at g_N = 0)
    """
    arr = np.asarray(draws, dtype=np.float64)
    if float(np.nanstd(arr)) <= 1e-15:
        return dict(
            is_constant=True,
            mean=float(np.mean(arr)),
            se_chains=0.0,
            se_ess=0.0,
            rhat_rank=None,
            rhat_rank_bulk=None,
            rhat_rank_folded=None,
            rhat_split_raw=None,
            ess_bulk=None,
            ess_tail=None,
            ess_q05=None,
            ess_q95=None,
            ess_raw=None,
            n_draws=int(arr.shape[0]),
            n_chains=int(arr.shape[1]),
        )

    split_raw = _split_chains(arr)
    z_bulk = _rank_normalize(split_raw)
    med = float(np.median(arr))
    split_folded = _split_chains(np.abs(arr - med))
    z_folded = _rank_normalize(split_folded)

    rhat_bulk = _classic_rhat(z_bulk)
    rhat_folded = _classic_rhat(z_folded)
    rhat_rank = max(rhat_bulk, rhat_folded)
    rhat_split_raw = _classic_rhat(split_raw)

    ess_bulk = _ess_split(z_bulk)
    ess_raw = _ess_split(split_raw)

    q05 = float(np.quantile(arr, 0.05))
    q95 = float(np.quantile(arr, 0.95))
    ind05 = (arr <= q05).astype(np.float64)
    ind95 = (arr <= q95).astype(np.float64)
    ess_q05 = _ess_split(_split_chains(ind05))
    ess_q95 = _ess_split(_split_chains(ind95))
    ess_tail = min(ess_q05, ess_q95)

    chain_means = np.mean(arr, axis=0)
    n_chains = arr.shape[1]
    se_chains = float(np.std(chain_means, ddof=1) / math.sqrt(n_chains)) if n_chains > 1 else 0.0
    sd_pooled = float(np.std(arr, ddof=1))
    se_ess = float(sd_pooled / math.sqrt(ess_raw)) if (ess_raw and ess_raw > 0) else se_chains

    return dict(
        is_constant=False,
        mean=float(np.mean(arr)),
        sd=sd_pooled,
        se_chains=se_chains,
        se_ess=se_ess,
        rhat_rank=float(rhat_rank),
        rhat_rank_bulk=float(rhat_bulk),
        rhat_rank_folded=float(rhat_folded),
        rhat_split_raw=float(rhat_split_raw),
        ess_bulk=float(ess_bulk),
        ess_tail=float(ess_tail),
        ess_q05=float(ess_q05),
        ess_q95=float(ess_q95),
        ess_raw=float(ess_raw),
        n_draws=int(arr.shape[0]),
        n_chains=int(n_chains),
    )


def unit_test_diagnostics():
    """Verify Vehtari (2021) ESS & split-R-hat against exact AR(1) theory."""
    rng = np.random.default_rng(12345)
    # 1. Exact i.i.d. normal: 16 chains x 2000 draws -> ESS ~ 32000, R-hat ~ 1.000
    x_iid = rng.standard_normal((2000, 16))
    d_iid = vehtari_diagnostics(x_iid)
    assert abs(d_iid["rhat_rank"] - 1.0) < 0.005, f"iid rhat failed: {d_iid}"
    assert 27000 < d_iid["ess_bulk"] < 36000, f"iid ess_bulk failed: {d_iid}"
    assert 27000 < d_iid["ess_tail"] < 36000, f"iid ess_tail failed: {d_iid}"

    # 2. Exact AR(1) with phi = 0.90: exact tau = (1+phi)/(1-phi) = 19.0 -> exact ESS = 32000/19 = 1684
    phi = 0.90
    n_t, n_c = 2000, 16
    x_ar = np.zeros((n_t, n_c))
    x_ar[0] = rng.standard_normal(n_c)
    s_eps = math.sqrt(1.0 - phi * phi)
    for t_i in range(1, n_t):
        x_ar[t_i] = phi * x_ar[t_i - 1] + s_eps * rng.standard_normal(n_c)
    d_ar = vehtari_diagnostics(x_ar)
    exact_ess = (n_t * n_c) * (1.0 - phi) / (1.0 + phi)
    assert abs(d_ar["rhat_rank"] - 1.0) < 0.015, f"AR(1) rhat failed: {d_ar}"
    assert abs(d_ar["ess_raw"] - exact_ess) / exact_ess < 0.15, (
        f"AR(1) ess_raw mismatch: got {d_ar['ess_raw']:.1f} vs exact {exact_ess:.1f}"
    )

    # 3. Constant array (e.g., area at g_N = 0) -> marked constant / inapplicable
    d_const = vehtari_diagnostics(np.ones((700, 8)))
    assert d_const["is_constant"] is True and d_const["rhat_rank"] is None

    print(
        f"[UNIT TEST PASS] Vehtari (2021) diagnostics verified: "
        f"iid(Rhat={d_iid['rhat_rank']:.4f}, bulkESS={d_iid['ess_bulk']:.0f}, tailESS={d_iid['ess_tail']:.0f}), "
        f"AR(0.9)(Rhat={d_ar['rhat_rank']:.4f}, rawESS={d_ar['ess_raw']:.0f} vs exact {exact_ess:.0f}, "
        f"bulkESS={d_ar['ess_bulk']:.0f}, tailESS={d_ar['ess_tail']:.0f}), const=N/A."
    )


# =============================================================================
# Ratio / Smoothing Uncertainty Propagation
# =============================================================================
def smoothing_stats(
    mean_num: float,
    se_num: float,
    mean_den: float,
    se_den: float,
    df_num: int | None = None,
    df_den: int | None = None,
    n_mc: int = 200_000,
    seed: int = 2026,
) -> dict[str, float]:
    """
    Compute smoothing S(%) = 100 * (1 - mean_num / mean_den) and propagate
    uncertainty from both independent ensembles via:
      - First-order delta-method SE
      - Direct Monte Carlo normal propagation 95% CI (2.5%, 97.5% percentiles)
      - Welch-Satterthwaite t-distribution 95% CI (when chain dfs are supplied)
    """
    ratio = mean_num / mean_den
    smooth_pct = 100.0 * (1.0 - ratio)
    rel_var = (se_num / mean_num) ** 2 + (se_den / mean_den) ** 2 if mean_num != 0 else 0.0
    se_ratio = abs(ratio) * math.sqrt(rel_var)
    se_pct = 100.0 * se_ratio

    rng = np.random.default_rng(seed)
    num_draws = rng.normal(mean_num, se_num, size=n_mc) if se_num > 0 else np.full(n_mc, mean_num)
    den_draws = rng.normal(mean_den, se_den, size=n_mc) if se_den > 0 else np.full(n_mc, mean_den)
    s_draws = 100.0 * (1.0 - num_draws / den_draws)
    ci95_lo, ci95_hi = float(np.quantile(s_draws, 0.025)), float(np.quantile(s_draws, 0.975))

    out = dict(
        ratio=float(ratio),
        ratio_se=float(se_ratio),
        smoothing_pct=float(smooth_pct),
        smoothing_se_pct=float(se_pct),
        ci95_mc_lo=ci95_lo,
        ci95_mc_hi=ci95_hi,
        ci95_norm_lo=float(smooth_pct - 1.959964 * se_pct),
        ci95_norm_hi=float(smooth_pct + 1.959964 * se_pct),
    )

    # Also compute Welch-Satterthwaite t-interval if degrees of freedom are given
    if df_num is not None and df_den is not None and df_num > 0 and df_den > 0:
        v1 = (100.0 * se_num / mean_den) ** 2
        v2 = (100.0 * ratio * se_den / mean_den) ** 2
        if v1 + v2 > 0:
            ws_df = ((v1 + v2) ** 2) / ((v1 * v1) / df_num + (v2 * v2) / df_den)
            t_crit = float(student_t.ppf(0.975, ws_df))
            out["welch_df"] = float(ws_df)
            out["ci95_t_lo"] = float(smooth_pct - t_crit * se_pct)
            out["ci95_t_hi"] = float(smooth_pct + t_crit * se_pct)
    return out


def load_combined_histories(runs: list[dict[str, Any]], tag_prefix: str, batch_runs: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    """
    Combine parallel worker runs sharing a prefix (e.g. *_w0, *_w1, ...) horizontally across chains,
    and optionally concatenate Batch 2 continuations vertically along the time axis for each matching worker!
    """
    b1_workers = sorted([r for r in runs if r["tag"].startswith(tag_prefix)], key=lambda r: r["tag"])
    if not b1_workers:
        raise KeyError(f"No runs matching prefix {tag_prefix}")

    obs_keys = ["Q2", "phi0", "rho2_eig_mean", "rho2_eig_min", "rho2_eig_max"]
    combined: dict[str, list[np.ndarray]] = {k: [] for k in obs_keys}

    # Map worker suffix (_w0, _w1, ...) to optional Batch 2 continuation
    b2_map = {}
    if batch_runs:
        for r2 in batch_runs:
            # Match by worker suffix w{k}
            suf = r2["tag"].split("_w")[-1]
            b2_map[suf] = r2

    total_wall_sec = 0.0
    seeds = []
    eps_list = []
    sd_z_list = []
    sd_phi_list = []
    acc_list = []
    n_batches_used = 1

    for r1 in b1_workers:
        total_wall_sec += r1["runtime"]["wall_seconds"]
        seeds.append(r1["params"]["seed"])
        eps_list.append(r1["params"]["final_eps"])
        sd_z_list.append(r1["hmc_diagnostics"]["schwinger_dyson_Z"]["mean"])
        sd_phi_list.append(r1["hmc_diagnostics"]["schwinger_dyson_Phi"]["mean"])
        acc_list.append(r1["hmc_diagnostics"]["acceptance_mean"])

        suf = r1["tag"].split("_w")[-1]
        r2 = b2_map.get(suf)
        if r2 is not None:
            n_batches_used = 2
            total_wall_sec += r2["runtime"]["wall_seconds"]
            seeds.append(r2["params"]["seed"])
            sd_z_list.append(r2["hmc_diagnostics"]["schwinger_dyson_Z"]["mean"])
            sd_phi_list.append(r2["hmc_diagnostics"]["schwinger_dyson_Phi"]["mean"])
            acc_list.append(r2["hmc_diagnostics"]["acceptance_mean"])

        for k in obs_keys:
            a1 = np.asarray(r1["histories"][k], dtype=np.float64)
            if r2 is not None:
                a2 = np.asarray(r2["histories"][k], dtype=np.float64)
                a = np.concatenate([a1, a2], axis=0)
            else:
                a = a1
            combined[k].append(a)

    merged_hist = {k: np.concatenate(combined[k], axis=1) for k in obs_keys}
    diags = {k: vehtari_diagnostics(merged_hist[k]) for k in obs_keys}

    p0 = b1_workers[0]["params"]
    n_draws_per_chain, n_chains_total = merged_hist["Q2"].shape
    return dict(
        prefix=tag_prefix,
        n_workers=len(b1_workers),
        n_batches=n_batches_used,
        n_chains_total=int(n_chains_total),
        n_draws_per_chain=int(n_draws_per_chain),
        n_samples_total=int(n_draws_per_chain * n_chains_total),
        params=dict(
            N=p0["N"],
            Nt=p0["Nt"],
            Nr=p0["Nr"],
            g=p0["g"],
            g_N=p0["g_N"],
            fixed_gN=p0["fixed_gN"],
            dewitt_measure=p0["dewitt_measure"],
            mu2=p0["mu2"],
            lam4=p0["lam4"],
            n_warmup=p0["n_warmup"],
            n_steps=p0["n_steps"],
            Q2_free=b1_workers[0]["observables"]["Q2_free"],
        ),
        seeds=seeds,
        final_eps_list=eps_list,
        acceptance_mean=float(np.mean(acc_list)),
        sd_z_mean=float(np.mean(sd_z_list)),
        sd_phi_mean=float(np.mean(sd_phi_list)),
        total_b200_gpu_hours=total_wall_sec / 3600.0,
        diagnostics=diags,
    )


def analyze_all(base_dir: str):
    fin_path = os.path.join(base_dir, "modal_hmc_final_results.json")
    b1_path = os.path.join(base_dir, "modal_hmc_followup_batch1_results.json")
    b2_n32_path = os.path.join(base_dir, "modal_hmc_followup_b2_n32_results.json")
    b2_n48_path = os.path.join(base_dir, "modal_hmc_followup_b2_n48_results.json")
    b2_both_path = os.path.join(base_dir, "modal_hmc_followup_b2_both_results.json")

    fin = json.load(open(fin_path))
    b1 = json.load(open(b1_path))
    fin_by_tag = {r["tag"]: r for r in fin["runs"]}
    b1_runs = b1["runs"]

    b2_runs_n32, b2_runs_n48 = [], []
    b2_gpu_hours = 0.0
    for p in [b2_n32_path, b2_n48_path, b2_both_path]:
        if os.path.exists(p):
            d2 = json.load(open(p))
            b2_gpu_hours += d2.get("total_b200_gpu_hours", 0.0)
            for r in d2["runs"]:
                if "N32" in r["tag"]:
                    b2_runs_n32.append(r)
                elif "N48" in r["tag"]:
                    b2_runs_n48.append(r)

    # 1. Four matched zero-coupling controls (full chain histories)
    ctrl = {}
    for N in [8, 16, 24, 32]:
        ctrl[N] = load_combined_histories(b1_runs, f"followup_ctrl_gN0.00_N{N}")

    # 2. Primary 4-rank roughness smoothing table (measured g_N=0 control vs analytic Q2_free)
    primary_rows = []
    for N in [8, 16, 24, 32]:
        r2 = fin_by_tag[f"final_UV_gN2.00_N{N}"]
        c0 = ctrl[N]
        q2_0_mean = c0["diagnostics"]["Q2"]["mean"]
        q2_0_se_ch = c0["diagnostics"]["Q2"]["se_chains"]
        q2_0_se_ess = c0["diagnostics"]["Q2"]["se_ess"]
        # Use conservative max(se_chains, se_ess) to guard against lucky-tight spread across C=4 chains (e.g. N=32)
        q2_0_se = max(q2_0_se_ch, q2_0_se_ess)
        q2_2_mean = r2["observables"]["Q2_anisotropy"]["mean"]
        q2_2_se = r2["observables"]["Q2_anisotropy"]["se"]
        q2_free = r2["observables"]["Q2_free"]
        assert abs(q2_free - c0["params"]["Q2_free"]) < 1e-12

        df0 = c0["n_chains_total"] - 1
        df2 = r2["params"]["n_chains"] - 1

        # Primary smoothing relative to matched measured g_N = 0 (with lam4 = 0.1)
        st_meas = smoothing_stats(q2_2_mean, q2_2_se, q2_0_mean, q2_0_se, df_num=df2, df_den=df0, seed=100 + N)
        # Analytic free-field baseline (g_N = 0, lam4 = 0)
        st_free = smoothing_stats(q2_2_mean, q2_2_se, q2_free, 0.0, seed=200 + N)
        # Pure scalar quartic shift at g_N = 0 (lam4 = 0.1 vs lam4 = 0)
        st_quartic = smoothing_stats(q2_0_mean, q2_0_se, q2_free, 0.0, seed=300 + N)

        primary_rows.append(
            dict(
                N=N,
                lattice=f"{r2['params']['Nt']}x{r2['params']['Nr']}",
                n_chains=r2["params"]["n_chains"],
                n_meas=r2["params"]["n_meas"],
                seed_g0=c0["seeds"][0],
                seed_g2=r2["params"]["seed"],
                Q2_free=q2_free,
                Q2_g0_mean=q2_0_mean,
                Q2_g0_se=q2_0_se,
                Q2_g0_se_chains=q2_0_se_ch,
                Q2_g0_se_ess=q2_0_se_ess,
                Q2_g2_mean=q2_2_mean,
                Q2_g2_se_chains=q2_2_se,
                smoothing_vs_measured_g0=st_meas,
                smoothing_vs_analytic_free=st_free,
                quartic_only_shift_vs_free=st_quartic,
            )
        )

    # 3. Full g_N = 0.5, 1.0, 1.5, 2.0 cross-check against BOTH denominators
    crosscheck_grid = []
    for N in [8, 16, 24, 32]:
        c0 = ctrl[N]
        q2_0_mean = c0["diagnostics"]["Q2"]["mean"]
        q2_0_se = max(c0["diagnostics"]["Q2"]["se_chains"], c0["diagnostics"]["Q2"]["se_ess"])
        q2_free = c0["params"]["Q2_free"]
        for gN in [0.0, 0.5, 1.0, 1.5, 2.0]:
            if gN == 0.0:
                m, se, src = q2_0_mean, q2_0_se, "followup_ctrl_full_history"
            else:
                rg = fin_by_tag[f"final_UV_gN{gN:.2f}_N{N}"]
                m = rg["observables"]["Q2_anisotropy"]["mean"]
                se = rg["observables"]["Q2_anisotropy"]["se"]
                src = "archived_final_summary"
            st_m = smoothing_stats(m, se, q2_0_mean, q2_0_se if gN > 0 else 0.0, seed=400 + int(10 * gN) + N)
            st_f = smoothing_stats(m, se, q2_free, 0.0, seed=500 + int(10 * gN) + N)
            crosscheck_grid.append(
                dict(
                    N=N,
                    g_N=gN,
                    source=src,
                    Q2_mean=m,
                    Q2_se=se,
                    Q2_g0_mean=q2_0_mean,
                    Q2_free=q2_free,
                    ratio_vs_measured_g0=m / q2_0_mean,
                    ratio_vs_measured_g0_se=st_m["ratio_se"],
                    smoothing_vs_measured_g0_pct=st_m["smoothing_pct"],
                    smoothing_vs_measured_g0_se_pct=st_m["smoothing_se_pct"],
                    ratio_vs_analytic_free=m / q2_free,
                    ratio_vs_analytic_free_se=se / q2_free,
                    smoothing_vs_analytic_free_pct=st_f["smoothing_pct"],
                )
            )

    # 4. Larger-lattice 12x8 matched pair at N = 24
    lat_g0 = load_combined_histories(b1_runs, "followup_lat12x8_gN0.00_N24_w")
    lat_g2 = load_combined_histories(b1_runs, "followup_lat12x8_gN2.00_N24_w")
    lat_g0_se = max(lat_g0["diagnostics"]["Q2"]["se_chains"], lat_g0["diagnostics"]["Q2"]["se_ess"])
    lat_g2_se = max(lat_g2["diagnostics"]["Q2"]["se_chains"], lat_g2["diagnostics"]["Q2"]["se_ess"])
    st_lat12x8 = smoothing_stats(
        lat_g2["diagnostics"]["Q2"]["mean"],
        lat_g2_se,
        lat_g0["diagnostics"]["Q2"]["mean"],
        lat_g0_se,
        df_num=lat_g2["n_chains_total"] - 1,
        df_den=lat_g0["n_chains_total"] - 1,
        seed=601,
    )
    st_lat12x8_free = smoothing_stats(
        lat_g2["diagnostics"]["Q2"]["mean"],
        lat_g2_se,
        lat_g2["params"]["Q2_free"],
        0.0,
        seed=602,
    )
    row24_8x4 = next(r for r in primary_rows if r["N"] == 24)
    s_12x8 = st_lat12x8["smoothing_pct"]
    se_12x8 = st_lat12x8["smoothing_se_pct"]
    s_8x4 = row24_8x4["smoothing_vs_measured_g0"]["smoothing_pct"]
    se_8x4 = row24_8x4["smoothing_vs_measured_g0"]["smoothing_se_pct"]
    diff_pp = s_12x8 - s_8x4
    diff_se_pp = math.sqrt(se_12x8**2 + se_8x4**2)

    lattice_test = dict(
        N=24,
        note="Both lattice dimensions change simultaneously: 8x4 -> 12x8 (volume 32 -> 96).",
        lat_12x8_g0=lat_g0,
        lat_12x8_g2=lat_g2,
        smoothing_12x8_vs_measured_g0=st_lat12x8,
        smoothing_12x8_vs_analytic_free=st_lat12x8_free,
        smoothing_8x4_vs_measured_g0=row24_8x4["smoothing_vs_measured_g0"],
        diff_12x8_minus_8x4_pp=float(diff_pp),
        diff_se_pp=float(diff_se_pp),
        diff_ci95_lo_pp=float(diff_pp - 1.959964 * diff_se_pp),
        diff_ci95_hi_pp=float(diff_pp + 1.959964 * diff_se_pp),
    )

    # 5. Extended independent-chain area convergence runs (N=32 fixed-g_N=2, N=48 scaled g_P=3.3821)
    ext_n32_b1 = load_combined_histories(b1_runs, "followup_ext_UV_gN2.00_N32_b1_w")
    ext_n48_b1 = load_combined_histories(b1_runs, "followup_ext_IR_mu005_N48_b1_w")

    # Batch 2 runs resumed from stationary Batch 1 checkpoints with exact spacetime-local collective moves
    # (homogeneous stationary Markov chains of 16 chains x 2000 draws = 32,000 draws each).
    ext_n32_b2 = load_combined_histories(b2_runs_n32, "followup_ext_UV_gN2.00_N32_b2_w") if b2_runs_n32 else None
    ext_n48_b2 = load_combined_histories(b2_runs_n48, "followup_ext_IR_mu005_N48_b2_w") if b2_runs_n48 else None

    ext_n32 = ext_n32_b2 if ext_n32_b2 is not None else ext_n32_b1
    ext_n48 = ext_n48_b2 if ext_n48_b2 is not None else ext_n48_b1

    # Also cross-check N=32 roughness smoothing using the extended 16-chain g_N=2 run
    st_n32_ext = smoothing_stats(
        ext_n32["diagnostics"]["Q2"]["mean"],
        max(ext_n32["diagnostics"]["Q2"]["se_chains"], ext_n32["diagnostics"]["Q2"]["se_ess"]),
        ctrl[32]["diagnostics"]["Q2"]["mean"],
        max(ctrl[32]["diagnostics"]["Q2"]["se_chains"], ctrl[32]["diagnostics"]["Q2"]["se_ess"]),
        df_num=ext_n32["n_chains_total"] - 1,
        df_den=ctrl[32]["n_chains_total"] - 1,
        seed=701,
    )

    # Check whether area convergence targets (rhat_rank < 1.01 and bulk/tail ESS >= 400) are met
    def _area_gate(ext_obj: dict[str, Any]) -> dict[str, Any]:
        gates = {}
        all_pass = True
        for k in ["Q2", "phi0", "rho2_eig_mean", "rho2_eig_min", "rho2_eig_max"]:
            dk = ext_obj["diagnostics"][k]
            ok = (
                (dk["rhat_rank"] is not None)
                and (dk["rhat_rank"] < 1.01)
                and (dk["ess_bulk"] >= 400.0)
                and (dk["ess_tail"] >= 400.0)
            )
            gates[k] = dict(
                rhat_rank=dk["rhat_rank"],
                ess_bulk=dk["ess_bulk"],
                ess_tail=dk["ess_tail"],
                meets_gate=ok,
            )
            if not ok:
                all_pass = False
        return dict(all_observables_pass=all_pass, per_observable=gates)

    # Archived summaries (explicitly labeled as archived summaries without raw chain histories)
    archived_summaries = {}
    for tag in [
        "final_UV_gN2.00_N8",
        "final_UV_gN2.00_N16",
        "final_UV_gN2.00_N24",
        "final_UV_gN2.00_N32",
        "final_IR_mu005_N48",
    ]:
        r = fin_by_tag[tag]
        archived_summaries[tag] = dict(
            source="archived_summary_only_no_chain_histories",
            params=r["params"],
            observables={
                k: r["observables"][k]
                for k in ["Q2_anisotropy", "phi0_tadpole", "rho2_eig_mean", "rho2_eig_min", "rho2_eig_max"]
            },
        )

    summary = dict(
        v5_reference_commit="f527e17520167b04494fd137c381757ac501b5de",
        followup_base_commit=b1.get("base_commit"),
        total_followup_b200_gpu_hours=b1["total_b200_gpu_hours"] + b2_gpu_hours,
        batch1_b200_gpu_hours=b1["total_b200_gpu_hours"],
        batch2_b200_gpu_hours=b2_gpu_hours,
        primary_4rank_table=primary_rows,
        crosscheck_gN_grid=crosscheck_grid,
        lattice_12x8_test=lattice_test,
        controls_full_history={str(N): ctrl[N] for N in [8, 16, 24, 32]},
        extended_n32_batch1=ext_n32_b1,
        extended_n32_batch2=ext_n32_b2,
        extended_n32_fixed_gN2=ext_n32,
        extended_n32_area_gate=_area_gate(ext_n32),
        extended_n32_smoothing_crosscheck=st_n32_ext,
        extended_n48_batch1=ext_n48_b1,
        extended_n48_batch2=ext_n48_b2,
        extended_n48_scaled_mu005=ext_n48,
        extended_n48_area_gate=_area_gate(ext_n48),
        archived_summaries=archived_summaries,
    )

    out_path = os.path.join(base_dir, "modal_hmc_followup_analysis.json")
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)

    # Print crisp terminal report
    print("\n================ FOLLOW-UP ANALYSIS SUMMARY ================")
    print(f"Total B200 GPU-hours: {summary['total_followup_b200_gpu_hours']:.3f} (B1={summary['batch1_b200_gpu_hours']:.3f}, B2={summary['batch2_b200_gpu_hours']:.3f})")
    print("\n--- 1. Primary 4-Rank Horizon Roughness Smoothing (g_N=2 vs matched measured g_N=0, lam4=0.1) ---")
    for row in primary_rows:
        sm = row["smoothing_vs_measured_g0"]
        sf = row["smoothing_vs_analytic_free"]
        sq = row["quartic_only_shift_vs_free"]
        print(
            f"N={row['N']:2d} ({row['lattice']:4s}) | "
            f"Q2(0)={row['Q2_g0_mean']:.6f}±{row['Q2_g0_se_chains']:.6f} | "
            f"Q2(2)={row['Q2_g2_mean']:.6f}±{row['Q2_g2_se_chains']:.6f} | "
            f"Q2_free={row['Q2_free']:.6f} | "
            f"Smooth(meas g=0)={sm['smoothing_pct']:.3f}% [95% MC: {sm['ci95_mc_lo']:.3f}%, {sm['ci95_mc_hi']:.3f}%] | "
            f"Smooth(free)={sf['smoothing_pct']:.3f}% | QuarticOnly={sq['smoothing_pct']:.3f}%"
        )

    print("\n--- 2. Larger-Lattice Test at N=24 (12x8 vs 8x4, both dims changed) ---")
    print(
        f"12x8: Q2(0)={lat_g0['diagnostics']['Q2']['mean']:.6f}±{lat_g0['diagnostics']['Q2']['se_chains']:.6f}, "
        f"Q2(2)={lat_g2['diagnostics']['Q2']['mean']:.6f}±{lat_g2['diagnostics']['Q2']['se_chains']:.6f}, "
        f"Smooth(12x8)={s_12x8:.3f}% ± {se_12x8:.3f}% [95% MC: {st_lat12x8['ci95_mc_lo']:.3f}%, {st_lat12x8['ci95_mc_hi']:.3f}%]"
    )
    print(
        f"8x4:  Smooth(8x4) ={s_8x4:.3f}% ± {se_8x4:.3f}% | "
        f"Diff(12x8 - 8x4) = {diff_pp:+.3f} ± {diff_se_pp:.3f} pp "
        f"[95% CI: {lattice_test['diff_ci95_lo_pp']:+.3f}, {lattice_test['diff_ci95_hi_pp']:+.3f} pp]"
    )

    print("\n--- 3. Convergence Diagnostics (Vehtari 2021 rank-split R-hat, bulk-ESS, tail-ESS) ---")
    diag_list = [
        ("ctrl_g0_N8", ctrl[8]),
        ("ctrl_g0_N16", ctrl[16]),
        ("ctrl_g0_N24", ctrl[24]),
        ("ctrl_g0_N32", ctrl[32]),
        ("lat12x8_g0_N24", lat_g0),
        ("lat12x8_g2_N24", lat_g2),
        ("ext_UV_g2_N32_B1", ext_n32_b1),
        ("ext_IR_mu005_N48_B1", ext_n48_b1),
    ]
    if ext_n32_b2 is not None:
        diag_list.append(("ext_UV_g2_N32_B2", ext_n32_b2))
    if ext_n48_b2 is not None:
        diag_list.append(("ext_IR_mu005_N48_B2", ext_n48_b2))
    for label, obj in diag_list:
        print(f"  [{label}] chains={obj['n_chains_total']} x draws={obj['n_draws_per_chain']} (N_tot={obj['n_samples_total']}, SD_Z={obj['sd_z_mean']:.5f}, SD_Phi={obj['sd_phi_mean']:.5f})")
        for k in ["Q2", "phi0", "rho2_eig_mean", "rho2_eig_min", "rho2_eig_max"]:
            dk = obj["diagnostics"][k]
            if dk["is_constant"]:
                print(f"    {k:14s}: constant={dk['mean']:.4f} (N/A at g_N=0)")
            else:
                gate = "PASS" if (dk["rhat_rank"] < 1.01 and dk["ess_bulk"] >= 400 and dk["ess_tail"] >= 400) else "CHECK"
                print(
                    f"    {k:14s}: mean={dk['mean']:.6f}±{dk['se_chains']:.6f} | "
                    f"rhat_rank={dk['rhat_rank']:.4f} (bulk={dk['rhat_rank_bulk']:.4f}, fold={dk['rhat_rank_folded']:.4f}) | "
                    f"ESS_bulk={dk['ess_bulk']:.1f}, ESS_tail={dk['ess_tail']:.1f} [{gate}]"
                )

    print(f"\nArea convergence gates: N=32 -> {summary['extended_n32_area_gate']['all_observables_pass']}, N=48 -> {summary['extended_n48_area_gate']['all_observables_pass']}")
    print(f"Extended N=32 smoothing cross-check (vs measured g_N=0): {st_n32_ext['smoothing_pct']:.3f}% [95% MC: {st_n32_ext['ci95_mc_lo']:.3f}%, {st_n32_ext['ci95_mc_hi']:.3f}%]")
    print(f"Saved full analysis to {out_path}")


if __name__ == "__main__":
    unit_test_diagnostics()
    if "--unit-only" not in sys.argv:
        here = os.path.dirname(os.path.abspath(__file__))
        analyze_all(here)
