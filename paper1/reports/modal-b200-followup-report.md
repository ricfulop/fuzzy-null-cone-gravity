# Targeted B200 Follow-Up: Matched Zero-Coupling Roughness Controls, Larger-Lattice Test, and Separate Local-Area Convergence

**Repository:** `ricfulop/fuzzy-null-cone-gravity`  
**v5 Reference Commit:** `f527e17520167b04494fd137c381757ac501b5de` (bit-for-bit identical `modal_fsg_hmc.py` and `modal_hmc_final_results.json` on `main` at `589a99b9a5fa5695ebaf20cd3082b76f94a02346`)  
**Hardware:** Modal NVIDIA B200 GPUs (FP64)  
**Actual Compute Used:** `1.774 B200 GPU-hours` total (`1.054 GPU-hrs` Batch 1 across 16 GPUs + `0.721 GPU-hrs` Batch 2 across 8 GPUs; ~\$11.09 USD at \$6.25/GPU-hr)

---

## 1. Executive Summary

Holding the quartic scalar self-interaction ($\lambda_4 = 0.1$), area stiffness ($\mu^2 = 0.25$), lattice geometry, and sampling procedure strictly fixed, turning on the spin-2 graviton coupling from $g_N = 0$ to $g_N = 2.0$ suppresses traceless horizon roughness $Q_2 = \frac{1}{V}\sum_x \frac{1}{N}\mathrm{Tr}_N(\Phi_{\mathrm{trless}}^2)$ by **$47.68\%$–$48.02\%$** across $N \in \{8, 16, 24, 32\}$. The quartic scalar self-interaction alone accounts for only **$0.010\%$–$0.166\%$** of the reduction relative to the analytic free-field baseline $Q_{2,\mathrm{free}}$ ($\lambda_4 = 0$), confirming that the earlier $\approx 47.8\%$ estimate was not a quartic-baseline artifact (departures from the analytic-free denominator are at most $0.087$ percentage points at $N=8$ and $\le 0.005$ pp for $N \ge 24$).

On the larger **$12 \times 8$ lattice at $N = 24$** (where both lattice dimensions increase simultaneously from $8 \times 4$, tripling spacetime volume from $32$ to $96$), the measured smoothing is **$47.621\%$** ($95\%$ MC interval $[47.572\%, 47.670\%]$), differing from the $8 \times 4$ result ($47.681\%$) by **$-0.060 \pm 0.047$ percentage points** ($95\%$ CI $[-0.151, +0.032]$ pp). The suppression persists unchanged on the larger lattice.

Evaluating convergence separately for roughness ($Q_2$, global $\phi_0$) and local area ($\widehat{\rho}^2(0,0)$ mean and extreme eigenvalues) via Vehtari et al. (2021) rank-normalized split-$\hat{R}$ and bulk/tail ESS shows that all roughness ensembles converge cleanly ($\hat{R}_{\mathrm{rank}} \le 1.0099$, bulk/tail $\mathrm{ESS} \ge 516$). For local area, Batch 1 revealed genuine $k_\parallel \neq 0$ local-monopole autocorrelation ($\hat{R}_{\mathrm{rank}} = 1.108$ at $N=32$ and $1.033$ at $N=48$). Adding an exact spacetime-local Fourier-wave collective scale move ($\Phi(x) \mapsto \Phi(x) + \frac{\delta(x)}{\sqrt{N}} I_N$, $Z(x) \mapsto e^{-\frac{g_N}{2\sqrt{N}}\delta(x)} Z(x)$, leaving $S_{z,\mathrm{ang}}$ and $-\Delta S_{\mathrm{DW}} + \ln J$ identically zero to $< 3\times 10^{-13}$) in Batch 2 converged **both** extended 16-chain ensembles decisively past the target gates ($\hat{R}_{\mathrm{rank}} < 1.01$ and bulk/tail $\mathrm{ESS} \ge 400$):
- **$N=32$, fixed $g_N=2.0$, $\mu^2=0.25$:** local mean-area $\hat{R}_{\mathrm{rank}}$ dropped from `1.416` (archived) $\to$ `1.108` (B1) $\to$ **`1.0098`** (B2), with bulk $\mathrm{ESS} = 2{,}576$ and tail $\mathrm{ESS} = 5{,}624$. The converged local area mean shifts modestly downward from the unconverged archived value ($1.331 \pm 0.111$, $\lambda \in [0.374, 3.353]$) to **$1.2107 \pm 0.0057$** (**$\lambda \in [0.342, 3.038]$**), remaining strictly positive on all $32{,}000$ draws.
- **$N=48$, scaled $g_P=3.3821$, $\mu^2=0.05$:** local mean-area $\hat{R}_{\mathrm{rank}}$ dropped from `1.181` (archived) $\to$ `1.033` (B1) $\to$ **`1.0021`** (B2), with bulk $\mathrm{ESS} = 4{,}020$ and tail $\mathrm{ESS} = 9{,}123$, yielding **$\langle\widehat{\rho}^2\rangle = 1.0142 \pm 0.0009$** and **$\lambda \in [0.726, 1.381]$**.

---

## 2. Primary Four-Rank Comparison Table (Explicit Denominators)

All ensembles use `fixed_gN=True`, `dewitt_measure=True`, $\mu^2 = 0.25$, $\lambda_4 = 0.1$, `n_warmup=220`, `n_steps=10`, `init_eps=0.10`. Controls use fresh independent seeds (`9101–9104`) and exact matched `n_chains` / `n_meas` (`final_UV_gN2.00_N{N}` seeds `8675, 8742, 8809, 8876`). Conservative SE $\max(\mathrm{SE}_{\mathrm{chains}}, \mathrm{SE}_{\mathrm{ESS}})$ is used for full-history controls so the 4-chain $N=32$ control is not artificially narrowed by small-$C$ sampling luck.

| Rank $N$ (Lattice) | $\overline{Q}_2(g_N=0)$ (measured, $\lambda_4=0.1$) | $\overline{Q}_2(g_N=2)$ (archived summary) | $Q_{2,\mathrm{free}}$ (analytic, $\lambda_4=0$) | Primary Smoothing $100[1 - \overline{Q}_2(2)/\overline{Q}_2(0)]$ | 95% MC Interval (both ensembles) | Analytic-Free Smoothing $100[1 - \overline{Q}_2(2)/Q_{2,\mathrm{free}}]$ | Quartic-Only Shift $100[1 - \overline{Q}_2(0)/Q_{2,\mathrm{free}}]$ |
|---|---:|---:|---:|---:|---|---:|---:|
| **$N=8$** ($8\times 4$, 8 chains) | $0.339394 \pm 0.000370$ | $0.176430 \pm 0.000066$ | $0.339959$ | **$48.016\%$** | $[47.899\%, 48.134\%]$ | $48.103 \pm 0.019\%$ | $+0.166 \pm 0.109\%$ |
| **$N=16$** ($8\times 4$, 8 chains) | $0.253634 \pm 0.000140$ | $0.132504 \pm 0.000080$ | $0.253780$ | **$47.758\%$** | $[47.674\%, 47.841\%]$ | $47.788 \pm 0.032\%$ | $+0.058 \pm 0.055\%$ |
| **$N=24$** ($8\times 4$, 4 chains) | $0.202589 \pm 0.000144$ | $0.105993 \pm 0.000027$ | $0.202612$ | **$47.681\%$** | $[47.603\%, 47.758\%]$ | $47.686 \pm 0.013\%$ | $+0.011 \pm 0.071\%$ |
| **$N=32$** ($12\times 4$, 4 chains) | $0.169826 \pm 0.000141^\dagger$ | $0.088697 \pm 0.000037$ | $0.169842$ | **$47.772\%$** | $[47.676\%, 47.866\%]$ | $47.777 \pm 0.022\%$ | $+0.010 \pm 0.083\%$ |

$^\dagger$ At $N=32$ ($C=4$), $\mathrm{SE}_{\mathrm{chains}} = 0.000018$ while $\mathrm{SE}_{\mathrm{ESS}} = 0.000141$; we propagate the conservative $\mathrm{SE}_{\mathrm{ESS}}$ above (using $\mathrm{SE}_{\mathrm{chains}}$ alone would yield $[47.727\%, 47.816\%]$). Cross-checking $N=32$ with the extended 16-chain Batch 2 $g_N=2$ run ($\overline{Q}_2 = 0.088668 \pm 0.000012$) gives **$47.789\%$** ($95\%$ MC $[47.702\%, 47.875\%]$).

### Cross-Check Across $g_N \in \{0.0, 0.5, 1.0, 1.5, 2.0\}$ (Both Denominators)

| $g_N$ | $N=8$ ($\overline{Q}_2(0)$ / $Q_{2,\mathrm{free}}$) | $N=16$ ($\overline{Q}_2(0)$ / $Q_{2,\mathrm{free}}$) | $N=24$ ($\overline{Q}_2(0)$ / $Q_{2,\mathrm{free}}$) | $N=32$ ($\overline{Q}_2(0)$ / $Q_{2,\mathrm{free}}$) |
|---:|---|---|---|---|
| `0.0` | $0.00 \pm 0.11\%$ / $0.17\%$ | $0.00 \pm 0.06\%$ / $0.06\%$ | $0.00 \pm 0.07\%$ / $0.01\%$ | $0.00 \pm 0.08\%$ / $0.01\%$ |
| `0.5` | $8.05 \pm 0.12\%$ / $8.20\%$ | $9.75 \pm 0.08\%$ / $9.80\%$ | $11.08 \pm 0.10\%$ / $11.09\%$ | $12.30 \pm 0.09\%$ / $12.31\%$ |
| `1.0` | $23.24 \pm 0.10\%$ / $23.37\%$ | $25.23 \pm 0.05\%$ / $25.27\%$ | $26.57 \pm 0.06\%$ / $26.58\%$ | $27.79 \pm 0.06\%$ / $27.80\%$ |
| `1.5` | $37.00 \pm 0.08\%$ / $37.10\%$ | $37.93 \pm 0.04\%$ / $37.96\%$ | $38.64 \pm 0.06\%$ / $38.64\%$ | $39.30 \pm 0.05\%$ / $39.31\%$ |
| `2.0` | $48.02 \pm 0.06\%$ / $48.10\%$ | $47.76 \pm 0.04\%$ / $47.79\%$ | $47.68 \pm 0.04\%$ / $47.69\%$ | $47.77 \pm 0.05\%$ / $47.78\%$ |

---

## 3. Larger-Lattice Test at $N = 24$ ($12 \times 8$ vs. $8 \times 4$)

Both lattice dimensions change simultaneously ($N_\tau: 8 \to 12$, $N_r: 4 \to 8$; volume $32 \to 96$), holding $N=24$, $\mu^2=0.25$, $\lambda_4=0.1$, `fixed_gN=True`, `dewitt_measure=True` fixed (8 independent chains × 700 post-warmup draws per coupling; seeds `9201–9202` at $g_N=0$, `9211–9212` at $g_N=2$):

- **$12 \times 8$ measured means:** $\overline{Q}_2(g_N=0) = 0.202480 \pm 0.000080$, $\overline{Q}_2(g_N=2) = 0.106057 \pm 0.000028$ ($Q_{2,\mathrm{free}}^{12\times 8} = 0.202454$).
- **$12 \times 8$ smoothing:** **$47.621\% \pm 0.025\%$** ($95\%$ MC interval $[47.572\%, 47.670\%]$).
- **$8 \times 4$ smoothing at $N=24$:** **$47.681\% \pm 0.039\%$** ($95\%$ MC interval $[47.603\%, 47.758\%]$).
- **Difference ($12\times 8$ minus $8\times 4$):** **$-0.060 \pm 0.047$ percentage points** ($95\%$ CI $[-0.151, +0.032]$ pp).

---

## 4. Separate Convergence Diagnostics: Roughness vs. Local Area

### 4A. New Full-History Ensembles (Vehtari et al. 2021 Rank-Split $\hat{R}$, Bulk-ESS, Tail-ESS)

| Ensemble (Chains × Draws, Seeds) | $\langle\mathrm{SD}_Z\rangle$ / $\langle\mathrm{SD}_\Phi\rangle$ | $Q_2$: $\hat{R}_{\mathrm{rank}}$ / Bulk / Tail ESS | $\phi_0$: $\hat{R}_{\mathrm{rank}}$ / Bulk / Tail ESS | $\langle\widehat{\rho}^2\rangle$: Mean ± SE | $\langle\widehat{\rho}^2\rangle$: $\hat{R}_{\mathrm{rank}}$ / Bulk / Tail ESS | $\lambda_{\min}$ / $\lambda_{\max}$: $\hat{R}_{\mathrm{rank}}$ (min ESS) |
|---|---|---|---|---:|---|---|
| `ctrl_g0_N8` ($8\times 700$, `9101`) | 0.9992 / 1.0000 | `1.0016` / 3274 / 2810 | `1.0000` / 4514 / 4887 | `1.0` (const) | N/A ($g_N=0$) | N/A ($g_N=0$) |
| `ctrl_g0_N16` ($8\times 700$, `9102`) | 1.0001 / 1.0005 | `1.0019` / 2919 / 2895 | `1.0010` / 4062 / 4759 | `1.0` (const) | N/A ($g_N=0$) | N/A ($g_N=0$) |
| `ctrl_g0_N24` ($4\times 700$, `9103`) | 0.9999 / 1.0006 | `1.0048` / 1197 / 1075 | `1.0039` / 1631 / 2174 | `1.0` (const) | N/A ($g_N=0$) | N/A ($g_N=0$) |
| `ctrl_g0_N32` ($4\times 600$, `9104`) | 0.9998 / 1.0002 | `1.0099` / 516 / 676 | `1.0027` / 1402 / 1817 | `1.0` (const) | N/A ($g_N=0$) | N/A ($g_N=0$) |
| `lat12x8_g0_N24` ($8\times 700$, `9201–9202`) | 0.9997 / 1.0001 | `1.0040` / 1355 / 2489 | `1.0040` / 3279 / 3687 | `1.0` (const) | N/A ($g_N=0$) | N/A ($g_N=0$) |
| `lat12x8_g2_N24` ($8\times 700$, `9211–9212`) | 1.0000 / 1.0004 | `1.0041` / 896 / 2087 | `1.0011` / 2886 / 3658 | $1.287 \pm 0.059$ | `1.3268` / 20 / 42 *(unconverged)* | `1.262` / `1.274` (22) |
| `ext_UV_g2_N32_B1` ($16\times 2000$, `9301–9304`, global only) | 0.9999 / 1.0000 | `1.0040` / 4588 / 8394 | `1.0010` / 15505 / 19379 | $1.187 \pm 0.018$ | `1.1082` / 96 / 372 *(below gate)* | `1.094` / `1.100` (103) |
| `ext_IR_mu005_N48_B1` ($16\times 2000$, `9401–9404`, global only) | 1.0000 / 1.0001 | `1.0027` / 5447 / 10128 | `1.0007` / 16230 / 20692 | $1.019 \pm 0.003$ | `1.0325` / 369 / 1661 *(below gate)* | `1.030` / `1.027` (396) |
| **`ext_UV_g2_N32_B2`** ($16\times 2000$, `9501–9504`, +local coll.) | **1.0000 / 1.0000** | **`1.0032` / 4575 / 9053** | **`1.0006` / 14560 / 19282** | **$1.2107 \pm 0.0057$** | **`1.0098` / 2576 / 5624 [PASS]** | **`1.0083` / `1.0091` (2738) [PASS]** |
| **`ext_IR_mu005_N48_B2`** ($16\times 2000$, `9601–9604`, +local coll.) | **1.0000 / 1.0000** | **`1.0021` / 5158 / 9816** | **`1.0005` / 16624 / 20513** | **$1.0142 \pm 0.0009$** | **`1.0021` / 4020 / 9123 [PASS]** | **`1.0022` / `1.0024` (4371) [PASS]** |

### 4B. Archived Summaries (Original Chain Histories Unavailable — Classic Unsplit $\hat{R}$ Only)

| Archived Tag (Seed) | $Q_2$: classic $\hat{R}$ ($\tau_{\mathrm{int}}$) | $\phi_0$: classic $\hat{R}$ | $\langle\widehat{\rho}^2\rangle$: Mean ± SE | $\langle\widehat{\rho}^2\rangle$: classic $\hat{R}$ ($\tau_{\mathrm{int}}$) |
|---|---|---|---:|---|
| `final_UV_gN2.00_N8` (`8675`) | `1.0000` (1.02) | `1.0009` | $1.5047 \pm 0.0308$ | `1.0118` (3.62) |
| `final_UV_gN2.00_N16` (`8742`) | `1.0032` (1.44) | `1.0002` | $1.3579 \pm 0.0421$ | `1.0503` (16.89) |
| `final_UV_gN2.00_N24` (`8809`) | `1.0000` (2.06) | `1.0000` | $1.2856 \pm 0.0488$ | `1.0685` (18.21) |
| `final_UV_gN2.00_N32` (`8876`) | `1.0034` (4.51) | `1.0014` | $1.3305 \pm 0.1111$ | `1.4163` (35.07) — *resolved by `ext_UV_g2_N32_B2`* |
| `final_IR_mu005_N48` (`5492`) | `1.0211` (2.98) | `1.0050` | $1.0206 \pm 0.0191$ | `1.1807` (11.85) — *resolved by `ext_IR_mu005_N48_B2`* |

---

## 5. Actual B200 Compute Accounting

| Group | GPUs | Wall-Sum (s) | B200 GPU-Hours | Output File |
|---|---:|---:|---:|---|
| Matched $g_N=0$ controls ($N=8,16,24,32$) | 4 | 188.5 | 0.0524 | `modal_hmc_followup_batch1_results.json` |
| $12\times 8$ matched pair ($N=24$, $g_N=0,2$) | 4 | 460.8 | 0.1280 | `modal_hmc_followup_batch1_results.json` |
| Extended $N=32$ Batch 1 (16 chains, global coll.) | 4 | 1376.5 | 0.3824 | `modal_hmc_followup_batch1_results.json` |
| Extended $N=48$ Batch 1 (16 chains, global coll.) | 4 | 1767.0 | 0.4908 | `modal_hmc_followup_batch1_results.json` |
| Extended $N=32$ Batch 2 (16 chains, +local coll.) | 4 | 1020.8 | 0.2836 | `modal_hmc_followup_b2_both_results.json` |
| Extended $N=48$ Batch 2 (16 chains, +local coll.) | 4 | 1574.2 | 0.4373 | `modal_hmc_followup_b2_both_results.json` |
| **Total Follow-Up** | **24** | **6387.8** | **1.7744** | **`modal_hmc_followup_analysis.json`** |
