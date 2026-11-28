"""Independent algebra controls for the original H6 recursion.

Inert on import. The tracked driver calls run_unreduced_algebra_controls.
Analytic targets use manufactured TT fields, exact noncommuting matrices,
generic rational polynomial multiplication of the original Lichnerowicz
equation, isolated energy impulses, and oriented surface fluxes. These
pointwise controls do not solve any constraint or validate an integrated
finite-ball identity, an exterior completion, or a gravitational H6 value.
"""

from collections.abc import Mapping
from fractions import Fraction as F
from math import fsum, isfinite

from sixth_unreduced_algebra import (
    SPATIAL_FAMILIES, a2_flux, assemble_integrated_energy, energy_densities,
    fourth_sources, scalar_adjoint_flux, scalar_metric_operator, third_sources,
    york_norm_flux,
)


def _zero(rank):
    return 0 if rank == 0 else tuple(_zero(rank-1) for _ in range(3))


def _diag(a, b, c):
    return ((a, 0, 0), (0, b, 0), (0, 0, c))


def _jet(value=0, gradient=None, hessian=None):
    return {"value": value, "gradient": _zero(1) if gradient is None else gradient,
            "hessian": _zero(2) if hessian is None else hessian}


def _third_inputs():
    return dict(q=_diag(1, -1, 0), s=((0, 1, 0), (1, 0, 0), (0, 0, 0)),
                t0=_zero(2), r2=0, r3=0, b2=_zero(1), a1w1=_zero(1),
                u_jet=_jet(), inverse1=_diag(-1, 1, 0), dinverse1=_zero(3))


def _fourth_inputs():
    return dict(q=_diag(1, -1, 0), s=((0, 1, 0), (1, 0, 0), (0, 0, 0)),
                t0=_zero(2), t1=_zero(2), l0w2=_zero(2), ec=_zero(2),
                r2=0, r4=0, u_jet=_jet(), u3_jet=_jet(),
                inverse1=_diag(-1, 1, 0), dinverse1=_zero(3),
                inverse2=_diag(F(1, 2), F(1, 2), 0), dinverse2=_zero(3))


def _energy_inputs():
    return dict(q=_zero(2), s=_zero(2),
                sources={"r": _zero(1), "B2": 0, "B3": 0, "B4": 0, "s2": 0},
                w1=_zero(1), w2=_zero(1), b3=_zero(1), a2w1=_zero(1),
                l6=0, u=0, u3=0, u4=0)


def _product(left, right, degree=4):
    """Ordinary rational Taylor coefficients, independent of source formulas."""
    result = [F(0) for _ in range(degree+1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            if i+j <= degree:
                result[i+j] += a*b
    return result


def _negative_power(perturbation, power, degree=4):
    result = [F(1)]+[F(0)]*degree
    monomial = list(result)
    binomial = F(1)
    for k in range(1, degree+1):
        monomial = _product(monomial, perturbation, degree)
        binomial *= F(power-k+1, k)
        result = [a+binomial*b for a, b in zip(result, monomial)]
    return result


def _lichnerowicz_residual(*, u, u3, u4, lap_u, lap_u3, lap_u4,
                           d1u, d1u3, d2u, r2, r3, r4, s2, s3, s4):
    """Expand -8 Delta_bar psi + R psi - sigma^2 psi^-7 to order four."""
    perturbation = [F(0), F(0), F(u), F(u3), F(u4)]
    psi = [F(1), F(0), F(u), F(u3), F(u4)]
    curvature = [F(0), F(0), F(r2), F(r3), F(r4)]
    sigma_squared = [F(0), F(0), F(s2), F(s3), F(s4)]
    laplacian = [F(0), F(0), F(lap_u), F(lap_u3)+F(d1u),
                 F(lap_u4)+F(d1u3)+F(d2u)]
    rpsi = _product(curvature, psi)
    spower = _product(sigma_squared, _negative_power(perturbation, -7))
    return tuple(-8*a+b-c for a, b, c in zip(laplacian, rpsi, spower))


def run_unreduced_algebra_controls(config, progress=None):
    """Return complete JSON rows; any failed target fails the aggregate.

    Completed rows are emitted before the next control. Unexpected exceptions
    retain a named failed row via the callback and are then re-raised. No
    backend inverse, Fourier, HEAT, or gravitational geometry is constructed.
    """
    rows = []
    current = "settings"
    try:
        if not isinstance(config, Mapping):
            raise TypeError("config must be a mapping")
        raw_tolerance = config.get("unreduced_algebra_scaled_tolerance", 2e-12)
        if isinstance(raw_tolerance, bool):
            raise ValueError("unreduced_algebra_scaled_tolerance must be finite and positive")
        tolerance = float(raw_tolerance)
        if not isfinite(tolerance) or tolerance <= 0:
            raise ValueError("unreduced_algebra_scaled_tolerance must be finite and positive")

        def record(name, actual, reference, **metadata):
            nonlocal current
            current = name
            if set(actual) != set(reference) or not actual:
                raise ValueError("Control actual/reference keys must agree and be nonempty")
            comparisons = {}
            for key, target in reference.items():
                value, target_float = float(actual[key]), float(target)
                if not isfinite(value) or not isfinite(target_float):
                    raise ValueError("Nonfinite algebra control value")
                residual = value-target_float
                scaled = abs(residual)/max(1., abs(target_float))
                comparisons[key] = {"actual": value, "reference": target_float,
                                    "reference_exact": str(target), "residual": residual,
                                    "scaled_error": scaled, "passed": scaled <= tolerance}
            row = {"name": name, "checks": comparisons,
                   "maximum": max(item["scaled_error"] for item in comparisons.values()),
                   "passed": all(item["passed"] for item in comparisons.values()), **metadata}
            rows.append(row)
            if progress is not None:
                progress({"stage": "unreduced_algebra_control_complete", "control": name,
                          "passed": row["passed"], "control_result": row})

        # Q=[[y,x,0],[x,-y,0],[0,0,0]] is trace free and divergence free.
        # At (1,2,0), div(Q^2/2)=(1,2,0), which must not be dropped.
        inverse1 = ((-2, -1, 0), (-1, 2, 0), (0, 0, 0))
        dinverse1 = (((0, -1, 0), (-1, 0, 0), (0, 0, 0)),
                     ((-1, 0, 0), (0, 1, 0), (0, 0, 0)), _zero(2))
        inverse2 = _diag(F(5, 2), F(5, 2), 0)
        dinverse2 = (_diag(1, 1, 0), _diag(2, 2, 0), _zero(2))
        for name, jet, target in (
            ("metric_operator_linear", _jet(5, (1, 2, 0)), {"D1": 0, "D2": 5}),
            ("metric_operator_quadratic", _jet(15, (4, 13, 0), ((2, 1, 0), (1, 6, 0), (0, 0, 0))),
             {"D1": 6, "D2": 50}),
        ):
            current = name
            record(name, {"D1": scalar_metric_operator(inverse1, dinverse1, jet),
                          "D2": scalar_metric_operator(inverse2, dinverse2, jet)}, target)

        current = "noncommuting_raw_tangent"
        inputs3 = _third_inputs()
        inputs3.update(b2=(1, 2, 3), a1w1=(4, -1, 2))
        third = third_sources(**inputs3)
        fourth = fourth_sources(**_fourth_inputs())
        u1 = third["U1"]
        trace_square = fsum(u1[i][j]*u1[j][i] for i in range(3) for j in range(3))
        inputs = _energy_inputs()
        inputs.update(q=_diag(1, -1, 0), s=_third_inputs()["s"])
        density = energy_densities(**inputs)
        record(current, {"s2": third["s2"], "s3": third["s3"], "s4": fourth["s4"],
                         "trace_U1_squared": trace_square, "raw_gram": density["raw_gram"],
                         **{"r_"+str(i): value for i, value in enumerate(third["r"])}},
               {"s2": 2, "s3": 0, "s4": F(2, 3), "trace_U1_squared": -2,
                "raw_gram": -F(4, 45), "r_0": -3, "r_1": 3, "r_2": 1})

        current = "canonical_correction_sign"
        inputs = _fourth_inputs()
        k = F(3, 7)
        inputs["ec"] = tuple(tuple(k*v for v in row) for row in inputs["s"])
        changed = fourth_sources(**inputs)
        record(current, {"s4": changed["s4"], "change": changed["s4"]-fourth["s4"]},
               {"s4": F(2, 3)-4*k, "change": -4*k})

        current = "longitudinal_tangent_sign"
        inputs3, inputs4 = _third_inputs(), _fourth_inputs()
        inputs3["t0"], inputs4["t0"] = inputs3["s"], inputs4["s"]
        response_variations = {}
        for field, coefficient in (("t1", F(1, 5)), ("l0w2", F(1, 7))):
            varied = _fourth_inputs()
            varied[field] = tuple(tuple(coefficient*value for value in row) for row in varied["s"])
            response_variations[field+"_only_s4"] = fourth_sources(**varied)["s4"]
        record(current, {"s3": third_sources(**inputs3)["s3"], "s4": fourth_sources(**inputs4)["s4"],
                         **response_variations},
               {"s3": -4, "s4": F(8, 3), "t1_only_s4": -F(2, 15), "l0w2_only_s4": F(2, 21)})

        current = "original_lichnerowicz_series"
        inputs3, inputs4 = _third_inputs(), _fourth_inputs()
        inputs3["t0"], inputs4["t0"] = inputs3["s"], inputs4["s"]
        u, u3 = F(2, 7), -F(3, 8)
        uj = _jet(u, (1, -2, 3), ((2, 1, 0), (1, 6, 0), (0, 0, -1)))
        u3j = _jet(u3, (-3, 4, 1), ((-3, 2, 0), (2, 5, 0), (0, 0, 7)))
        inputs3.update(r2=F(3, 5), r3=-F(7, 11), u_jet=uj)
        inputs4.update(r2=F(3, 5), r4=F(5, 13), u_jet=uj, u3_jet=u3j)
        source3, source4 = third_sources(**inputs3), fourth_sources(**inputs4)
        residual_inputs = dict(u=u, u3=u3, u4=F(5, 9), lap_u=7, lap_u3=0, lap_u4=0,
                               d1u=4, d1u3=8, d2u=4, r2=F(3, 5), r3=-F(7, 11), r4=F(5, 13),
                               s2=2, s3=-4, s4=F(8, 3))
        residual = _lichnerowicz_residual(**residual_inputs)
        # Nonzero trial Laplacians also test the inverse's -8 coefficient.
        residual_inputs.update(lap_u3=F(2, 9), lap_u4=-F(1, 6))
        trial = _lichnerowicz_residual(**residual_inputs)
        record(current, {"F3": source3["F3"], "F4": source4["F4"],
                         "B2": source3["B2"], "B3": source3["B3"], "B4": source4["B4"],
                         "cubic_trial_residual": -8*float(F(2, 9))-source3["F3"],
                         "quartic_trial_residual": 8*float(F(1, 6))-source4["F4"]},
               {"F3": -residual[3], "F4": -residual[4],
                "B2": F(73, 5), "B3": -F(315, 11), "B4": F(743, 39),
                "cubic_trial_residual": trial[3], "quartic_trial_residual": trial[4]},
               reference_method="Exact rational convolution of original PDE; no inverse or source formula as target")

        impulses = (
            ("raw_gram", -F(4, 45)), ("york_second_response", -2),
            ("york_third_source", -4), ("york_A2_strong", 2), ("curvature", -1),
            ("scalar_B4_u", -1), ("scalar_B3_u3", -1), ("scalar_B2_u4", -1),
            ("scalar_s2_u2", 28),
        )
        for family, target in impulses:
            current = "energy_impulse_"+family
            inputs = _energy_inputs()
            if family == "raw_gram":
                inputs.update(q=_diag(1, -1, 0), s=_third_inputs()["s"])
            elif family == "york_second_response":
                inputs["sources"]["r"] = (1, 0, 0)
                inputs["w2"] = (1, 0, 0)
            elif family == "york_third_source":
                inputs.update(b3=(1, 0, 0), w1=(1, 0, 0))
            elif family == "york_A2_strong":
                inputs.update(a2w1=(1, 0, 0), w1=(1, 0, 0))
            elif family == "curvature":
                inputs["l6"] = 1
            elif family == "scalar_B4_u":
                inputs["sources"]["B4"], inputs["u"] = 1, 1
            elif family == "scalar_B3_u3":
                inputs["sources"]["B3"], inputs["u3"] = 1, 1
            elif family == "scalar_B2_u4":
                inputs["sources"]["B2"], inputs["u4"] = 1, 1
            else:
                inputs["sources"]["s2"], inputs["u"] = 1, 1
            record(current, energy_densities(**inputs),
                   {name: target if name == family else 0 for name in SPATIAL_FAMILIES},
                   scope="Independent ledger-input impulse; not a physical constraint solution")

        current = "global_gram_square_once"
        assembled = assemble_integrated_energy({name: 1 for name in SPATIAL_FAMILIES}, 5)
        record(current, {"total": assembled["total"], "ledger_entries": len(assembled["ledger"]),
                         "global_count": assembled["global_gram_square_count"],
                         "global_term": assembled["ledger"]["gram_square"]},
               {"total": 14, "ledger_entries": 10, "global_count": 1, "global_term": 5})

        current = "oriented_boundary_fluxes"
        shear = ((0, 1, 0), (1, 0, 0), (0, 0, 0))
        record(current, {
            "A2": a2_flux(q=_diag(2, -1, -1), w1=(1, 0, 0), t0=shear, t1=shear, t2=shear, normal=(0, 1, 0)),
            "York_norm": york_norm_flux((1, 0, 0), shear, (0, 1, 0)),
            "scalar_adjoint": scalar_adjoint_flux(_jet(2, (-2, 1, 0)), _jet(3, (1, -1, 2)), (1, 0, 0)),
        }, {"A2": 5, "York_norm": 2, "scalar_adjoint": 64},
               integrated_flux_identity_checked=False)

        current = "ledger_and_finite_input_guards"
        rejections = {}
        for name, supplied in (
            ("missing_family", {key: 0 for key in SPATIAL_FAMILIES[1:]}),
            ("duplicate_global_key", {**{key: 0 for key in SPATIAL_FAMILIES}, "gram_square": 1}),
            ("unreduced_canonical_cross", {**{key: 0 for key in SPATIAL_FAMILIES}, "canonical_cross": 1}),
        ):
            try:
                assemble_integrated_energy(supplied, 1)
            except ValueError:
                rejections[name] = 1
            else:
                rejections[name] = 0
        try:
            scalar_metric_operator(_diag(1, 1, 1), _zero(3), _jet(gradient=(float("nan"), 0, 0)))
        except ValueError:
            rejections["nonfinite_jet"] = 1
        else:
            rejections["nonfinite_jet"] = 0
        record(current, rejections, {key: 1 for key in rejections})

    except Exception as exc:
        row = {"name": current, "passed": False, "exception_type": type(exc).__name__, "exception": str(exc)}
        if progress is not None:
            progress({"stage": "unreduced_algebra_control_exception", "control": current,
                      "passed": False, "control_result": row})
        raise
    return {"passed": bool(rows) and all(row["passed"] for row in rows), "rows": rows,
            "maximum": max(row["maximum"] for row in rows), "control_count": len(rows),
            "check_count": sum(len(row["checks"]) for row in rows),
            "failure_count": sum(not row["passed"] for row in rows),
            "settings": {"scaled_tolerance": tolerance, "scale": "max(1, abs(reference))"},
            "accuracy_certified": False, "full_H6_validation": False,
            "integrated_flux_validation": False,
            "scope": "Manufactured pointwise recursion, nine spatial ledger entries, separate global term and oriented flux densities"}
