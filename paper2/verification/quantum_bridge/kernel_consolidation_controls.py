"""Tracked exact-representation and independent integral checks for regrouping.

Imports are inert. Equal-endpoint canonical kernels represent the same
integral modulo termwise variable swaps: compare symmetrized integrands,
not their potentially different pointwise representatives. No quadrature
variable or measure is changed by the exact consolidation.
"""

from contextlib import contextmanager
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from math import isfinite

from gaussian_fourier_bilinears import (
    GaussianBilinearKernel, compile_scalar_bilinear, compile_york_bilinear,
)
from gaussian_green_fields import HeatField
from gaussian_heat_pairings import PairingAccumulator
from gaussian_heat_sources import from_heat
from gaussian_kernel_consolidation import consolidate_kernel
from reduced_gaunt_exact import multiply


def kernel_record(kernel):
    """Lossless exact representation for the immutable run's output files."""
    return {"statistics": kernel.statistics(), "blocks": [
        {"endpoints": [{"beta": str(beta), "heat_rate": None if a is None else str(a)}
                       for beta, a in block],
         "terms": [{"powers": list(term), "coefficient": str(c)} for term, c in terms.items()]}
        for block, terms in kernel.blocks.items()]}


def _rational_block(block, terms, t, u):
    """The exact block integrand multiplied by its common sqrt(A+B)."""
    A = block[0][0] if block[0][1] is None else block[0][0]+block[0][1]*t*t
    B = block[1][0] if block[1][1] is None else block[1][0]+block[1][1]*u*u
    axes = (t, u, A, B, A+B)
    powers = []
    for axis, value in enumerate(axes):
        exponents = {term[axis] for term in terms}
        if axis < 2:
            powers.append({e: F(1) if e == -1 else value**(2*e) for e in exponents})
        elif axis == 4:
            powers.append({e: value**(-e) for e in exponents})
        else:
            powers.append({e: value**e for e in exponents})
    total = F(0)
    for term, c in terms.items():
        for axis in range(5):
            c *= powers[axis][term[axis]]
        total += c
    return total


def check_consolidation(raw, compact, *, points=None, max_exact_term_visits=1_000_000):
    """Held-out exact rational checks, with explicit canonical-swap semantics.

    These sampled equalities complement the exact algebraic rewrite proof;
    they do not by themselves prove a rational-function identity globally.
    """
    if points is None:
        points = ((0, 0), (F(1, 4), F(2, 3)), (F(2, 5), F(1, 7)), (F(3, 4), F(1, 4)), (1, 1))
    points = tuple(tuple(F(x) for x in point) for point in points)
    if not points or any(len(p) != 2 or any(x < 0 or x > 1 for x in p) for p in points):
        raise ValueError("Nonempty rational points in the unit square are required")
    if type(max_exact_term_visits) is not int or max_exact_term_visits < 1:
        raise ValueError("max_exact_term_visits must be positive")
    if not set(compact.blocks).issubset(raw.blocks):
        raise ValueError("Consolidation introduced a different endpoint descriptor")
    planned = sum((len(terms)+len(compact.blocks.get(block, {})))*len(points)
                  *(2 if block[0] == block[1] else 1) for block, terms in raw.blocks.items())
    if planned > max_exact_term_visits:
        raise ValueError("Exact representation checks exceed their declared visit budget")
    rows = []
    for index, (block, terms) in enumerate(raw.blocks.items()):
        new_terms = compact.blocks.get(block, {})
        symmetric = block[0] == block[1]
        for t, u in points:
            old_direct = _rational_block(block, terms, t, u)
            new_direct = _rational_block(block, new_terms, t, u)
            old, new = old_direct, new_direct
            if symmetric:
                old = (old+_rational_block(block, terms, u, t))/2
                new = (new+_rational_block(block, new_terms, u, t))/2
            rows.append({"block": index, "point": [str(t), str(u)], "symmetrized": symmetric,
                         "raw_rational": str(old), "compact_rational": str(new),
                         "residual": str(new-old), "passed": new == old,
                         "unsymmetrized_representatives_differ": new_direct != old_direct})
    return {"passed": all(row["passed"] for row in rows), "rows": rows,
            "check_count": len(rows), "planned_exact_term_visits": planned,
            "scope": "Exact sampled stripped block integrands; symmetrized for equal endpoints"}


def run_consolidation_controls(config, progress=None):
    """Run small exact controls and genuine two-parameter integral comparisons."""
    caps = dict(config.get("consolidation_caps", {}))
    options = dict(config["fourier_options"])
    comparison_tolerance = D(str(config.get("consolidation_comparison_tolerance", "1e-10")))
    spatial_tolerance = float(config.get("consolidation_spatial_tolerance", 1e-9))
    if (not comparison_tolerance.is_finite() or comparison_tolerance <= 0
            or not isfinite(spatial_tolerance) or spatial_tolerance <= 0):
        raise ValueError("Require positive finite control tolerances")
    rows = []
    def record(row, *, fail=True):
        rows.append(row)
        if progress:
            progress({"stage": "kernel_consolidation_control", "completed": len(rows), "of": 10,
                      "control_result": row})
        if fail and not row["passed"]:
            raise AssertionError("Kernel consolidation control failed: "+row["name"])

    @contextmanager
    def control(name):
        try:
            yield
        except Exception as exc:
            if not rows or rows[-1]["name"] != name:
                record({"name": name, "passed": False,
                        "failure": {"type": type(exc).__name__, "message": str(exc)}}, fail=False)
            raise

    fixed = ((F(3, 2), None), (F(5, 4), None))
    one = ((F(3, 2), None), (F(1, 2), F(2)))
    unequal = ((F(1, 2), F(1, 3)), (F(3, 4), F(5, 4)))
    equal = ((F(1, 2), F(1)), (F(1, 2), F(1)))
    cases = [
        ("empty_kernel", GaussianBilinearKernel(), True),
        ("fixed_gaussian_widths", GaussianBilinearKernel({fixed: {
            (-1,-1,2,-1,3): F(3,7), (-1,-1,-2,1,0): F(-2), (-1,-1,0,0,1): F(5)}}), False),
        ("one_parameter_unequal_rate", GaussianBilinearKernel({one: {
            (-1,2,-2,-1,2): F(3), (-1,0,1,2,0): F(-1,3), (-1,1,0,0,1): F(2)}}), False),
        ("two_parameter_unequal_rates", GaussianBilinearKernel({unequal: {
            (3,1,-2,1,2): F(5), (1,2,2,-3,1): F(-2,7), (0,0,0,0,3): F(1)}}), False),
        ("equal_endpoint_orientation", GaussianBilinearKernel({equal: {(1,0,-1,2,1): F(1)}}), False),
        ("sum_width_identity_zero", GaussianBilinearKernel({unequal: {
            (0,0,0,0,0): F(1), (0,0,1,0,1): F(-1), (0,0,0,1,1): F(-1)}}), True),
        ("heat_weight_identity_zero", GaussianBilinearKernel({unequal: {
            (1,0,0,0,0): F(1,3), (0,0,1,0,0): F(-1), (0,0,0,0,0): F(1,2)}}), True),
    ]
    for name, raw, expect_zero in cases:
        with control(name):
            original = {block: dict(terms) for block, terms in raw.blocks.items()}
            compact, stats = consolidate_kernel(raw, **caps)
            identity = check_consolidation(raw, compact)
            repeated, repeated_stats = consolidate_kernel(compact, **caps)
            immutable = raw.blocks == original
            idempotent = repeated.blocks == compact.blocks
            passed = (identity["passed"] and immutable and idempotent
                      and (not compact.blocks if expect_zero else bool(compact.blocks)))
            if name == "equal_endpoint_orientation":
                passed = passed and any(row["unsymmetrized_representatives_differ"] for row in identity["rows"])
            record({"name": name, "passed": passed, "identity": identity, "statistics": stats,
                    "raw_map_unchanged": immutable, "compact_idempotent": idempotent,
                    "repeat_statistics": repeated_stats,
                    "raw": kernel_record(raw), "compact": kernel_record(compact)})

    with control("exact_work_boundary"):
        raw = cases[3][1]
        compact, stats = consolidate_kernel(raw, **caps)
        work = stats["work"]["total"]
        if work < 2:
            raise AssertionError("Nontrivial exact-budget control needs at least two operations")
        exact_caps = {**caps, "max_expansion_products": work}
        repeated, repeated_stats = consolidate_kernel(raw, **exact_caps)
        rejected = False
        try:
            consolidate_kernel(raw, **{**caps, "max_expansion_products": work-1})
        except (ValueError, RuntimeError) as exc:
            rejected = "budget" in str(exc).lower() or "max_expansion_products" in str(exc)
            rejection = {"type": type(exc).__name__, "message": str(exc)}
        else:
            rejection = None
        record({"name": "exact_work_boundary", "passed": rejected and compact.blocks == repeated.blocks
                and repeated_stats["work"]["total"] == work, "work": work, "rejection": rejection})

    # Independent identity removes the inverse, leaving a spatial two-HEAT
    # integral. Unequal rates/moments test both parameter axes and derivatives.
    with control("two_parameter_setup"):
        hf = HeatField(F(1,3), {((0,0,0),1): F(1)})
        hg = HeatField(F(5,4), {((0,0,0),2): F(1)})
        pf = {(0,0,0): F(1), (1,0,0): F(1)}
        pg = {(0,0,0): F(1), (1,0,0): F(1), (0,2,0): F(1)}
        f, g = from_heat(hf,pf,F(1,2)), from_heat(hg,pg,F(3,4))
        df, dg = tuple(f.derivative(i) for i in range(3)), tuple(g.derivative(i) for i in range(3))
        scalar = GaussianBilinearKernel()
        for left, right in zip(df,dg):
            scalar = scalar.add(compile_scalar_bilinear(left,right))
        york = compile_york_bilinear(df,dg)
        spatial = PairingAccumulator().add(multiply(pf,pg),F(5,4),fields=(hf,hg)).integrate(
            quadrature_orders=(16,24,32),absolute_tolerance=1e-12,relative_tolerance=1e-10)
    for name, raw, factor in (("two_parameter_scalar_identity",scalar,F(1,8)),
                              ("two_parameter_york_identity",york,-F(3,4))):
        with control(name):
            compact, stats = consolidate_kernel(raw, **caps)
            identity = check_consolidation(raw,compact)
            old, new = raw.integrate(**options), compact.integrate(**options)
            with localcontext() as ctx:
                ctx.prec = options["check_precision"]
                old_value, new_value = D(old["normalized_value_decimal"]), D(new["normalized_value_decimal"])
                difference = abs(old_value-new_value)/max(D(1),abs(old_value))
            expected = float(factor)*spatial["value"]
            full = new["value"]
            error = None if full is None else abs(full-expected)/max(1,abs(expected))
            passed = (identity["passed"] and old["converged"] and new["converged"] and spatial["converged"]
                      and difference <= comparison_tolerance and error is not None and error <= spatial_tolerance)
            record({"name": name, "passed": passed, "identity": identity, "statistics": stats,
                    "raw_integration": old, "compact_integration": new, "spatial_reference": spatial,
                    "expected_physical_value": expected, "spatial_scaled_error": error,
                    "raw_compact_scaled_difference_decimal": str(difference)})
    return {"passed": all(row["passed"] for row in rows), "rows": rows, "cases": len(rows),
            "settings": {"caps": caps, "fourier_options": options,
                         "comparison_tolerance": str(comparison_tolerance), "spatial_tolerance": spatial_tolerance},
            "accuracy_certified": False,
            "scope": "Exact rewrite and small independent integral controls; production H6 remains separate"}
