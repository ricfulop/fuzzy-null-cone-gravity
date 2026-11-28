"""Complete twenty-entry reduced H6 assembly for supplied raw TT directions.

This is an implementation of the audited integrated formula, not a claim
that its full numerical or independent mixed-functional checks have passed.
All construction, compilation and integration require explicit calls.
The three backends have different normalized units. Only full physical
values and physical error indicators enter the final energy ledger.
"""

from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from math import fsum, isfinite, nextafter, pi
from sys import float_info

from gaussian_green_fields import HeatField, scalar_field, york_field
from gaussian_fourier_bilinears import compile_scalar_bilinear, compile_york_bilinear
from higher_local_geometry import local_geometry
from reduced_gaunt_exact import add, tensor_square
from reduced_sixth_sources import prepare_sixth_sources
from reduced_sixth_local import local_families, merge_families, triple_tail
from sixth_canonical_terms import canonical_terms
from triple_heat_integration import integrate_triple_heat


def prepare_sixth_data(q, s, progress=None):
    """Compute R4/L6/b3 once and construct the same first responses as H5.

    The combined geometry routine shares the connection for scalar and
    York data. No earlier lower-order geometry calculation is repeated.
    Returned fields retain the established H5 dictionary convention.
    """
    geometry = local_geometry(q, s, local_R_order=4, density_order=6,
                              york_order=3, progress=progress)
    york = {key: geometry[key] for key in ("E", "b", "york_order")}
    s2, r2 = tensor_square(s), geometry["R"][2]
    u = scalar_field(add(s2, r2, -1), F(1))
    v = scalar_field(add(r2, s2, 7), F(1))
    W = york_field(tuple(york["b"][1]), F(1))
    dW = [[W[k].derivative(a) for k in range(3)] for a in range(3)]
    divergence = HeatField(F(1), {})
    for k in range(3):
        divergence = divergence.add(dW[k][k])
    T = [[dW[i][j].add(dW[j][i]).add(divergence, -F(2, 3) if i == j else F(0))
          for j in range(3)] for i in range(3)]
    return {"q": q, "s": s, "geometry": geometry, "york": york, "s2": s2,
            "u": u, "v": v, "w": W, "dw": dW, "divergence_w": divergence,
            "t0": T, "grad_u": u.gradient(), "grad_v": v.gradient()}


def compile_nested_sources(sources, *, compiler_options=None, progress=None):
    """Compile both production inverse terms with their outer signs exactly once."""
    options = dict(compiler_options or {})
    def callback(name):
        return None if progress is None else lambda row: progress({"family": name, **row})
    scalar = compile_scalar_bilinear(sources["K"], sources["F"],
        progress=callback("scalar_nested"), **options).scale(-1)
    york = compile_york_bilinear(sources["J_L"], sources["J_R"],
        progress=callback("york_nested"), **options).scale(-2)
    return {"scalar_nested": scalar, "york_nested": york}


def assemble_sixth(data, basis, labels=None, *, sources=None, nested=None,
                   compiler_options=None, canonical_options=None,
                   max_pair_products=20_000_000, progress=None):
    """Return exact local/nested/triple objects and the two canonical entries.

    Supplying sources/nested avoids recompilation during a tracked profile.
    The caller must ensure they belong to this data. No invariant projection,
    physical angular limit, or omission of higher auxiliary modes is made.
    """
    if sources is None:
        sources = prepare_sixth_sources(data)
    if nested is None:
        nested = compile_nested_sources(sources, compiler_options=compiler_options, progress=progress)
    if set(nested) != {"scalar_nested", "york_nested"}:
        raise ValueError("Both complete nested inverse families are required")
    canonical = canonical_terms(data, basis, labels=labels,
                                progress=progress, **dict(canonical_options or {}))
    local = local_families(data, sources, max_pair_products=max_pair_products, progress=progress)
    triple, triple_statistics = triple_tail(data, max_pair_products=max_pair_products, progress=progress)
    return {"local": local, "local_merged": merge_families(local, max_pair_products),
            "nested": nested, "canonical": canonical, "triple": triple,
            "sources": sources,
            "statistics": {"local": {key: value.statistics for key, value in local.items()},
                           "nested": {key: value.statistics() for key, value in nested.items()},
                           "triple_components": triple_statistics},
            "conventions": {"local": "full physical spatial integral",
                "nested": "normalized Decimal times pi^(3/2); outer -1/-2 already included",
                "triple": "normalized Decimal times 4*pi; minus sign and tensor multiplicities already included",
                "canonical": "exact rational times pi^(3/2)/2; no additional canonical coordinate rescaling"}}


def _finite(value, name):
    value = float(value)
    if not isfinite(value):
        raise ArithmeticError(name+" is not a finite physical value")
    return value


def physical_indicator(values, factor=1, name="error indicator"):
    """Sum/scale nonnegative diagnostics before upward binary64 conversion.

    Decimal inputs may be smaller than the binary64 subnormal range. Such
    positive physical indicators round UP to the least positive float, never
    disappear. This is a conversion safeguard, not error certification.
    """
    def decimal(value):
        return D.from_float(value) if isinstance(value, float) else D(value)
    terms, scale = tuple(decimal(value) for value in values), decimal(factor)
    if not scale.is_finite() or scale <= 0 or any(not term.is_finite() or term < 0 for term in terms):
        raise ArithmeticError(name+" requires finite nonnegative terms and a positive factor")
    with localcontext() as ctx:
        ctx.prec = max(128, max((len(term.as_tuple().digits) for term in (*terms, scale)), default=1)+32)
        # Ceiling also preserves a positive contribution too small to change
        # an already positive sum at the declared Decimal precision.
        ctx.rounding = "ROUND_CEILING"
        total = sum(terms, D(0))*scale
        converted = _finite(total, name)
        if total > D.from_float(converted):
            converted = nextafter(converted, float("inf"))
        return _finite(converted, name)


def integrate_sixth(assembly, *, local_options=None, fourier_options=None,
                    triple_options=None, absolute_tolerance=1e-7,
                    relative_tolerance=1e-7, roundoff_safety_factor=100,
                    progress=None):
    """Integrate every family and gate a summed physical absolute indicator.

    The indicator sums measured refinement, tail and arithmetic diagnostics;
    it is not a certified error bound. A numerical gate does not replace the
    independent mixed W2/u4 or spatial-reference validation of the functional.
    """
    atol, rtol = _finite(absolute_tolerance, "absolute_tolerance"), _finite(relative_tolerance, "relative_tolerance")
    safety = _finite(roundoff_safety_factor, "roundoff_safety_factor")
    if atol <= 0 or rtol < 0 or safety < 1:
        raise ValueError("Require positive atol, nonnegative rtol and roundoff factor >=1")
    ledger = []
    def append(name, value, error, converged, details):
        row = {"family": name, "physical_value": _finite(value, name),
               "absolute_error_indicator": _finite(error, name+" error"),
               "converged": bool(converged), "details": details}
        if row["absolute_error_indicator"] < 0:
            raise ValueError("Error indicators must be nonnegative")
        ledger.append(row)
        if progress:
            progress({"stage": "sixth_energy_family", "completed": len(ledger), "of": 20, **row})

    for name, accumulator in assembly["local"].items():
        result = accumulator.integrate(**dict(local_options or {}))
        bound = result["tail"]["absolute_bound"]
        if bound is None:
            raise ArithmeticError("Unrepresentable local tail bound")
        arithmetic = physical_indicator([result["machine_epsilon_times_bin_l1"]], safety)
        error = physical_indicator([result["termwise_refinement_change"], bound, arithmetic])
        append(name, result["value"], error, result["converged"], result)
    canonical = assembly["canonical"]
    for name in ("gram_square", "canonical_cross"):
        value = canonical[name+"_value"]
        append(name, value, physical_indicator([abs(value)], D.from_float(safety)*D.from_float(float_info.epsilon)), True,
               {"exact_factor": canonical[name+"_exact"], "unit": "pi^(3/2)/2"})
    for name, kernel in assembly["nested"].items():
        result = kernel.integrate(**dict(fourier_options or {}))
        error = physical_indicator([result[key] for key in (
            "parameter_termwise_change_decimal", "precision_termwise_change_decimal", "working_roundoff_proxy_decimal")], pi**1.5)
        value = result["value"]
        if value is None:
            raise ArithmeticError("Nested value is outside physical binary64 range")
        arithmetic = physical_indicator([abs(value)], D.from_float(safety)*D.from_float(float_info.epsilon))
        error = physical_indicator([error, arithmetic])
        append(name, value, error, result["converged"], result)
    triple = integrate_triple_heat(assembly["triple"], progress=progress, **dict(triple_options or {}))
    value = triple["value"]
    if value is None:
        raise ArithmeticError("Triple value is outside physical binary64 range")
    error = physical_indicator([triple["combined_indicator_normalized_decimal"]], 4*pi, "triple indicator")
    arithmetic = physical_indicator([abs(value)], D.from_float(safety)*D.from_float(float_info.epsilon))
    error = physical_indicator([error, arithmetic])
    append("triple_tail", value, error, triple["converged"], triple)
    if len(ledger) != 20 or len({row["family"] for row in ledger}) != 20:
        raise AssertionError("Full H6 ledger must contain twenty distinct signed entries")
    total = _finite(fsum(row["physical_value"] for row in ledger), "H6 total")
    absolute_sum = fsum(abs(row["physical_value"]) for row in ledger)
    arithmetic = physical_indicator([absolute_sum], D.from_float(safety)*D.from_float(float_info.epsilon))
    error = physical_indicator([row["absolute_error_indicator"] for row in ledger]+[arithmetic])
    budget = atol+rtol*abs(total)
    return {"value": total, "ledger": ledger, "absolute_family_sum": absolute_sum,
            "total_absolute_error_indicator": error, "total_acceptance_budget": budget,
            "converged": all(row["converged"] for row in ledger) and error <= budget,
            "accuracy_certified": False, "independent_mixed_functional_validation": False,
            "settings": {"absolute_tolerance_physical": atol, "relative_tolerance": rtol,
                         "roundoff_safety_factor": safety,
                         "indicator_conversion": "Decimal sum/scaling, upward binary64 rounding; positive underflow retained as least subnormal"},
            "scope": "Numerical assembly of the reduced H6 formula; independent functional validation and four-point loop remain separate"}
