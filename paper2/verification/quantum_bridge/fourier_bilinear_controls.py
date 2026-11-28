"""Independent normalization, adjointness and derivative controls for H6 kernels.

All calculations occur inside run_fourier_controls. No source, kernel,
quadrature or test is constructed at import time. These are backend controls,
not gravitational amplitudes or a physical cutoff study.
"""

from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from math import pi

from gaussian_heat_sources import from_heat, from_polynomial
from gaussian_green_fields import HeatField, scalar_field, york_field
from gaussian_heat_pairings import PairingAccumulator
from gaussian_fourier_bilinears import compile_scalar_bilinear, compile_york_bilinear
from reduced_gaunt_exact import ONE, multiply, gaussian_moment


def _decimal(value):
    value = F(value)
    return D(value.numerator)/D(value.denominator)


def _pg_inner_normalized(left, right, rate, precision):
    """Direct Cartesian Gaussian moments, divided by pi^(3/2)."""
    with localcontext() as ctx:
        ctx.prec = precision
        rate = _decimal(rate)
        value = D(0)
        for powers, coefficient in multiply(left, right).items():
            moment = gaussian_moment(powers)
            if moment:
                value += _decimal(coefficient*moment)/(rate**(sum(powers)//2+1)*rate.sqrt())
        return value


def _sum_kernels(kernels):
    kernels = tuple(kernels)
    if not kernels:
        raise ValueError("Need at least one kernel")
    total = kernels[0]
    for kernel in kernels[1:]:
        total = total.add(kernel)
    return total


def run_fourier_controls(config, progress=None):
    options = dict(config.get("fourier_options", {}))
    options.setdefault("orders", [8, 16, 24])
    options.setdefault("precision", 64)
    options.setdefault("check_precision", 96)
    options.setdefault("absolute_tolerance", "1e-12")
    options.setdefault("relative_tolerance", "1e-10")
    options.setdefault("max_term_evaluations", 20_000_000)
    caps = dict(config.get("fourier_compile_caps", {}))
    exact_tolerance = D(str(config.get("fourier_exact_control_tolerance", "1e-40")))
    spatial_tolerance = float(config.get("fourier_spatial_control_tolerance", 1e-9))
    spatial_options = dict(config.get("fourier_spatial_options", {}))
    spatial_options.setdefault("quadrature_orders", [16, 24, 32])
    spatial_options.setdefault("absolute_tolerance", 1e-12)
    spatial_options.setdefault("relative_tolerance", 1e-10)
    precision = options["check_precision"]
    rows = []

    def record(name, kernel, expected=None, spatial=None, factor=F(1)):
        result = kernel.integrate(**options)
        with localcontext() as ctx:
            ctx.prec = precision
            actual = D(result["normalized_value_decimal"])
            if spatial is None:
                expected_value = expected
                error = abs(actual-expected_value)/max(D(1), abs(expected_value))
                passed = result["converged"] and error < exact_tolerance
                reference = {"normalized_decimal": str(expected_value), "method": "direct analytic Gaussian moments"}
            else:
                reference = spatial.integrate(**spatial_options)
                expected_full = float(factor)*reference["value"]
                full = result["value"]
                if full is None:
                    raise ArithmeticError("Control result is outside display/assembly range")
                error = abs(full-expected_full)/max(1.0, abs(expected_full))
                passed = result["converged"] and reference["converged"] and error < spatial_tolerance
                reference = {"spatial_integral": reference, "factor": str(factor),
                             "expected_full_value": expected_full,
                             "method": "independent Cartesian angular/HEAT radial pairing"}
            row = {"name": name, "kernel": result, "reference": reference,
                   "scaled_error": str(error), "passed": passed}
        rows.append(row)
        if progress:
            progress({**row, "completed": len(rows)})
        if not passed:
            raise AssertionError(row)
        return result

    a, b = F(1), F(3, 2)
    f, g = from_polynomial(ONE, a), from_polynomial(ONE, b)
    with localcontext() as ctx:
        ctx.prec = precision
        base = 1/(_decimal(a)*_decimal(b)*_decimal(a+b).sqrt())
        gaussian_scalar = base/16
        c, d = (F(1), F(2), F(-1)), (F(-2), F(1), F(3))
        gaussian_york = -D(11)*_decimal(sum(x*y for x, y in zip(c, d)))*base/24
    record("scalar_gaussian", compile_scalar_bilinear(f, g, **caps), expected=gaussian_scalar)
    record("york_gaussian_vectors", compile_york_bilinear(
        tuple(f.scale(x) for x in c), tuple(g.scale(x) for x in d), **caps), expected=gaussian_york)

    pf = {(0, 0, 0): F(1), (1, 1, 0): F(1), (0, 0, 2): F(1, 3)}
    pg = {(0, 0, 0): F(2), (1, 0, 1): F(-1), (2, 0, 0): F(1, 2)}
    f, g = from_polynomial(pf, a), from_polynomial(pg, b)
    df, dg = tuple(f.derivative(i) for i in range(3)), tuple(g.derivative(i) for i in range(3))
    ordinary = _pg_inner_normalized(pf, pg, a+b, precision)
    with localcontext() as ctx:
        ctx.prec = precision
        scalar_expected, york_expected = ordinary/8, -3*ordinary/4
    record("pure_gaussian_gradient_scalar", _sum_kernels(
        compile_scalar_bilinear(x, y, **caps) for x, y in zip(df, dg)), expected=scalar_expected)
    record("pure_gaussian_gradient_york", compile_york_bilinear(df, dg, **caps), expected=york_expected)

    # A one-heat endpoint paired with a pure Gaussian can be checked by
    # moving the inverse to the pure endpoint and using known first fields.
    hg = HeatField(F(5, 4), {((0, 0, 0), 2): F(1), ((1, 0, 0), 1): F(1, 3)})
    pg_heat = {(0, 0, 0): F(1), (0, 1, 0): F(1, 2)}
    heat_source = from_heat(hg, pg_heat, F(3, 4))
    scalar_spatial = PairingAccumulator()
    scalar_spatial.add(pg_heat, F(3, 4), fields=(scalar_field(pf, a), hg))
    k = compile_scalar_bilinear(f, heat_source, **caps)
    reverse = compile_scalar_bilinear(heat_source, f, **caps)
    assert k.blocks == reverse.blocks, "Full endpoint reciprocity failed"
    assert not k.add(k, factor=-1).blocks, "Exact signed kernel merge failed"
    record("one_heat_scalar_adjoint", k, spatial=scalar_spatial)
    bp = (pf, {(1, 0, 0): F(1)}, {(0, 1, 1): F(-1)})
    bsource = tuple(from_polynomial(poly, a) for poly in bp)
    heat_vector = (heat_source, heat_source.scale(-2), heat_source.scale(F(1, 2)))
    inverse_vector = york_field(bp, a)
    vector_spatial = PairingAccumulator()
    for field, coefficient in zip(inverse_vector, (F(1), F(-2), F(1, 2))):
        vector_spatial.add(pg_heat, F(3, 4), fields=(field, hg), factor=coefficient)
    record("one_heat_york_adjoint", compile_york_bilinear(bsource, heat_vector, **caps), spatial=vector_spatial)

    # A genuine 2D parameter control: unequal moments/rates and nonradial
    # polynomials. Spatial derivative identities remove the inverse exactly.
    hf = HeatField(F(1, 3), {((0, 0, 0), 1): F(1)})
    hg2 = HeatField(F(5, 4), {((0, 0, 0), 2): F(1)})
    p2f = {(0, 0, 0): F(1), (1, 0, 0): F(1)}
    p2g = {(0, 0, 0): F(1), (1, 0, 0): F(1), (0, 2, 0): F(1)}
    f2, g2 = from_heat(hf, p2f, F(1, 2)), from_heat(hg2, p2g, F(3, 4))
    df2, dg2 = tuple(f2.derivative(i) for i in range(3)), tuple(g2.derivative(i) for i in range(3))
    product_spatial = PairingAccumulator()
    product_spatial.add(multiply(p2f, p2g), F(5, 4), fields=(hf, hg2))
    record("two_heat_gradient_scalar", _sum_kernels(
        compile_scalar_bilinear(x, y, **caps) for x, y in zip(df2, dg2)),
        spatial=product_spatial, factor=F(1, 8))
    record("two_heat_gradient_york", compile_york_bilinear(df2, dg2, **caps),
           spatial=product_spatial, factor=-F(3, 4))

    # The inverse has negative eigenvalues between -1/k² and -3/(4k²).
    vector = (f2, g2, f2.add(g2, factor=-F(1, 3)))
    energy_kernel = _sum_kernels(compile_scalar_bilinear(x, x, **caps) for x in vector).scale(8)
    energy = energy_kernel.integrate(**options)
    york = compile_york_bilinear(vector, vector, **caps).integrate(**options)
    with localcontext() as ctx:
        ctx.prec = precision
        e, y = D(energy["normalized_value_decimal"]), D(york["normalized_value_decimal"])
        passed = energy["converged"] and york["converged"] and e > 0 and -e <= y <= -3*e/4
    row = {"name": "york_negative_spectral_bounds", "energy": energy, "york": york,
           "passed": passed, "inequality": "-E <= York <= -3E/4, for this real self pairing"}
    rows.append(row)
    if progress:
        progress({**row, "completed": len(rows)})
    if not passed:
        raise AssertionError(row)
    return {"passed": all(row["passed"] for row in rows), "rows": rows,
            "cases": len(rows), "settings": {"fourier_options": options,
                "compile_caps": caps, "spatial_options": spatial_options,
                "exact_control_tolerance": str(exact_tolerance),
                "spatial_control_tolerance": spatial_tolerance},
            "scope": "Normalization, derivative and adjoint backend controls; no gravitational H6 or four-point result"}
