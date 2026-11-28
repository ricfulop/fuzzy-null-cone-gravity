"""Requested-only controls of the independent noncompact scalar/York adapter.

run_noncompact_york_controls(config, budget, progress=None) uses the caller's
cumulative spherical_york_reference.WorkBudget. Import performs no source
construction, inverse solve or evaluation. This module imports no production
HEAT, Fourier, pairing, geometry, or reduced-Hamiltonian implementation.

Fixtures are C exp(-3 r^2/2), Y(c exp(-r^2)) with c=(1,2,-1), and
Y(grad exp(-r^2)). Both York cases actually compute V=P f, h=div V, chi=P h.
Targets never supply computed h, matching constants, or boundary moments.
For U=P exp(-beta r^2)=I_0(beta r^2)/(2 beta), independent positive Decimal
Gauss-Legendre moments I_m(z)=int_0^1 t^(2m) exp(-z t^2) dt give

  Y(c g)_i = -c_i (7 I_0+I_1)/(16 beta)
              -x_i(c.x)(I_1-I_2)/8,
  chi(c g) = -(c.x)(I_0-I_1)/(4 beta),
  V(grad g)=grad U, h=-g, chi=-U, Y(grad g)=-3 grad U/4.

The target Cartesian derivatives use product rules, independently of the
solid-harmonic jet bridge. GL refinement is a diagnostic, not a rigorous GL
error bound; the optional omitted positive-Gaussian tail is bounded separately.
All values must be finite. Structural/reference/stub rows are required at
every level; physical comparisons only at the finest declared radial level,
and the latest refinement for each physical quantity. Coarse raw failures
are retained. Single-level calls require all their physical comparisons.

Callbacks: noncompact_control_row contains row; source_result contains case,
intervals and the complete adapter event under source; level_result contains
the complete level under level, including on a stopped/exception path. Named
setup failures also emit preparation_result with preparation. Other adapter
progress is forwarded as noncompact_solver_progress. The caller owns
wall-time/memory limits, persistent checkpoints and any angular/radius slices.
This suite supplies primitive checks, not gravitational H6 validation.
"""

from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import product
from math import ceil, fsum, isfinite, sqrt
from traceback import format_exc

from run import gauss_nodes
from gaussian_radial_reference import GaussianRadialSource
from solid_harmonic_reference import AngularGrid
from spherical_green_reference import ExteriorSource, OriginSeries, RadialGrid
from spherical_york_reference import (
    ScalarSourceSpec, SumProviders, WorkBudget, solve_scalar_source,
    solve_york_source,
)


def _real(value, name, positive=False):
    if isinstance(value, (bool, str, bytes)):
        raise ValueError(name + " must be a real number")
    value = float(value)
    if not isfinite(value) or (positive and value <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return value


def _integer(value, name, minimum=1, maximum=None):
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        raise ValueError(name + " is outside its integer range")
    return value


def _settings(config):
    levels = tuple(_integer(n, "noncompact_intervals", 24)
                   for n in config.get("noncompact_intervals", (128, 256, 512)))
    orders = tuple(_integer(n, "noncompact_quadrature_orders", 8, 256)
                   for n in config.get("noncompact_quadrature_orders", (24, 48)))
    if not levels or any(a >= b for a, b in zip(levels, levels[1:])):
        raise ValueError("noncompact_intervals must be nonempty and strictly increasing")
    if len(orders) < 2 or any(a >= b for a, b in zip(orders, orders[1:])):
        raise ValueError("noncompact_quadrature_orders requires >=2 increasing orders")
    settings = {"intervals": levels, "quadrature_orders": orders}
    for name, default in (("radius", 4), ("origin_radius", .5),
                          ("value_tolerance", 3e-4), ("component_tolerance", 3e-3),
                          ("structural_tolerance", 1e-10), ("reference_tolerance", 1e-22)):
        settings[name] = _real(config.get("noncompact_"+name, default), name, True)
    for name, default, lower, upper in (
        ("n_mu", 8, 4, 128), ("n_phi", 16, 7, 256),
        ("h_lmax", 2, 2, 8), ("guard_lmax", 3, 3, 9),
        ("origin_degree", 16, 4, 32), ("decimal_precision", 64, 56, 256),
        ("max_reference_nodes", 200_000, 1, None),
        ("max_reference_radii", 128, 1, None),
    ):
        settings[name] = _integer(config.get("noncompact_"+name, default), name, lower, upper)
    if settings["origin_radius"] != .5 or not 2 <= settings["radius"] <= 16:
        raise ValueError("controls require fixed patch .5 and 2<=R<=16")
    if (settings["guard_lmax"] <= settings["h_lmax"]
            or settings["n_mu"] < settings["guard_lmax"]+1
            or settings["n_phi"] < 2*settings["guard_lmax"]+1):
        raise ValueError("angular grid must cover a strictly larger full guard band")
    provider = {}
    for name, default, lower in (("precision", 64, 32), ("max_products", 5_000_000, 1),
            ("max_terms", 100_000, 1), ("max_moment_steps", 2_000_000, 1),
            ("max_series_terms", 2048, 16), ("max_point_evaluations", 100_000, 1),
            ("max_cached_radii", 256, 1), ("max_z", 10_000, 1)):
        provider[name] = _integer(config.get("noncompact_provider_"+name, default), name, lower)
    settings["provider"] = provider
    settings.update({"source_rates": {"scalar_gaussian": "3/2", "vector_gaussian": "1",
                                      "gradient_gaussian": "1"},
                     "vector": (1, 2, -1), "target_mmax": 4,
                     "target_max_panel_width": 1.5, "origin_derivative_order": 4})
    return settings


def _plain(value):
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(item) for item in value]
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, Decimal):
        return str(value)
    if isinstance(value, float):
        return _real(value, "JSON value")
    return value


def _emit(progress, stage, **payload):
    if progress is not None:
        progress({"stage": stage, **_plain(payload)})


def _flatten(value):
    if isinstance(value, (tuple, list)):
        return [number for item in value for number in _flatten(item)]
    return [_real(value, "comparison value")]


def _row(rows, progress, case, n, name, actual, target, tolerance,
         *, role="structural", required=True, **extra):
    left, right = _flatten(actual), _flatten(target)
    if len(left) != len(right):
        raise ValueError("comparison lengths differ: "+name)
    errors = [abs(a-b)/max(1.0, abs(a), abs(b)) for a, b in zip(left, right)]
    error = _real(max(errors, default=0.0), "comparison residual")
    row = {"case": case, "intervals": n, "check": name, "actual": _plain(actual),
           "target": _plain(target), "scaled_residual": error,
           "tolerance": tolerance, "passed": error <= tolerance,
           "role": role, "acceptance_required": required, **_plain(extra)}
    rows.append(row)
    _emit(progress, "noncompact_control_row", row=row)
    return row


class _MomentReference:
    """Positive Decimal integration; no recurrence or production field target."""
    def __init__(self, settings, budget, rows, progress):
        self.settings, self.budget, self.rows, self.progress = settings, budget, rows, progress
        self.cache, self.nodes = {}, {}
        self.node_evaluations = 0

    def _integrate(self, z, order):
        precision = self.settings["decimal_precision"]
        if order not in self.nodes:
            self.nodes[order] = gauss_nodes(order, precision)
        if z <= 1:
            upper, rate, factor = Decimal(1), z, Decimal(1)
        else:
            root = z.sqrt()
            upper, rate, factor = min(root, Decimal(12)), Decimal(1), 1/root
        panels = max(1, ceil(float(upper)/1.5))
        width = upper/panels
        needed = panels*order
        if self.node_evaluations+needed > self.settings["max_reference_nodes"]:
            raise RuntimeError("noncompact positive-GL reference node cap exceeded")
        self.budget.charge("point_evaluations", needed)
        self.node_evaluations += needed
        integrals = [Decimal(0)]*5
        for panel in range(panels):
            middle = width*(Decimal(panel)+Decimal("0.5"))
            for node, weight in self.nodes[order]:
                u = middle+width*node/2
                u2 = u*u
                term = width*weight*(-rate*u2).exp()
                for m in range(5):
                    integrals[m] += term
                    term *= u2
        moments, tails = [], []
        boundary = upper*(-upper*upper).exp()
        for m in range(5):
            moments.append(integrals[m]*factor)
            truncated = z > 1 and upper < z.sqrt()
            tails.append(boundary*factor/(2*(upper*upper-m)) if truncated else Decimal(0))
            if z > 1:
                factor /= z
            boundary *= upper*upper
        return tuple(moments), tuple(tails), panels

    def moments(self, beta, point):
        with localcontext() as ctx:
            ctx.prec = self.settings["decimal_precision"]
            b = Decimal(beta.numerator)/Decimal(beta.denominator)
            x = tuple(Decimal.from_float(float(value)) for value in point)
            z = b*sum((value*value for value in x), Decimal(0))
            key = str(z)
            if key not in self.cache:
                if len(self.cache) >= self.settings["max_reference_radii"]:
                    raise RuntimeError("noncompact positive-GL reference radius cap exceeded")
                records = [self._integrate(z, order) for order in self.settings["quadrature_orders"]]
                coarse, fine = records[-2][0], records[-1][0]
                if any(value <= 0 for value in fine):
                    raise ArithmeticError("positive reference moment lost positivity")
                difference = tuple(abs(a-b)/b for a, b in zip(coarse, fine))
                tail = tuple(a/b for a, b in zip(records[-1][1], fine))
                for name, errors in (("quadrature_refinement", difference), ("omitted_tail", tail)):
                    _row(self.rows, self.progress, "independent_target", None, name,
                         tuple(float(v) for v in errors), (0.0,)*5,
                         self.settings["reference_tolerance"], role="reference", z=str(z),
                         moments=[str(value) for value in fine],
                         errors_decimal=[str(value) for value in errors],
                         panels=[record[2] for record in records])
                self.cache[key] = fine
            return x, b, self.cache[key]

    def statistics(self):
        return {"positive_gl_nodes": self.node_evaluations, "distinct_z": len(self.cache),
                "orders": self.settings["quadrature_orders"], "mmax": 4,
                "charged_to_shared_point_budget": True,
                "quadrature_difference_is_not_a_rigorous_error_bound": True}


def _radial_jet(x, derivatives):
    g, g1, g2, g3 = derivatives
    return {"value": g, "gradient": tuple(2*t*g1 for t in x),
            "hessian": tuple(tuple(2*int(a == b)*g1+4*x[a]*x[b]*g2 for b in range(3))
                             for a in range(3)),
            "third": tuple(tuple(tuple(4*(int(a == b)*x[c]+int(a == c)*x[b]+int(b == c)*x[a])*g2
                                       +8*x[a]*x[b]*x[c]*g3 for c in range(3))
                                 for b in range(3)) for a in range(3))}


def _linear_radial_jet(x, vector, derivatives):
    radial = _radial_jet(x, derivatives)
    dot = sum((v*t for v, t in zip(vector, x)), Decimal(0))
    return {"value": dot*radial["value"],
            "gradient": tuple(vector[a]*radial["value"]+dot*radial["gradient"][a] for a in range(3)),
            "hessian": tuple(tuple(vector[a]*radial["gradient"][b]+vector[b]*radial["gradient"][a]
                                    +dot*radial["hessian"][a][b] for b in range(3)) for a in range(3)),
            "third": tuple(tuple(tuple(dot*radial["third"][a][b][c]
                         +vector[a]*radial["hessian"][b][c]+vector[b]*radial["hessian"][a][c]
                         +vector[c]*radial["hessian"][a][b] for c in range(3))
                         for b in range(3)) for a in range(3))}


def _numeric_tree(value):
    if isinstance(value, dict):
        return {key: _numeric_tree(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return tuple(_numeric_tree(item) for item in value)
    return _real(value, "target value")


def _scaled_tree(value, factor):
    if isinstance(value, dict):
        return {key: _scaled_tree(item, factor) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return tuple(_scaled_tree(item, factor) for item in value)
    return factor*value


def _targets(case, beta, point, reference):
    with localcontext() as ctx:
        ctx.prec = reference.settings["decimal_precision"]
        x, b, moments = reference.moments(beta, point)
        u = _radial_jet(x, tuple((-b)**k*moments[k]/(2*b) for k in range(4)))
        if case == "scalar_gaussian":
            return {"field": _numeric_tree(_scaled_tree(u, Decimal(1)/8))}
        if case == "gradient_gaussian":
            v = tuple({"value": u["gradient"][i],
                       "gradient": tuple(u["hessian"][a][i] for a in range(3)),
                       "hessian": tuple(tuple(u["third"][a][b][i] for b in range(3)) for a in range(3))}
                      for i in range(3))
            w = {"value": tuple(-Decimal(3)/4*u["gradient"][i] for i in range(3)),
                 "gradient": tuple(tuple(-Decimal(3)/4*u["hessian"][a][i] for i in range(3)) for a in range(3)),
                 "hessian": tuple(tuple(tuple(-Decimal(3)/4*u["third"][a][b][i] for i in range(3))
                                          for b in range(3)) for a in range(3))}
            gaussian = (-b*sum((t*t for t in x), Decimal(0))).exp()
            h = _radial_jet(x, tuple(-(-b)**k*gaussian for k in range(4)))
            return _numeric_tree({"field": w, "V": v, "h": h, "chi": _scaled_tree(u, -1)})
        vector = (Decimal(1), Decimal(2), Decimal(-1))
        dot = sum((v*t for v, t in zip(vector, x)), Decimal(0))
        A = tuple(-(-b)**k*(7*moments[k]+moments[k+1])/(16*b) for k in range(3))
        B = tuple(-(-b)**k*(moments[k+1]-moments[k+2])/8 for k in range(3))
        w = {"value": tuple(vector[i]*A[0]+x[i]*dot*B[0] for i in range(3)),
             "gradient": tuple(tuple(2*vector[i]*x[a]*A[1]
                 +(int(i == a)*dot+x[i]*vector[a])*B[0]+2*x[i]*dot*x[a]*B[1]
                 for i in range(3)) for a in range(3)),
             "hessian": tuple(tuple(tuple(vector[i]*(2*int(a == b0)*A[1]+4*x[a]*x[b0]*A[2])
                 +(int(i == a)*vector[b0]+int(i == b0)*vector[a])*B[0]
                 +2*((int(i == a)*dot+x[i]*vector[a])*x[b0]
                     +(int(i == b0)*dot+x[i]*vector[b0])*x[a])*B[1]
                 +x[i]*dot*(2*int(a == b0)*B[1]+4*x[a]*x[b0]*B[2])
                 for i in range(3)) for b0 in range(3)) for a in range(3))}
        chi = _linear_radial_jet(x, vector,
                 tuple(-(-b)**k*(moments[k]-moments[k+1])/(4*b) for k in range(4)))
        h = {"value": sum((vector[i]*u["gradient"][i] for i in range(3)), Decimal(0)),
             "gradient": tuple(sum((vector[i]*u["hessian"][a][i] for i in range(3)), Decimal(0)) for a in range(3)),
             "hessian": tuple(tuple(sum((vector[i]*u["third"][a][b0][i] for i in range(3)), Decimal(0))
                                    for b0 in range(3)) for a in range(3))}
        return _numeric_tree({"field": w, "V": tuple(_scaled_tree(u, c) for c in vector), "h": h, "chi": chi})


class _StubProvider:
    """Protocol fixture only; never supplied to a differential-equation solve."""
    def __init__(self, *, radius=1, modes=((0, 0),), center=0.0, error=0.0, bounds=None):
        self.radius_exact = Fraction(radius)
        self.radius = float(radius)
        self.mode_keys, self.center, self.error = tuple(modes), center, error
        self.bounds = dict(bounds) if bounds is not None else {mode: 0.0 for mode in self.mode_keys}

    def outside_response_regular_jet(self, ell, m, radius, alpha, order=3, side=None):
        return (self.center,)+(0.0,)*order

    def response_error_bounds(self, ell, m, radius, alpha, order=3, side=None):
        return (self.error,)+(0.0,)*order

    def remainder_response_regular_jet(self, ell, m, radius, alpha, order=3):
        return (self.center,)+(0.0,)*order

    def origin_shift(self, ell, m, alpha):
        return self.center

    def algebraic_response_terms(self, ell, m, alpha):
        return {(ell+1, 0): self.center}

    def response_remainder_derivative(self, axis, alpha=1):
        return self

    def source_regular_jet(self, ell, m, radius, order=1):
        return (0.0,)*(order+1)

    def weighted_outer_bounds(self):
        return dict(self.bounds)

    def provenance(self):
        return {"kind": "exact_binary64_stub_only", "radius_exact": str(self.radius_exact)}


def _stub_checks(settings, rows, progress):
    def spec(policy, bounds, provider=None, radius=1.0):
        return ScalarSourceSpec({mode: (0.0, 0.0, 0.0) for mode in bounds},
                {mode: OriginSeries(radius/2, (0.0,), (0.0,)*3) for mode in bounds},
                ExteriorSource(radius, {}, bounds), residual_policy=policy,
                residual_provider=provider)

    def rejection(name, function):
        provenance = {}
        try:
            function()
            rejected = False
        except (TypeError, ValueError) as error:
            rejected = True
            provenance = {"exception": type(error).__name__, "message": str(error)}
        _row(rows, progress, "adapter_stub", None, name, int(rejected), 1, 0.0, **provenance)

    a, b = (0, 0), (1, 0)
    rejection("policy_zero_rejects_positive_omission", lambda: spec("zero", {a: 1.0}))
    rejection("policy_exact_requires_provider", lambda: spec("exact_correction", {a: 1.0}))
    rejection("policy_bounded_cannot_attach_central_correction",
              lambda: spec("bounded_omission", {a: 1.0}, _StubProvider()))
    rejection("policy_unknown_rejected", lambda: spec("implicit", {a: 0.0}))
    rejection("positive_bound_mode_coverage", lambda: spec("exact_correction", {a: 1.0, b: 1.0},
              _StubProvider(bounds={a: 1.0})))
    rejection("provider_bounds_keys", lambda: spec("exact_correction", {a: 1.0}, _StubProvider(bounds={})))
    rejection("provider_negative_bound", lambda: spec("exact_correction", {a: 1.0},
              _StubProvider(bounds={a: -1.0})))
    rejection("provider_nonfinite_bound", lambda: spec("exact_correction", {a: 1.0},
              _StubProvider(bounds={a: float("inf")})))
    rejection("provider_bound_must_be_covered", lambda: spec("exact_correction", {a: 1.0},
              _StubProvider(bounds={a: 2.0})))
    rejection("exact_anchor_mismatch", lambda: spec("exact_correction", {a: 0.0},
              _StubProvider(radius=Fraction(1, 10)), .1))
    rejection("sum_exact_anchor_mismatch", lambda: SumProviders(
              (_StubProvider(radius=Fraction(1, 10)), _StubProvider(radius=.1))))
    for name, candidate in (
        ("policy_zero_kept", spec("zero", {a: 0.0})),
        ("policy_bounded_omission_kept", spec("bounded_omission", {a: 1.0})),
        ("exact_binary64_anchor_attaches", spec("exact_correction", {a: 0.0}, _StubProvider(radius=.1), .1)),
    ):
        _row(rows, progress, "adapter_stub", None, name, 1, 1, 0.0,
             residual_policy=candidate.residual_policy,
             declared_bounds=[{"mode": key, "bound": value} for key, value in candidate.exterior.remainder_outer_bounds.items()])
    tiny = 2.0**-53
    summed = SumProviders((_StubProvider(center=1.0), _StubProvider(center=tiny)))
    center = summed.outside_response_regular_jet(0, 0, 2.0, 1, order=3)
    bounds = summed.response_error_bounds(0, 0, 2.0, 1, order=3)
    exact_error = abs(Fraction.from_float(center[0])-(Fraction(1)+Fraction.from_float(tiny)))
    _row(rows, progress, "adapter_stub", None, "central_fsum_rounding_in_bound",
         int(Fraction.from_float(bounds[0]) >= exact_error and bounds[0] > 0), 1, 0.0,
         central=center, bounds=bounds, exact_rounding_error=str(exact_error))
    summed = SumProviders((_StubProvider(error=1.0, bounds={a: 1.0}),
                           _StubProvider(error=tiny, bounds={a: tiny})))
    target = Fraction(1)+Fraction.from_float(tiny)
    error = summed.response_error_bounds(0, 0, 2.0, 1, order=0)[0]
    outer = summed.weighted_outer_bounds()[a]
    _row(rows, progress, "adapter_stub", None, "positive_bounds_round_up",
         (int(Fraction.from_float(error) >= target), int(Fraction.from_float(outer) >= target)), (1, 1), 0.0,
         response_bound=error, weighted_outer_bound=outer, exact_positive_sum=str(target))
    rejection("nonfinite_central_origin_rejected", lambda: SumProviders(
              (_StubProvider(center=float("nan")),)).origin_shift(0, 0, 1))


def _fixture_sources(case, settings):
    if case == "scalar_gaussian":
        beta, polynomials = Fraction(3, 2), ({(0, 0, 0): Fraction(1)},)
    elif case == "vector_gaussian":
        beta, polynomials = Fraction(1), tuple({(0, 0, 0): Fraction(c)} for c in settings["vector"])
    else:
        beta = Fraction(1)
        polynomials = tuple({tuple(int(a == i) for a in range(3)): Fraction(-2)} for i in range(3))
    return beta, tuple(GaussianRadialSource(poly, beta, ell_max=1, **settings["provider"]) for poly in polynomials)


def _source_spec(source, grid, settings, budget):
    count = len(source.mode_keys)*len(grid.radii)
    budget.charge("point_evaluations", count)
    budget.charge("stored_values", count)
    samples = source.sample(grid.radii, max_values=count)
    tail = source.exterior(grid.radius)
    origin, records, limits = {}, [], {}
    patch = settings["origin_radius"]
    index = grid.radii.index(patch)
    for mode in source.mode_keys:
        record = source.origin_data(*mode, degree=settings["origin_degree"],
                  patch_radius=Fraction.from_float(patch), derivative_order=settings["origin_derivative_order"])
        origin[mode] = OriginSeries(record["patch_radius"], record["coefficients"], record["error_bounds_float_s"])
        records.append({"mode": mode, **record})
        limits[mode] = {patch: (samples[mode][index], samples[mode][index])}
    spec = ScalarSourceSpec(samples, origin,
            ExteriorSource(grid.radius, {}, tail.weighted_outer_bounds()),
            residual_policy="exact_correction", residual_provider=tail, source_limits=limits,
            provenance={"kind": "synthetic_PG_fixture", "beta": str(source.beta),
                        "origin_data": _plain(records), "source_statistics": source.statistics()})
    return spec


def _points(settings):
    patch, radius = settings["origin_radius"], settings["radius"]
    return (("origin", (0.0, 0.0, 0.0), None),
            ("inside_patch", (.125, -.25, .125), None),
            ("patch_left", (patch, 0.0, 0.0), "left"),
            ("patch_right", (patch, 0.0, 0.0), "right"),
            ("interior_mixed", (.37, -.41, .29), None),
            ("interior_far", (1.2, -.9, .6), None),
            ("R_x_left", (radius, 0.0, 0.0), "left"),
            ("R_x_right", (radius, 0.0, 0.0), "right"),
            ("R_y_left", (0.0, radius, 0.0), "left"),
            ("R_y_right", (0.0, radius, 0.0), "right"),
            ("outside_axis", (1.5*radius, 0.0, 0.0), None),
            ("outside_mixed", (.8*radius, -.9*radius, .7*radius), None))


def _laurent_derivatives(field, order, budget):
    fields = {(): field}
    for degree in range(1, order+1):
        for axes in product(range(3), repeat=degree):
            previous = fields[axes[:-1]]
            reservation = 3*len(previous.terms)
            budget.charge("exact_products", reservation)
            fields[axes] = previous.derivative(axes[-1], max_products=max(1, reservation))
            budget.charge("stored_values", len(fields[axes].terms))
    return fields


def _laurent_jet(fields, point, order, budget):
    budget.charge("point_evaluations")
    result = {"value": fields[()].evaluate(point)}
    for degree, name in ((1, "gradient"), (2, "hessian"), (3, "third")):
        if degree > order:
            break
        def nest(prefix):
            if len(prefix) == degree:
                return fields[prefix].evaluate(point)
            return tuple(nest(prefix+(axis,)) for axis in range(3))
        result[name] = nest(())
    return result


def _sum_tree(left, right):
    if isinstance(left, (tuple, list)):
        return tuple(_sum_tree(a, b) for a, b in zip(left, right))
    return fsum((left, right))


def _vector_jet(scalar_jets, order):
    result = {"value": tuple(jet["value"] for jet in scalar_jets)}
    if order >= 1:
        result["gradient"] = tuple(tuple(jet["gradient"][a] for jet in scalar_jets) for a in range(3))
    if order >= 2:
        result["hessian"] = tuple(tuple(tuple(jet["hessian"][a][b] for jet in scalar_jets)
                                       for b in range(3)) for a in range(3))
    return result


def _level(case, n, sources, beta, settings, angular, budget, reference, rows, progress, finest):
    record = {"case": case, "intervals": n, "radius": settings["radius"],
              "source_checkpoints": [], "physical_samples": {}, "rows": [],
              "status": "started", "work_before": budget.snapshot()}
    start = len(rows)

    def checkpoint(event):
        stage = event.get("stage", "unknown")
        if stage in ("scalar_source_origin_gate", "york_exact_support_gate"):
            _row(rows, progress, case, n, stage, int(bool(event["passed"])), 1, 0.0,
                 adapter_event=event)
        if stage in ("scalar_source_base_result", "scalar_source_result", "york_source_result", "york_field_result"):
            record["source_checkpoints"].append(_plain(event))
            _emit(progress, "source_result", case=case, intervals=n, source=event)
        else:
            _emit(progress, "noncompact_solver_progress", case=case, intervals=n, event=event)

    def comparison(label, actual, target, *, physical=True, **extra):
        for name in target:
            tolerance = (settings["value_tolerance"] if name == "value" else settings["component_tolerance"]) if physical else settings["structural_tolerance"]
            key = label+"/"+name
            _row(rows, progress, case, n, key, actual[name], target[name], tolerance,
                 role="physical" if physical else "structural", required=finest if physical else True, **extra)
            if physical:
                record["physical_samples"][key] = {"value": _plain(actual[name]), "tolerance": tolerance}

    phase = "radial_grid_and_source_specs"
    try:
        # Refuse a grid which cannot fit before allocating its radii/samples.
        solve_modes = sum(len(source.mode_keys) for source in sources)
        if case != "scalar_gaussian":
            solve_modes += (settings["h_lmax"]+1)**2
        if solve_modes*(n+1) > budget.remaining("fd_unknowns"):
            raise RuntimeError("declared level cannot fit the remaining cumulative FD cap")
        grid = RadialGrid(settings["radius"], n, (settings["origin_radius"],))
        specs = tuple(_source_spec(source, grid, settings, budget) for source in sources)
        declared_bounds = [value for spec in specs for origin in spec.origin_series.values() for value in origin.error_bounds]
        _row(rows, progress, case, n, "source_origin_error_declarations", declared_bounds,
             [0.0]*len(declared_bounds), settings["structural_tolerance"])
        phase = "inverse_solves"
        if case == "scalar_gaussian":
            field = solve_scalar_source(grid, specs[0], inverse="C", budget=budget, progress=checkpoint)
            scalar_fields, scalar_names = (field,), ("field",)
            york_algebraic = None
        else:
            field = solve_york_source(grid, specs, angular, h_lmax=settings["h_lmax"],
                    guard_lmax=settings["guard_lmax"], support_tolerance=settings["structural_tolerance"],
                    budget=budget, progress=checkpoint)
            scalar_fields, scalar_names = (*field.V, field.chi), ("V0", "V1", "V2", "chi")
            york_algebraic = tuple(_laurent_derivatives(item, 2, budget) for item in field.algebraic)
            if case == "gradient_gaussian":
                values = field.h_source.samples.get((0, 0), ())
                _row(rows, progress, case, n, "localized_h00_preserved", int(any(value != 0 for value in values)), 1, 0.0)
        phase = "algebraic_derivative_preparation"
        algebraic = tuple(_laurent_derivatives(item.algebraic, 3, budget) for item in scalar_fields)
        for label, point, side in _points(settings):
            phase = "point_checks/"+label
            targets = _targets(case, beta, point, reference)
            actual = field.jet(point, order=3 if case == "scalar_gaussian" else 2, side=side)
            comparison(label+"/field", actual, targets["field"], point=point, side=side)
            if case != "scalar_gaussian":
                for i, first in enumerate(field.V):
                    comparison(label+"/V"+str(i), first.jet(point, order=2, side=side),
                               {key: value for key, value in targets["V"][i].items() if key != "third"}, point=point, side=side)
                comparison(label+"/h", field.divergence_jet(point, order=2, side=side),
                           {key: value for key, value in targets["h"].items() if key != "third"}, point=point, side=side)
                comparison(label+"/chi", field.chi.jet(point, order=3, side=side), targets["chi"], point=point, side=side)
            radius = sqrt(fsum(x*x for x in point))
            for scalar, name in zip(scalar_fields, scalar_names):
                bounds = [bound for mode in scalar.mode_keys for bound in scalar.correction_error_bounds(*mode, radius, order=3, side=side)]
                _row(rows, progress, case, n, label+"/"+name+"/correction_arithmetic_bounds",
                     bounds, [0.0]*len(bounds), settings["structural_tolerance"], point=point, side=side,
                     total_field_error_certified=False)
            if radius >= grid.radius and side != "left":
                for scalar, name, derivatives in zip(scalar_fields, scalar_names, algebraic):
                    full = scalar.jet(point, order=3, side=side)
                    remainder = scalar.exterior_remainder_jet(point, order=3)
                    power = _laurent_jet(derivatives, point, 3, budget)
                    reconstructed = {key: _sum_tree(power[key], remainder[key]) for key in power}
                    comparison(label+"/"+name+"/algebraic_plus_remainder", full, reconstructed, physical=False, point=point)
                if york_algebraic is not None:
                    power = _vector_jet(tuple(_laurent_jet(item, point, 2, budget) for item in york_algebraic), 2)
                    remainder = field.exterior_remainder_jet(point, order=2)
                    reconstructed = {key: _sum_tree(power[key], remainder[key]) for key in power}
                    comparison(label+"/York/algebraic_plus_remainder", actual, reconstructed, physical=False, point=point)
        phase = "final_diagnostics"
        record["diagnostics"] = field.diagnostics()
        record["source_statistics"] = [source.statistics() for source in sources]
        record["status"] = "completed"
    except Exception as error:
        record["status"] = "stopped"
        record["exception"] = {"type": type(error).__name__, "message": str(error),
                               "phase": phase, "traceback": format_exc()}
        _row(rows, progress, case, n, "level_exception", 1, 0, 0.0, exception=record["exception"])
    record["rows"] = rows[start:]
    record["passed"] = all(row["passed"] for row in record["rows"] if row["acceptance_required"])
    record["raw_passed"] = all(row["passed"] for row in record["rows"])
    record["work_after"] = budget.snapshot()
    _emit(progress, "level_result", level=record)
    return record


def _preparation_failure(case, phase, error, budget, rows, progress):
    """Preserve setup failures too; KeyboardInterrupt/SystemExit still propagate."""
    detail = {"type": type(error).__name__, "message": str(error),
              "phase": phase, "traceback": format_exc()}
    row = _row(rows, progress, case, None, "preparation_exception", 1, 0, 0.0, exception=detail)
    record = {"case": case, "intervals": None, "status": "stopped", "phase": phase,
              "exception": detail, "rows": [row], "passed": False, "raw_passed": False,
              "source_checkpoints": [], "physical_samples": {}, "work_after": budget.snapshot()}
    _emit(progress, "preparation_result", preparation=record)
    _emit(progress, "level_result", level=record)
    return record


def run_noncompact_york_controls(config, budget, progress=None):
    """Run only when called by the tracked driver; return JSON-compatible evidence.

    Config keys use the noncompact_ prefix; defaults and provider limits are
    returned in settings. budget must be the shared WorkBudget from the driver.
    Exceptions during a level become required failed rows, immediately emit a
    stopped level_result, and stop this suite rather than retrying beyond caps.
    Invalid configuration or a callback exception may propagate to the driver.
    """
    if not isinstance(budget, WorkBudget):
        raise TypeError("budget must be the caller's shared WorkBudget")
    if progress is not None and not callable(progress):
        raise TypeError("progress must be callable")
    settings = _settings(config)
    before = budget.snapshot()
    rows, levels, refinements = [], [], []
    reference, angular, stopped = None, None, False
    phase = "adapter_stub_guards"
    try:
        _stub_checks(settings, rows, progress)
        phase = "positive_GL_reference_preparation"
        reference = _MomentReference(settings, budget, rows, progress)
        phase = "full_m_angular_grid_preparation"
        angular = AngularGrid(settings["guard_lmax"], settings["n_mu"], settings["n_phi"])
    except Exception as error:
        levels.append(_preparation_failure("suite_preparation", phase, error, budget, rows, progress))
        stopped = True
    for case in ("scalar_gaussian", "vector_gaussian", "gradient_gaussian"):
        if stopped:
            break
        try:
            beta, sources = _fixture_sources(case, settings)
        except Exception as error:
            levels.append(_preparation_failure(case, "Gaussian_source_preparation", error, budget, rows, progress))
            stopped = True
            break
        previous = None
        for n in settings["intervals"]:
            fine = n == settings["intervals"][-1]
            level = _level(case, n, sources, beta, settings, angular, budget, reference, rows, progress, fine)
            levels.append(level)
            if level["status"] != "completed":
                stopped = True
                break
            if previous is not None:
                refinement_rows = []
                for key, current in level["physical_samples"].items():
                    prior = previous["physical_samples"][key]
                    refinement_rows.append(_row(rows, progress, case, n, "radial_refinement/"+key,
                        current["value"], prior["value"], current["tolerance"],
                        role="refinement", required=fine, previous_intervals=previous["intervals"]))
                result = {"case": case, "coarse_intervals": previous["intervals"], "fine_intervals": n,
                          "acceptance_required": fine, "rows": refinement_rows,
                          "passed": all(row["passed"] for row in refinement_rows)}
                refinements.append(result)
                _emit(progress, "refinement_result", refinement=result)
            previous = level
        if stopped:
            break
    required = [row for row in rows if row["acceptance_required"]]
    diagnostic = [row for row in rows if not row["acceptance_required"]]
    complete = len(levels) == 3*len(settings["intervals"]) and not stopped
    result = {"passed": complete and all(row["passed"] for row in required),
              "rows": rows, "levels": levels, "refinements": refinements,
              "settings": settings, "completed_all_levels": complete,
              "maximum": max((row["scaled_residual"] for row in rows), default=0.0),
              "maximum_required": max((row["scaled_residual"] for row in required), default=0.0),
              "required_failure_count": sum(not row["passed"] for row in required),
              "diagnostic_failure_count": sum(not row["passed"] for row in diagnostic),
              "work_before": before, "work_after": budget.snapshot(),
              "target_work": reference.statistics() if reference is not None else None,
              "angular": angular.statistics() if angular is not None else None,
              "primitive_scope": "Independent noncompact scalar and two-stage York adapter controls",
              "full_H6_validation": False,
              "limitations": [
                  "No gravitational source, W2/u3/u4 solve, or Hamiltonian integral is evaluated here.",
                  "Gaussian exterior companions and exact angular algebra are shared with the adapter under test; targets use independent Decimal GL moments.",
                  "FD and matching-a0 uncertainty require measured refinement; provider arithmetic bounds do not certify the total field error.",
                  "The local ODE supplies higher radial jets; its residual is not an independent PDE control.",
                  "This invocation refines only radius resolution; separate angular/radius slices belong to the driver.",
                  "No integrated pairing/norm/flux target is supplied by this bounded noncompact suite.",
                  "GL refinement is diagnostic; only the explicitly omitted positive-integrand tail has the stated analytic bound.",
                  "Provider internals have their own reported caps; wall-time and memory caps are enforced externally.",
              ]}
    return _plain(result)
