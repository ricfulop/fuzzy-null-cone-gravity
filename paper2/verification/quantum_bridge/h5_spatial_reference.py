"""Independent Cartesian/spherical quadrature of the full reduced H5 density.

Input is the dictionary returned by reduced_fifth_order.prepare_fields.
No geometry, source inversion, or amplitude construction is performed here.
Existing polynomial and HeatField jets are differentiated as needed, compiled
to binary64 once, and evaluated at Cartesian points on complete spherical
quadrature grids. Heat moments are shared by every field and angle at a fixed
radius. No axial symmetry, angular projection, or exact angular-pairing
backend is used in this reference.

integrate_h5_reference(data, radial_orders=(12,16,20),
    angular_orders=((8,16),(10,24),(12,32)), radius=24, panel_width=2,
    progress=None)

Pairs radial_orders[k] with angular_orders[k]=(n_mu,n_phi). Each level returns
eight signed entries (six H5 families, with both York and scalar-gradient
alternatives), total_strong, total_weak, total (=strong), york_forms,
scalar_gradient_forms, and successive-level
changes. There is no automatic convergence verdict or computed tail bound.
The caller supplies acceptance criteria and an exterior bound from its
independent pairing backend. Strong/weak disagreement at finite radius may
include a boundary term as well as numerical error. Jointly refined levels
do not separately identify radial and angular discretization errors.

The integrated strong density is
 -L5 -4 b2.W +2 W.A1W +R3*(v-u)
 +2*(S:T0)*(7*u+v) -8*grad(v).Q.grad(u).
The alternative York term is -2*T0:T1 +2*Tr(Q*T0^2). W is held
contravariant and fixed in york_series; both connection terms are retained.
The independent scalar-gradient strong term is +8*v*Q:Hess(u), since the
physical Q is flat TT. Both reported totals retain the weak scalar-gradient
term; their strong/weak names distinguish ONLY the York alternatives.
All q^5, q^3*p^2, q*p^4 terms are included for arbitrary supplied TT directions.

This is a finite-radius floating-point spatial reference for a homogeneous
Hamiltonian coefficient, not a four-point loop or a full scattering amplitude.
No construction, evaluation, residual checks, or runs occur on import.
"""

from functools import lru_cache
from math import ceil, cos, exp, fsum, isfinite, pi, sin, sqrt

from gaussian_green_fields import heat_integrals
from reduced_gaunt_vertices import pg
from run import gauss_nodes
from york_operator_fields import york_series


FAMILIES = ("curvature", "york_source", "york_strong", "york_weak",
            "scalar_curvature", "scalar_momentum", "scalar_gradient", "scalar_gradient_strong")


@lru_cache(None)
def _nodes(order):
    # run.gauss_nodes has NORMALIZED weights summing to one, not two.
    return tuple((float(x), float(w)) for x, w in gauss_nodes(order, 60))


class _CompiledFields:
    """Float coefficients and shared Cartesian monomials, built once per call."""

    def __init__(self):
        self.powers = []
        self.power_ids = {}
        self.specs = []
        self.poly_cache = {}
        self.heat_cache = {}
        self.rates = {}

    def _power(self, powers):
        if powers not in self.power_ids:
            self.power_ids[powers] = len(self.powers)
            self.powers.append(powers)
        return self.power_ids[powers]

    def polynomial(self, poly, rate):
        # rate is the full Gaussian rate: physical field = P exp(-rate*r^2/2).
        rate = float(rate) if poly else 0.0
        key = (rate, tuple(sorted(poly.items())))
        if key not in self.poly_cache:
            terms = tuple((self._power(powers), float(coefficient))
                          for powers, coefficient in key[1] if coefficient)
            self.poly_cache[key] = len(self.specs)
            self.specs.append(("polynomial", rate, terms))
        return self.poly_cache[key]

    def heat(self, field):
        key = (field.a, tuple(sorted(field.terms.items())))
        if key not in self.heat_cache:
            grouped = {}
            for (powers, m), coefficient in key[1]:
                grouped.setdefault(self._power(powers), []).append((m, float(coefficient)))
            terms = tuple((power_id, tuple(coefficients)) for power_id, coefficients in grouped.items())
            self.heat_cache[key] = len(self.specs)
            self.specs.append(("heat", field.a, terms))
            old_max = self.rates.get(field.a, (float(field.a), 0))[1]
            self.rates[field.a] = (float(field.a), max(old_max, field.max_m))
        return self.heat_cache[key]

    def angular_grid(self, n_mu, n_phi):
        maxima = [max((powers[i] for powers in self.powers), default=0) for i in range(3)]
        azimuth = [(cos(2 * pi * j / n_phi), sin(2 * pi * j / n_phi)) for j in range(n_phi)]
        result = []
        for mu, weight in _nodes(n_mu):
            transverse = sqrt(max(0.0, 1 - mu * mu))
            for cp, sp in azimuth:
                coordinates = (transverse * cp, transverse * sp, mu)
                tables = []
                for value, maximum in zip(coordinates, maxima):
                    row = [1.0]
                    for _ in range(maximum):
                        row.append(row[-1] * value)
                    tables.append(row)
                monomials = tuple(tables[0][p[0]] * tables[1][p[1]] * tables[2][p[2]]
                                  for p in self.powers)
                # Full sphere integral: 4*pi times normalized mu weight/n_phi.
                result.append((4 * pi * weight / n_phi, monomials))
        return result

    def at_radius(self, radius):
        """Shared heat moments and radial monomial factors for all angular nodes."""
        r2 = radius * radius
        moments = {key: heat_integrals(value * r2, maximum)
                   for key, (value, maximum) in self.rates.items()}
        maximum_degree = max((sum(powers) for powers in self.powers), default=0)
        radial_powers = [1.0]
        for _ in range(maximum_degree):
            radial_powers.append(radial_powers[-1] * radius)
        radial = [radial_powers[sum(powers)] for powers in self.powers]
        envelopes = {rate: exp(-rate * r2 / 2)
                     for kind, rate, _ in self.specs if kind == "polynomial"}
        compiled = []
        for kind, rate, terms in self.specs:
            if kind == "polynomial":
                compiled.append(tuple((index, c * radial[index] * envelopes[rate]) for index, c in terms))
            else:
                table = moments[rate]
                compiled.append(tuple((index, radial[index] * fsum(c * table[m] for m, c in coeffs))
                                      for index, coeffs in terms))
        return compiled

    @staticmethod
    def evaluate(compiled, monomials):
        return [fsum(c * monomials[index] for index, c in terms) for terms in compiled]

    def statistics(self):
        return {"polynomial_fields": len(self.poly_cache), "heat_fields": len(self.heat_cache),
                "shared_cartesian_monomials": len(self.powers),
                "heat_moment_orders": {str(a): maximum for a, (_, maximum) in self.rates.items()}}


def _map_tree(tree, depth, function):
    if depth == 0:
        return function(tree)
    return [_map_tree(item, depth - 1, function) for item in tree]


def _values(refs, values):
    if refs is None:
        return None
    if isinstance(refs, int):
        return values[refs]
    return [_values(item, values) for item in refs]


def _compile(data):
    """Differentiate supplied fields only; do not call prepare_fields or an inverse."""
    geometry = data["geometry"]
    if geometry.get("q") != data["q"]:
        raise ValueError("prepared geometry and q do not agree")
    if 3 not in geometry["R"] or 5 not in geometry["L"] or len(data["york"]["b"]) < 3:
        raise ValueError("H5 reference needs R3, L5, and b2 in the supplied data")
    compiler = _CompiledFields()
    refs = {
        "q": _map_tree(data["q"], 2, lambda poly: compiler.polynomial(poly, 1)),
        "s": _map_tree(data["s"], 2, lambda poly: compiler.polynomial(poly, 1)),
        "r3": compiler.polynomial(geometry["R"][3], 3),
        "l5": compiler.polynomial(geometry["L"][5], 5),
        "b2": _map_tree(data["york"]["b"][2], 1, lambda poly: compiler.polynomial(poly, 3)),
        "u": compiler.heat(data["u"]), "v": compiler.heat(data["v"]),
        "grad_u": _map_tree(data["grad_u"], 1, compiler.heat),
        "grad_v": _map_tree(data["grad_v"], 1, compiler.heat),
        "w": _map_tree(data["w"], 1, compiler.heat),
        "dw": _map_tree(data["dw"], 2, compiler.heat),
    }
    ddw = [[[data["dw"][a][k].derivative(b) for k in range(3)] for b in range(3)] for a in range(3)]
    refs["ddw"] = _map_tree(ddw, 3, compiler.heat)
    # Six distinct analytic scalar Hessian fields; symmetry supplies the rest.
    refs["hess_u"] = [[None] * 3 for _ in range(3)]
    for a in range(3):
        for b in range(a, 3):
            ref = compiler.heat(data["grad_u"][a].derivative(b))
            refs["hess_u"][a][b] = refs["hess_u"][b][a] = ref
    refs.update({"metric": [], "dmetric": [], "ddmetric": [], "inverse": [], "gamma": [None]})
    for n in (0, 1):
        metric = geometry["metric"][n]
        derivative = [[[pg(metric[i][j], a, n) for j in range(3)] for i in range(3)] for a in range(3)]
        second = [[[[pg(derivative[a][i][j], b, n) for j in range(3)] for i in range(3)]
                   for b in range(3)] for a in range(3)]
        compile_poly = lambda poly: compiler.polynomial(poly, n)
        refs["metric"].append(_map_tree(metric, 2, compile_poly))
        refs["dmetric"].append(_map_tree(derivative, 3, compile_poly))
        refs["ddmetric"].append(_map_tree(second, 4, compile_poly))
        refs["inverse"].append(_map_tree(geometry["inverse"][n], 2, compile_poly))
    refs["gamma"].append(_map_tree(geometry["gamma"][1], 3, lambda poly: compiler.polynomial(poly, 1)))
    return compiler, refs


def _density(values, refs):
    fields = {name: _values(reference, values) for name, reference in refs.items()}
    york = york_series(fields["metric"], fields["dmetric"], fields["ddmetric"],
                       fields["inverse"], fields["gamma"], fields["w"], fields["dw"],
                       fields["ddw"], max_order=1)
    t0, t1 = york["L"]
    q, s, w = fields["q"], fields["s"], fields["w"]
    u, v = fields["u"], fields["v"]
    result = {
        "curvature": -fields["l5"],
        "york_source": -4 * fsum(fields["b2"][i] * w[i] for i in range(3)),
        "york_strong": 2 * fsum(w[i] * york["A"][1][i] for i in range(3)),
        "york_weak": -2 * fsum(t0[i][j] * t1[i][j] for i in range(3) for j in range(3))
                     + 2 * fsum(q[i][k] * t0[i][j] * t0[k][j]
                                for i in range(3) for j in range(3) for k in range(3)),
        "scalar_curvature": fields["r3"] * (v - u),
        "scalar_momentum": 2 * fsum(s[i][j] * t0[i][j] for i in range(3) for j in range(3)) * (7 * u + v),
        "scalar_gradient": -8 * fsum(fields["grad_v"][i] * q[i][j] * fields["grad_u"][j]
                                     for i in range(3) for j in range(3)),
        "scalar_gradient_strong": 8 * v * fsum(q[i][j] * fields["hess_u"][i][j]
                                               for i in range(3) for j in range(3)),
    }
    if not all(isfinite(value) for value in result.values()):
        raise ArithmeticError("nonfinite H5 spatial density")
    return result


def _change(old, new):
    delta = new - old
    return {"signed_delta": delta, "absolute_delta": abs(delta),
            "relative_delta": abs(delta) / abs(new) if new else None,
            "scaled_delta": abs(delta) / max(1.0, abs(new))}


def integrate_h5_reference(data, radial_orders=(12, 16, 20),
                           angular_orders=((8, 16), (10, 24), (12, 32)),
                           radius=24, panel_width=2, progress=None):
    """Return JSON-serializable finite-radius strong/weak H5 reference levels.

    progress(dict) is called after every completed radial panel. It reports
    one-based level/panel indices and cumulative quadrature point count in
    that level. Fractions are converted once during field compilation, not
    at individual quadrature points. The input dictionary is not modified.
    """
    radial_orders = tuple(radial_orders)
    angular_orders = tuple(tuple(pair) for pair in angular_orders)
    if not radial_orders or len(radial_orders) != len(angular_orders):
        raise ValueError("radial_orders and angular_orders must have the same nonzero length")
    if any(type(order) is not int or order < 2 for order in radial_orders):
        raise ValueError("radial orders must be integers at least two")
    if any(len(pair) != 2 or any(type(n) is not int or n < 2 for n in pair) for pair in angular_orders):
        raise ValueError("angular orders must be (n_mu,n_phi) integer pairs, each at least two")
    radius, panel_width = float(radius), float(panel_width)
    if not all(isfinite(value) and value > 0 for value in (radius, panel_width)):
        raise ValueError("radius and panel_width must be finite and positive")
    compiler, refs = _compile(data)
    panels = ceil(radius / panel_width)
    levels = []
    for level, (radial_order, (n_mu, n_phi)) in enumerate(zip(radial_orders, angular_orders), 1):
        angular = compiler.angular_grid(n_mu, n_phi)
        panel_values = {name: [] for name in FAMILIES}
        absolute_values = {name: [] for name in FAMILIES}
        points = 0
        for panel in range(panels):
            left, right = panel * panel_width, min(radius, (panel + 1) * panel_width)
            length, midpoint = right - left, (right + left) / 2
            radial_values = {name: [] for name in FAMILIES}
            radial_absolute = {name: [] for name in FAMILIES}
            for node, weight in _nodes(radial_order):
                r = midpoint + length * node / 2
                compiled = compiler.at_radius(r)
                angular_values = {name: [] for name in FAMILIES}
                angular_absolute = {name: [] for name in FAMILIES}
                for angular_weight, monomials in angular:
                    density = _density(compiler.evaluate(compiled, monomials), refs)
                    for name, value in density.items():
                        angular_values[name].append(angular_weight * value)
                        angular_absolute[name].append(angular_weight * abs(value))
                # Normalized GL weights require length, not length/2.
                radial_weight = length * weight * r * r
                for name in FAMILIES:
                    radial_values[name].append(radial_weight * fsum(angular_values[name]))
                    radial_absolute[name].append(radial_weight * fsum(angular_absolute[name]))
                points += len(angular)
            for name in FAMILIES:
                panel_values[name].append(fsum(radial_values[name]))
                absolute_values[name].append(fsum(radial_absolute[name]))
            if progress is not None:
                progress({"stage": "h5_spatial_reference", "level": level, "levels": len(radial_orders),
                          "panel": panel + 1, "panels": panels, "radial_order": radial_order,
                          "n_mu": n_mu, "n_phi": n_phi, "radius_completed": right,
                          "points_completed": points})
        families = {name: fsum(panel_values[name]) for name in FAMILIES}
        absolute_integrals = {name: fsum(absolute_values[name]) for name in FAMILIES}
        common = [families[name] for name in FAMILIES
                  if name not in ("york_strong", "york_weak", "scalar_gradient_strong")]
        strong = fsum(common + [families["york_strong"]])
        weak = fsum(common + [families["york_weak"]])
        difference = families["york_strong"] - families["york_weak"]
        scalar_difference = families["scalar_gradient_strong"] - families["scalar_gradient"]
        row = {"level": level, "radial_order": radial_order, "angular_orders": [n_mu, n_phi],
               "radius": radius, "panel_width": panel_width, "panels": panels, "points": points,
               "families": families, "total_strong": strong, "total_weak": weak, "total": strong,
               "york_forms": {"strong": families["york_strong"], "weak": families["york_weak"],
                              "difference": difference, "absolute_difference": abs(difference)},
               "scalar_gradient_forms": {"strong": families["scalar_gradient_strong"],
                                         "weak": families["scalar_gradient"],
                                         "difference": scalar_difference,
                                         "absolute_difference": abs(scalar_difference)},
               "absolute_density_integrals": absolute_integrals,
               "cancellation_ratios": {name: absolute_integrals[name] / abs(families[name])
                                       if families[name] else None for name in FAMILIES},
               "refinement_from_previous": None}
        if levels:
            previous = levels[-1]
            row["refinement_from_previous"] = {
                "families": {name: _change(previous["families"][name], families[name]) for name in FAMILIES},
                "total_strong": _change(previous["total_strong"], strong),
                "total_weak": _change(previous["total_weak"], weak),
            }
        levels.append(row)
    return {
        "scope": "Numerical finite-radius full H5 spatial reference; not a loop amplitude.",
        "config": {"radial_orders": list(radial_orders), "angular_orders": [list(pair) for pair in angular_orders],
                   "radius": radius, "panel_width": panel_width},
        "compilation": compiler.statistics(), "levels": levels, "latest": levels[-1],
        "tail_bound": None, "tail_bound_source": "caller supplies independent exterior bound",
        "convergence_decision": "caller supplied criteria; refinement differences are not certified errors",
        "arithmetic": "binary64; compensated summation, with once-per-radius shared heat moments",
    }
