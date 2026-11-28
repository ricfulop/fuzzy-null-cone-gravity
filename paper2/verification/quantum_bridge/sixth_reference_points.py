"""Direct pointwise bridge for the original mixed sixth-order recursion.

UnreducedPointOracle consumes supplied exact local_geometry output and
physical jets of independently computed dependent fields. It constructs no
geometry or inverse, imports no production Green/pairing backend, and never
substitutes target values for numerical fields. Import is inert.

The raw q,s polynomials multiply exp(-r^2/2). Geometry's integer rate n
means beta=n/2. Metric/inverse/Gamma and R[n] have rate n; b[n] has rate
n+1. R is pointwise curvature, while L[6] is the integrated curvature density.
All arrays use derivative-first spatial indices and ordinary Taylor orders.

The exact local data remain available through the deeply read-only .local
mapping for separate origin/exterior source construction. Neither sparse
point values nor this module supplies a source remainder bound. All work
occurs on explicitly requested construction or evaluation in a tracked run.
The point/term counters are declared work proxies, not CPU or memory bounds.
No numerical validation of this new bridge has yet been performed.
"""

from collections.abc import Mapping
from fractions import Fraction
from math import exp, fsum, isfinite
from types import MappingProxyType

from reduced_gaunt_exact import assert_tt
from reduced_gaunt_vertices import pg
from sixth_unreduced_algebra import third_sources, fourth_sources, energy_densities
from york_operator_fields import york_series


def _positive_int(value, name):
    if type(value) is not int or value <= 0:
        raise ValueError(name+" must be a positive integer")
    return value


def _finite(value, name):
    if isinstance(value, bool):
        raise TypeError(name+" must be a finite real number")
    result = float(value)
    if not isfinite(result):
        raise ValueError(name+" must be finite")
    return result


def _array(value, shape, name):
    if not shape:
        return _finite(value, name)
    if not isinstance(value, (tuple, list)) or len(value) != shape[0]:
        raise ValueError(name+" has an incorrect Cartesian shape")
    return tuple(_array(item, shape[1:], name) for item in value)


def _freeze(value):
    if isinstance(value, Mapping):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, (tuple, list)):
        return tuple(_freeze(item) for item in value)
    return value


def _indices(layout):
    if layout is None:
        return set()
    if isinstance(layout, Mapping):
        return set().union(*(_indices(item) for item in layout.values()))
    if isinstance(layout, tuple):
        return set().union(*(_indices(item) for item in layout))
    return {layout}


def _restore(layout, values):
    if layout is None:
        return None
    if isinstance(layout, Mapping):
        return {key: _restore(item, values) for key, item in layout.items()}
    if isinstance(layout, tuple):
        return tuple(_restore(item, values) for item in layout)
    return values[layout]


class UnreducedPointOracle:
    """Compile supplied exact geometry once; evaluate only supplied field jets.

    Required geometry: metric/inverse/Gamma through two, pointwise R2--R4,
    density L6, and b1--b3 from the same q,s. The caller retains independent
    geometry anchors. Compilation preserves exact polynomial/rate keys and
    deduplicates identical physical polynomials before floating evaluation.
    """

    def __init__(self, q, s, geometry, *, max_terms=1_000_000,
                 max_degree=64, max_compile_products=20_000_000,
                 max_point_evaluations=1_000_000,
                 max_term_evaluations=200_000_000, progress=None):
        self._caps = {name: _positive_int(value, name) for name, value in (
            ("terms", max_terms), ("degree", max_degree),
            ("compile_products", max_compile_products),
            ("point_evaluations", max_point_evaluations),
            ("term_evaluations", max_term_evaluations))}
        self._work = {"compile_products": 0, "stored_terms": 0,
                      "point_evaluations": 0, "term_evaluations": 0}
        self._compiled, self._lookup = [], {}
        if not isinstance(geometry, Mapping):
            raise TypeError("geometry must be supplied exact local data")
        if _freeze(geometry.get("q")) != _freeze(q):
            raise ValueError("geometry belongs to a different q")
        if not geometry.get("E") or _freeze(geometry["E"][0]) != _freeze(s):
            raise ValueError("geometry tangent does not describe supplied s")
        assert_tt(q)
        assert_tt(s)
        for family in ("metric", "inverse", "gamma"):
            if len(geometry.get(family, ())) < 3:
                raise ValueError(family+" coefficients through order two are required")
        for family, orders in (("R", (2, 3, 4)), ("L", (6,)), ("b", (1, 2, 3))):
            for n in orders:
                try:
                    geometry[family][n]
                except (KeyError, IndexError, TypeError) as exc:
                    raise ValueError("missing "+family+str(n)) from exc
        rates = geometry["gaussian_rates"]
        for family, orders in (("metric", range(3)), ("inverse", range(3)),
                               ("gamma", (1, 2)), ("R", (2, 3, 4)),
                               ("L", (6,)), ("b", (1, 2, 3))):
            for n in orders:
                expected = n+1 if family == "b" else n
                if rates[family][n] != expected:
                    raise ValueError("inconsistent physical Gaussian rate for "+family+str(n))
        self._local = _freeze({"q": q, "s": s, "geometry": geometry,
            "gaussian_rates": rates, "gaussian_rate_convention": "integer rate n means beta=n/2",
            "physical_betas": {"q": Fraction(1, 2), "s": Fraction(1, 2),
                "R": {n: Fraction(n, 2) for n in (2, 3, 4)},
                "L": {6: Fraction(3)}, "b": {n: Fraction(n+1, 2) for n in (1, 2, 3)}}})
        metric, dmetric, ddmetric, inverse, dinverse, gamma = [], [], [], [], [], [None]
        for n in range(3):
            g, gi = geometry["metric"][n], geometry["inverse"][n]
            dg = tuple(tuple(tuple(self._differentiate(g[i][j], a, n)
                        for j in range(3)) for i in range(3)) for a in range(3))
            ddg = tuple(tuple(tuple(tuple(self._differentiate(dg[a][i][j], b, n)
                         for j in range(3)) for i in range(3)) for b in range(3)) for a in range(3))
            dgi = tuple(tuple(tuple(self._differentiate(gi[i][j], a, n)
                         for j in range(3)) for i in range(3)) for a in range(3))
            metric.append(self._tree(g, (3, 3), n))
            dmetric.append(self._tree(dg, (3, 3, 3), n))
            ddmetric.append(self._tree(ddg, (3, 3, 3, 3), n))
            inverse.append(self._tree(gi, (3, 3), n))
            dinverse.append(self._tree(dgi, (3, 3, 3), n))
            if n:
                gamma.append(self._tree(geometry["gamma"][n], (3, 3, 3), n))
        self._layouts = {
            "q": self._tree(q, (3, 3), 1), "s": self._tree(s, (3, 3), 1),
            "metric": tuple(metric), "dmetric": tuple(dmetric), "ddmetric": tuple(ddmetric),
            "inverse": tuple(inverse), "dinverse": tuple(dinverse), "gamma": tuple(gamma),
            "R": {n: self._tree(geometry["R"][n], (), n) for n in (2, 3, 4)},
            "b": {n: self._tree(geometry["b"][n], (3,), n+1) for n in (1, 2, 3)},
            "L6": self._tree(geometry["L"][6], (), 6)}
        if progress is not None:
            progress({"stage": "unreduced_point_oracle_prepared", "statistics": self.statistics()})

    @property
    def local(self):
        """Deeply read-only exact source data; no fitted or computed Green fields."""
        return self._local

    def _charge_compile(self, count):
        if self._work["compile_products"]+count > self._caps["compile_products"]:
            raise ValueError("max_compile_products exceeded")
        self._work["compile_products"] += count

    def _differentiate(self, poly, axis, rate):
        self._charge_compile(2*len(poly))
        return pg(poly, axis, rate)

    def _tree(self, tree, shape, rate):
        if shape:
            if not isinstance(tree, (tuple, list)) or len(tree) != shape[0]:
                raise ValueError("local polynomial tensor has an incorrect shape")
            return tuple(self._tree(item, shape[1:], rate) for item in tree)
        if not isinstance(tree, Mapping):
            raise TypeError("local scalar must be an exact Cartesian polynomial")
        self._charge_compile(len(tree))
        terms = []
        for powers, coefficient in sorted(tree.items()):
            if (not isinstance(powers, tuple) or len(powers) != 3
                    or any(type(p) is not int or p < 0 for p in powers)
                    or sum(powers) > self._caps["degree"]):
                raise ValueError("invalid or oversized Cartesian monomial")
            if isinstance(coefficient, bool) or not isinstance(coefficient, (int, Fraction)):
                raise TypeError("local coefficients must remain exact rational numbers")
            if coefficient:
                terms.append((powers, Fraction(coefficient)))
        key = (rate, tuple(terms))
        if key not in self._lookup:
            if self._work["stored_terms"]+len(terms) > self._caps["terms"]:
                raise ValueError("max_terms exceeded")
            converted = tuple((powers, _finite(c, "coefficient")) for powers, c in terms)
            if any(c == 0 for _, c in converted):
                raise ValueError("nonzero rational coefficient underflows binary64")
            self._work["stored_terms"] += len(terms)
            self._lookup[key] = len(self._compiled)
            self._compiled.append((rate, converted))
        return self._lookup[key]

    def statistics(self):
        return {"caps": dict(self._caps), "work": dict(self._work),
                "unique_polynomials": len(self._compiled),
                "coefficient_arithmetic": "exact compilation, binary64 physical evaluation",
                "work_units": "coefficient visits/derivative upper bounds and evaluated monomial terms",
                "geometry_recomputed": False, "inverse_fields_constructed": False,
                "accuracy_certified": False}

    def local_point(self, point, *, order=3):
        """Evaluate physical local data for third, fourth, or energy (6) stage."""
        if type(order) is not int or order not in (3, 4, 6):
            raise ValueError("order must be 3, 4 or 6")
        point = _array(point, (3,), "point")
        radius_squared = _finite(fsum(x*x for x in point), "radius squared")
        nmax = 1 if order == 3 else 2
        layouts = {name: self._layouts[name] for name in ("q", "s")}
        layouts.update({name: self._layouts[name][:nmax+1]
                        for name in ("metric", "dmetric", "ddmetric", "inverse", "dinverse", "gamma")})
        layouts["R"] = {n: self._layouts["R"][n] for n in range(2, 4 if order == 3 else 5)}
        layouts["b"] = {n: self._layouts["b"][n] for n in range(1, 4 if order == 6 else 3)}
        if order == 6:
            layouts["L6"] = self._layouts["L6"]
        indices = sorted(_indices(layouts))
        term_count = sum(len(self._compiled[i][1]) for i in indices)
        if self._work["point_evaluations"]+1 > self._caps["point_evaluations"]:
            raise ValueError("max_point_evaluations exceeded")
        if self._work["term_evaluations"]+term_count > self._caps["term_evaluations"]:
            raise ValueError("max_term_evaluations exceeded")
        self._work["point_evaluations"] += 1
        self._work["term_evaluations"] += term_count
        maxima = tuple(max((p[a] for i in indices for p, _ in self._compiled[i][1]), default=0)
                       for a in range(3))
        powers = []
        for x, maximum in zip(point, maxima):
            row = [1.0]
            for _ in range(maximum):
                row.append(_finite(row[-1]*x, "Cartesian power"))
            powers.append(row)
        envelopes = {rate: exp(-rate*radius_squared/2) for rate, _ in
                     (self._compiled[i] for i in indices)}
        values = {}
        for index in indices:
            rate, terms = self._compiled[index]
            if terms and envelopes[rate] == 0:
                raise ValueError("physical Gaussian envelope underflows; use exterior bounds instead")
            value = fsum(c*powers[0][p[0]]*powers[1][p[1]]*powers[2][p[2]] for p, c in terms)
            values[index] = _finite(envelopes[rate]*value, "physical polynomial")
        return {"point": point, "order": order, **_restore(layouts, values)}

    @staticmethod
    def _york(local, w1_jet, order):
        return york_series(local["metric"], local["dmetric"], local["ddmetric"],
                           local["inverse"], local["gamma"], w1_jet["value"],
                           w1_jet["gradient"], w1_jet["hessian"], max_order=order)

    @staticmethod
    def _third(local, u_jet, york):
        return third_sources(q=local["q"], s=local["s"], t0=york["L"][0],
            r2=local["R"][2], r3=local["R"][3], b2=local["b"][2], a1w1=york["A"][1],
            u_jet=u_jet, inverse1=local["inverse"][1], dinverse1=local["dinverse"][1])

    @staticmethod
    def _fourth(local, u_jet, u3_jet, w2_jet, ec, york):
        dw = _array(w2_jet["gradient"], (3, 3), "W2 gradient")
        divergence = fsum(dw[k][k] for k in range(3))
        l0w2 = tuple(tuple(dw[j][i]+dw[i][j]-(2*divergence/3 if i == j else 0)
                          for j in range(3)) for i in range(3))
        result = fourth_sources(q=local["q"], s=local["s"], t0=york["L"][0],
            t1=york["L"][1], l0w2=l0w2, ec=ec, r2=local["R"][2], r4=local["R"][4],
            u_jet=u_jet, u3_jet=u3_jet, inverse1=local["inverse"][1],
            dinverse1=local["dinverse"][1], inverse2=local["inverse"][2],
            dinverse2=local["dinverse"][2])
        return {**result, "l0w2": l0w2}

    def third_point(self, point, *, u_jet, w1_jet):
        """Return F3/r and raw tangent data from supplied physical u/W1 jets."""
        local = self.local_point(point, order=3)
        york = self._york(local, w1_jet, 1)
        return {**self._third(local, u_jet, york), "t0": york["L"][0],
                "t1": york["L"][1], "a1w1": york["A"][1]}

    def fourth_point(self, point, *, u_jet, u3_jet, w1_jet, w2_jet, ec):
        """Return the original F4/s4, with explicit W2/u3 and supplied Ec."""
        local = self.local_point(point, order=4)
        york = self._york(local, w1_jet, 2)
        return {**self._third(local, u_jet, york),
                **self._fourth(local, u_jet, u3_jet, w2_jet, ec, york),
                "t0": york["L"][0], "t1": york["L"][1], "t2": york["L"][2],
                "a1w1": york["A"][1], "a2w1": york["A"][2]}

    def energy_point(self, point, *, u_jet, u3_jet, u4_jet, w1_jet, w2_jet, ec):
        """Nine original spatial densities and fields for separate flux checks.

        No canonical c-vector term is added. The global Gram square remains
        with assemble_integrated_energy, once after spatial integration.
        """
        local = self.local_point(point, order=6)
        york = self._york(local, w1_jet, 2)
        sources = {**self._third(local, u_jet, york),
                   **self._fourth(local, u_jet, u3_jet, w2_jet, ec, york)}
        densities = energy_densities(q=local["q"], s=local["s"], sources=sources,
            w1=w1_jet["value"], w2=w2_jet["value"], b3=local["b"][3],
            a2w1=york["A"][2], l6=local["L6"], u=u_jet["value"],
            u3=u3_jet["value"], u4=u4_jet["value"])
        return {"densities": densities, "sources": sources, "local": local,
                "york": york, "global_gram_square_included": False,
                "spatial_or_exterior_integration_performed": False}
