"""Inline wiring controls for the unreduced point oracle; inert on import.

The caller supplies an already prepared oracle in a tracked calculation.
Direct Decimal polynomial/product-rule jets do not use pg. The York reference
lowers W first and differentiates its covariant derivative, independently of
the Lie-derivative construction used by york_series. Original source-algebra
functions remain shared: their independent rational controls are a separate
required ingredient, not replaced by these wiring checks. No inverse, geometry
construction, quadrature or gravitational energy evaluation occurs here.
"""

from collections.abc import Mapping
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from math import isfinite

from sixth_reference_points import UnreducedPointOracle
from sixth_unreduced_algebra import third_sources, fourth_sources, energy_densities


def _decimal(value):
    if isinstance(value, D):
        return value
    if isinstance(value, float):
        return D.from_float(value)
    value = F(value)
    return D(value.numerator)/D(value.denominator)


def _map(tree, function):
    if isinstance(tree, Mapping):
        return {key: _map(value, function) for key, value in tree.items()}
    if isinstance(tree, (tuple, list)):
        return tuple(_map(value, function) for value in tree)
    return function(tree)


def _zero(rank):
    return D(0) if rank == 0 else tuple(_zero(rank-1) for _ in range(3))


def _flatten(value, prefix=""):
    if isinstance(value, Mapping):
        result = {}
        for key, item in value.items():
            result.update(_flatten(item, prefix+"."+str(key)))
        return result
    if isinstance(value, (tuple, list)):
        result = {}
        for i, item in enumerate(value):
            result.update(_flatten(item, prefix+"."+str(i)))
        return result
    return {prefix: value}


class _Direct:
    def __init__(self, point, work, limit):
        self.point = tuple(_decimal(x) for x in point)
        self.work, self.limit, self.cache = work, limit, {}

    def polynomial(self, polynomial, axes=()):
        """Monomial falling factorials, without a derivative-polynomial helper."""
        if self.work["reference_monomials"]+len(polynomial) > self.limit:
            raise ValueError("point_control_max_reference_monomials exceeded")
        self.work["reference_monomials"] += len(polynomial)
        orders = tuple(axes.count(i) for i in range(3))
        total = D(0)
        for powers, coefficient in polynomial.items():
            if any(p < d for p, d in zip(powers, orders)):
                continue
            term = _decimal(coefficient)
            for x, power, count in zip(self.point, powers, orders):
                for j in range(count):
                    term *= power-j
                exponent = power-count
                term *= D(1) if exponent == 0 else x**exponent
            total += term
        return total

    def scalar(self, polynomial, rate, axes=()):
        key = (rate, tuple(sorted(polynomial.items())), tuple(axes))
        if key in self.cache:
            return self.cache[key]
        n = D(rate)
        x = self.point
        envelope = (-n*sum(a*a for a in x)/2).exp()
        p = self.polynomial(polynomial)
        if not axes:
            result = envelope*p
        elif len(axes) == 1:
            a = axes[0]
            result = envelope*(self.polynomial(polynomial, (a,))-n*x[a]*p)
        else:
            a, b = axes
            result = envelope*(self.polynomial(polynomial, (a, b))
                -n*x[a]*self.polynomial(polynomial, (b,))
                -n*x[b]*self.polynomial(polynomial, (a,))
                +(n*n*x[a]*x[b]-(n if a == b else 0))*p)
        self.cache[key] = result
        return result

    def tensor(self, tree, rank, rate, axes=()):
        if rank == 0:
            return self.scalar(tree, rate, axes)
        return tuple(self.tensor(item, rank-1, rate, axes) for item in tree)


def _direct_local(oracle, point, work, limit):
    raw = oracle.local
    geometry = raw["geometry"]
    direct = _Direct(point, work, limit)
    result = {"q": direct.tensor(raw["q"], 2, 1),
              "s": direct.tensor(raw["s"], 2, 1)}
    for name in ("metric", "inverse"):
        result[name] = tuple(direct.tensor(geometry[name][n], 2, n) for n in range(3))
    for name, original in (("dmetric", "metric"), ("dinverse", "inverse")):
        result[name] = tuple(tuple(direct.tensor(geometry[original][n], 2, n, (a,))
                            for a in range(3)) for n in range(3))
    result["ddmetric"] = tuple(tuple(tuple(direct.tensor(geometry["metric"][n], 2, n, (a, b))
                                for b in range(3)) for a in range(3)) for n in range(3))
    result["R"] = {n: direct.scalar(geometry["R"][n], n) for n in (2, 3, 4)}
    result["b"] = {n: direct.tensor(geometry["b"][n], 1, n+1) for n in (1, 2, 3)}
    result["L6"] = direct.scalar(geometry["L"][6], 6)
    # Build connections from direct metric jets, not supplied Gamma polynomials.
    dg, ddg, gi, dgi = (result[key] for key in ("dmetric", "ddmetric", "inverse", "dinverse"))
    result["gamma"] = tuple(tuple(tuple(tuple(sum(
        gi[b][k][t]*(dg[n-b][i][t][j]+dg[n-b][j][t][i]-dg[n-b][t][i][j])/2
        for b in range(n+1) for t in range(3)) for j in range(3))
        for i in range(3)) for k in range(3)) for n in range(3))
    result["dgamma"] = tuple(tuple(tuple(tuple(tuple(sum(
        (dgi[b][a][k][t]*(dg[n-b][i][t][j]+dg[n-b][j][t][i]-dg[n-b][t][i][j])
         +gi[b][k][t]*(ddg[n-b][a][i][t][j]+ddg[n-b][a][j][t][i]-ddg[n-b][a][t][i][j]))/2
        for b in range(n+1) for t in range(3)) for j in range(3))
        for i in range(3)) for k in range(3)) for a in range(3)) for n in range(3))
    return result


def _covariant_york(local, jet):
    """Formal ordinary coefficients using lowered-vector covariant derivatives."""
    g, dg, ddg, inv, gamma, dgamma = (local[k] for k in
        ("metric", "dmetric", "ddmetric", "inverse", "gamma", "dgamma"))
    w, dw, ddw = (_map(jet[k], _decimal) for k in ("value", "gradient", "hessian"))
    cov = tuple(tuple(sum(g[n][i][k]*w[k] for k in range(3)) for i in range(3)) for n in range(3))
    dcov = tuple(tuple(tuple(sum(dg[n][a][i][k]*w[k]+g[n][i][k]*dw[a][k]
                      for k in range(3)) for i in range(3)) for a in range(3)) for n in range(3))
    ddcov = tuple(tuple(tuple(tuple(sum(ddg[n][a][b][i][k]*w[k]+dg[n][a][i][k]*dw[b][k]
        +dg[n][b][i][k]*dw[a][k]+g[n][i][k]*ddw[a][b][k]
        for k in range(3)) for i in range(3)) for b in range(3)) for a in range(3)) for n in range(3))
    nabla = tuple(tuple(tuple(dcov[n][i][j]-sum(gamma[b][k][i][j]*cov[n-b][k]
        for b in range(n+1) for k in range(3)) for j in range(3)) for i in range(3)) for n in range(3))
    dnabla = tuple(tuple(tuple(tuple(ddcov[n][a][i][j]-sum(
        dgamma[b][a][k][i][j]*cov[n-b][k]+gamma[b][k][i][j]*dcov[n-b][a][k]
        for b in range(n+1) for k in range(3)) for j in range(3)) for i in range(3))
        for a in range(3)) for n in range(3))
    div = sum(dw[k][k] for k in range(3))
    ddiv = tuple(sum(ddw[a][k][k] for k in range(3)) for a in range(3))
    tensors = tuple(tuple(tuple(nabla[n][i][j]+nabla[n][j][i]-D(2)*g[n][i][j]*div/3
        for j in range(3)) for i in range(3)) for n in range(3))
    derivatives = tuple(tuple(tuple(tuple(dnabla[n][a][i][j]+dnabla[n][a][j][i]
        -D(2)*(dg[n][a][i][j]*div+g[n][i][j]*ddiv[a])/3
        for j in range(3)) for i in range(3)) for a in range(3)) for n in range(3))
    applied = tuple(tuple(sum(inv[b][j][a]*(derivatives[n-b][a][i][j]-sum(
        gamma[c][k][a][i]*tensors[n-b-c][k][j]+gamma[c][k][a][j]*tensors[n-b-c][i][k]
        for c in range(n-b+1) for k in range(3)))
        for b in range(n+1) for j in range(3) for a in range(3)) for i in range(3)) for n in range(3))
    return {"L": tensors, "A": applied}


def _fixture_jets():
    scalar = {"value": .375, "gradient": (1., -.5, .25),
              "hessian": ((2., -.25, .5), (-.25, -1., .75), (.5, .75, .125))}
    vector = {"value": (.5, -.75, .25),
              "gradient": tuple(tuple((3*a-2*i+1)/8 for i in range(3)) for a in range(3)),
              "hessian": tuple(tuple(tuple((a+b+2)*(i+1)/16 for i in range(3))
                              for b in range(3)) for a in range(3))}
    return {"u_jet": scalar, "u3_jet": _map(scalar, lambda x: -2*x),
            "u4_jet": _map(scalar, lambda x: 3*x), "w1_jet": vector,
            "w2_jet": _map(vector, lambda x: -1.5*x),
            "ec": ((.25, -.5, .125), (-.5, -.125, .75), (.125, .75, -.125))}


def _homogeneous_oracle(oracle, qscale, sscale, caps):
    def polynomial_tree(tree, scale):
        return _map(tree, lambda c: c*scale)
    raw = oracle.local
    g = dict(raw["geometry"])
    q = polynomial_tree(raw["q"], qscale)
    s = polynomial_tree(raw["s"], sscale)
    g["q"] = q
    for name in ("metric", "inverse", "gamma", "R", "L", "E", "b"):
        value = g[name]
        def scale_item(n, item):
            return None if item is None else polynomial_tree(item, qscale**n*(sscale if name in ("E", "b") else 1))
        g[name] = ({n: scale_item(n, item) for n, item in value.items()} if isinstance(value, Mapping)
                   else tuple(scale_item(n, item) for n, item in enumerate(value)))
    return UnreducedPointOracle(q, s, g, **caps)


def run_reference_point_controls(oracle, config=None, progress=None):
    """Return full inline comparisons, with named exception checkpoints.

    This checks the supplied direction with fixed artificial field jets. It
    neither recomputes its geometry nor solves physical dependent fields.
    Exact homogeneous transforms test wiring, not new gravitational solutions.
    """
    rows, current = [], "settings"
    work = {"reference_monomials": 0}
    reversed_s = doubled = None
    try:
        config = {} if config is None else config
        if not isinstance(config, Mapping):
            raise TypeError("config must be a mapping")
        tolerance = float(config.get("point_control_tolerance", 2e-10))
        precision = config.get("point_control_precision", 80)
        limit = config.get("point_control_max_reference_monomials", 20_000_000)
        points = tuple(tuple(p) for p in config.get("point_control_points",
                       ((0., 0., 0.), (.2, -.15, .1), (.7, -.3, .2))))
        if not isfinite(tolerance) or tolerance <= 0 or type(precision) is not int or not 32 <= precision <= 128:
            raise ValueError("finite positive tolerance and precision32..128 required")
        if type(limit) is not int or limit <= 0 or not 1 <= len(points) <= 8:
            raise ValueError("positive reference work limit and one through eight points required")
        if any(len(p) != 3 or any(not isfinite(float(x)) or abs(x) > 16 for x in p) for p in points):
            raise ValueError("bounded finite Cartesian points required")

        def record(name, actual, reference, *, exact_integers=False, **metadata):
            nonlocal current
            current = name
            a, b = _flatten(actual), _flatten(reference)
            if not a or set(a) != set(b):
                raise ValueError("comparison shape/key mismatch")
            checks = {}
            for key in a:
                if exact_integers and (type(a[key]) is not int or type(b[key]) is not int):
                    raise TypeError("Work-plan counts must be exact integers")
                value, target = float(a[key]), float(b[key])
                if not isfinite(value) or not isfinite(target):
                    raise ValueError("nonfinite point-oracle comparison")
                residual = value-target
                error = abs(residual)/max(1., abs(target))
                checks[key] = {"actual": a[key] if exact_integers else value,
                               "reference": b[key] if exact_integers else target,
                               "residual": residual, "scaled_error": error,
                               "passed": a[key] == b[key] if exact_integers else error <= tolerance}
            row = {"name": name, "checks": checks, "maximum": max(x["scaled_error"] for x in checks.values()),
                   "passed": all(x["passed"] for x in checks.values()),
                   "acceptance": "exact integer equality" if exact_integers else "scaled floating comparison",
                   **metadata}
            rows.append(row)
            if progress:
                progress({"stage": "point_control_complete", "control_result": row})

        current = "homogeneous_oracle_preparation"
        base_caps = oracle.statistics()["caps"]
        caps = {"max_"+key: value for key, value in base_caps.items()}
        reversed_s = _homogeneous_oracle(oracle, F(1), F(-1), caps)
        doubled = _homogeneous_oracle(oracle, F(2), F(2), caps)
        inputs = _fixture_jets()
        reversed_inputs = {name: _map(value, lambda x: -x) if name in ("w1_jet", "w2_jet", "ec") else value
                           for name, value in inputs.items()}
        degrees = {"u_jet": 2, "w1_jet": 2, "u3_jet": 3, "w2_jet": 3, "ec": 3, "u4_jet": 4}
        scaled_inputs = {name: _map(value, lambda x, n=degrees[name]: x*2**n) for name, value in inputs.items()}
        with localcontext() as context:
            context.prec = precision
            for index, point in enumerate(points):
                measured_plans, expected_plans = {}, {}

                def planned_call(name, order, action):
                    before = oracle.statistics()["work"]
                    value = action()
                    after = oracle.statistics()["work"]
                    plan = oracle.evaluation_plan(order=order)
                    measured_plans[name] = {k: after[k]-before[k] for k in ("point_evaluations", "term_evaluations")}
                    expected_plans[name] = {k: plan[k] for k in measured_plans[name]}
                    return value

                current = "point_"+str(index)+"_direct_local"
                local = _direct_local(oracle, point, work, limit)
                actual = planned_call("local6", 6, lambda: oracle.local_point(point, order=6))
                reference = {key: value for key, value in local.items() if key != "dgamma"}
                reference["gamma"] = (None,)+reference["gamma"][1:]
                # Order-zero Gamma is the public None sentinel, not a scalar.
                selected = {key: actual[key] for key in reference if key != "gamma"}
                expected = {key: reference[key] for key in reference if key != "gamma"}
                selected["gamma"] = actual["gamma"][1:]
                expected["gamma"] = reference["gamma"][1:]
                record(current, selected, expected, point=point,
                       reference_method="Decimal direct monomial/product-rule jets; metric-derived connection")
                current = "point_"+str(index)+"_covariant_york"
                y = _covariant_york(local, inputs["w1_jet"])
                full = planned_call("energy6", 6, lambda: oracle.energy_point(point, **inputs))
                record(current, {key: full["york"][key] for key in ("L", "A")}, y,
                       point=point, reference_method="Independent lowered-vector covariant derivative series")
                # Shared source functions check adapter wiring only; algebra has separate controls.
                current = "point_"+str(index)+"_source_wiring"
                fl = _map(local, float)
                fy = _map(y, float)
                source3 = third_sources(q=fl["q"], s=fl["s"], t0=fy["L"][0], r2=fl["R"][2],
                    r3=fl["R"][3], b2=fl["b"][2], a1w1=fy["A"][1], u_jet=inputs["u_jet"],
                    inverse1=fl["inverse"][1], dinverse1=fl["dinverse"][1])
                dw = inputs["w2_jet"]["gradient"]
                div = sum(dw[k][k] for k in range(3))
                lw2 = tuple(tuple(dw[i][j]+dw[j][i]-(2*div/3 if i == j else 0)
                                  for j in range(3)) for i in range(3))
                source4 = fourth_sources(q=fl["q"], s=fl["s"], t0=fy["L"][0], t1=fy["L"][1],
                    l0w2=lw2, ec=inputs["ec"], r2=fl["R"][2], r4=fl["R"][4], u_jet=inputs["u_jet"],
                    u3_jet=inputs["u3_jet"], inverse1=fl["inverse"][1], dinverse1=fl["dinverse"][1],
                    inverse2=fl["inverse"][2], dinverse2=fl["dinverse"][2])
                third = planned_call("third3", 3, lambda: oracle.third_point(
                    point, u_jet=inputs["u_jet"], w1_jet=inputs["w1_jet"]))
                fourth = planned_call("fourth4", 4, lambda: oracle.fourth_point(
                    point, **{k: v for k, v in inputs.items() if k != "u4_jet"}))
                record(current, {
                    "third": {k: third[k] for k in source3},
                    "fourth": {k: fourth[k] for k in source4},
                    "energy_sources": {k: full["sources"][k] for k in (*source3, *source4)}},
                    {"third": source3, "fourth": source4, "energy_sources": {**source3, **source4}}, point=point)
                record("point_"+str(index)+"_work_plans", measured_plans, expected_plans,
                       exact_integers=True,
                       reference_method="Public structural plan versus charged point/term counter increments")
                current = "point_"+str(index)+"_energy_wiring"
                densities = energy_densities(q=fl["q"], s=fl["s"], sources={**source3, **source4},
                    w1=inputs["w1_jet"]["value"], w2=inputs["w2_jet"]["value"], b3=fl["b"][3],
                    a2w1=fy["A"][2], l6=fl["L6"], u=inputs["u_jet"]["value"],
                    u3=inputs["u3_jet"]["value"], u4=inputs["u4_jet"]["value"])
                record(current, full["densities"], densities, point=point)
                current = "point_"+str(index)+"_canonical_impulse"
                delta_ec = ((.5, .25, -.125), (.25, -.25, .375), (-.125, .375, -.25))
                changed = dict(inputs)
                changed["ec"] = tuple(tuple(inputs["ec"][i][j]+delta_ec[i][j]
                                            for j in range(3)) for i in range(3))
                shifted = oracle.energy_point(point, **changed)
                trace_sd = float(sum(local["s"][i][j]*_decimal(delta_ec[j][i])
                                     for i in range(3) for j in range(3)))
                impulse = {"s4": -2*trace_sd, "F4": -2*trace_sd, "B4": -14*trace_sd}
                density_impulse = {key: 0. for key in densities}
                density_impulse["scalar_B4_u"] = 14*inputs["u_jet"]["value"]*trace_sd
                record(current, {"sources": {k: shifted["sources"][k]-full["sources"][k] for k in impulse},
                                 "densities": {k: shifted["densities"][k]-full["densities"][k] for k in densities}},
                       {"sources": impulse, "densities": density_impulse}, point=point,
                       reference_method="Direct Ec impulse; no source or energy formula called for target")
                current = "point_"+str(index)+"_u4_impulse"
                changed = dict(inputs)
                changed["u4_jet"] = dict(inputs["u4_jet"])
                changed["u4_jet"]["value"] += .5
                shifted = oracle.energy_point(point, **changed)
                density_impulse = {key: 0. for key in densities}
                density_impulse["scalar_B2_u4"] = -.5*source3["B2"]
                record(current, {k: shifted["densities"][k]-full["densities"][k] for k in densities},
                       density_impulse, point=point, reference_method="Direct u4 impulse")
                record("point_"+str(index)+"_ledger_exclusions",
                       {"gram_square_included": int(full["global_gram_square_included"]),
                        "integration_performed": int(full["spatial_or_exterior_integration_performed"]),
                        "density_count": len(full["densities"])},
                       {"gram_square_included": 0, "integration_performed": 0, "density_count": 9})
                current = "point_"+str(index)+"_momentum_parity"
                flipped = reversed_s.energy_point(point, **reversed_inputs)
                record(current, {
                    "densities": flipped["densities"], "F3": flipped["sources"]["F3"],
                    "F4": flipped["sources"]["F4"], "s4": flipped["sources"]["s4"], "r": flipped["sources"]["r"]},
                    {"densities": full["densities"], "F3": full["sources"]["F3"],
                     "F4": full["sources"]["F4"], "s4": full["sources"]["s4"], "r": _map(full["sources"]["r"], lambda x: -x)})
                current = "point_"+str(index)+"_joint_degree"
                scaled = doubled.energy_point(point, **scaled_inputs)
                record(current, {
                    "densities": scaled["densities"], "F3": scaled["sources"]["F3"],
                    "F4": scaled["sources"]["F4"], "s4": scaled["sources"]["s4"], "r": scaled["sources"]["r"]},
                    {"densities": _map(full["densities"], lambda x: 64*x), "F3": 8*full["sources"]["F3"],
                     "F4": 16*full["sources"]["F4"], "s4": 16*full["sources"]["s4"], "r": _map(full["sources"]["r"], lambda x: 8*x)})
        return {"passed": all(row["passed"] for row in rows), "rows": rows,
                "check_count": sum(len(row["checks"]) for row in rows),
                "failure_count": sum(not check["passed"] for row in rows for check in row["checks"].values()),
                "maximum": max(row["maximum"] for row in rows), "work": work,
                "oracle_statistics": {"base": oracle.statistics(), "reversed_s": reversed_s.statistics(),
                                      "doubled": doubled.statistics()},
                "settings": {"points": points, "tolerance": tolerance, "precision": precision,
                             "max_reference_monomials": limit, "oracle_caps_per_instance": caps},
                "limitations": ["Artificial field jets; no inverses or spatial integration",
                    "Local PG and York targets independent; source/energy algebra shared with separately controlled functions",
                    "Homogeneous transformations test wiring, not physical solution accuracy"],
                "independent_mixed_functional_validation": False, "accuracy_certified": False}
    except Exception as exc:
        if progress:
            statistics = {name: instance.statistics() for name, instance in
                          (("base", oracle), ("reversed_s", reversed_s), ("doubled", doubled))
                          if instance is not None}
            progress({"stage": "point_control_exception", "control_result": {
                "name": current, "passed": False, "exception_type": type(exc).__name__,
                "exception": str(exc), "completed_rows": len(rows), "work": dict(work),
                "oracle_statistics": statistics}})
        raise
