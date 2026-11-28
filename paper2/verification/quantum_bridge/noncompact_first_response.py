"""Independent first inverses for the actual mixed gravitational direction.

Q=E(diag(-1/2,-1/2,1)), S=B(XZ), in raw exp(-r^2/2) units. The
independent all-m radial path receives local polynomial sources only. The
production first-HEAT fields supply comparison targets after preparation;
they do not supply any radial, origin, boundary, or source coefficient.

This checks u=C(s2-R2), v=C(R2+7s2), W1=Y b1 and spatial derivatives.
There is no W2/u3/u4 solve, mixed H6 value, or loop calculation here.
Import is inert; the tracked driver owns the complete shared WorkBudget.
"""

from collections.abc import Mapping
from fractions import Fraction as F
from math import isfinite

from gaussian_radial_reference import GaussianRadialSource
from higher_local_geometry import local_geometry, compare_legacy
from reduced_gaunt_exact import add, gaussian_moment, tensor_square
from reduced_gaunt_vertices import tensor_from_labels
from spherical_green_reference import RadialGrid, OriginSeries, ExteriorSource
from spherical_york_reference import ScalarSourceSpec, solve_scalar_source, solve_york_source
from solid_harmonic_reference import AngularGrid


def _finite(value, name):
    if isinstance(value, bool):
        raise ValueError(name+" must be a finite number")
    value = float(value)
    if not isfinite(value):
        raise ValueError(name+" must be finite")
    return value


def _flat(value):
    if isinstance(value, (tuple, list)):
        return tuple(x for part in value for x in _flat(part))
    return (_finite(value, "comparison value"),)


def _comparison(actual, reference, tolerance):
    left, right = _flat(actual), _flat(reference)
    if len(left) != len(right) or not left:
        raise ValueError("nonempty comparison shapes must agree")
    changes = [abs(a-b)/max(1.0, abs(b)) for a, b in zip(left, right)]
    return {"actual": actual, "reference": reference, "scaled_errors": changes,
            "maximum": max(changes), "tolerance": tolerance,
            "passed": max(changes) <= tolerance}


def _settings(config):
    if not isinstance(config, Mapping):
        raise TypeError("config must be a mapping")
    result = {key: config.get("first_response_"+key, default) for key, default in {
        "intervals": [256, 512, 1024], "radius": 4, "origin_radius": .5,
        "origin_degree": 24, "n_mu": 8, "n_phi": 16,
        "angular_n_mu": 10, "angular_n_phi": 20,
        "scalar_ell_max": 4, "vector_ell_max": 5, "h_lmax": 6, "guard_lmax": 7,
        "value_tolerance": 3e-4, "component_tolerance": 3e-3,
        "support_tolerance": 1e-10, "provider_precision": 64,
        "provider_max_products": 5_000_000, "provider_max_terms": 100_000,
        "provider_max_moment_steps": 2_000_000, "provider_max_series_terms": 2048,
        "provider_max_point_evaluations": 300_000,
    }.items()}
    intervals = tuple(result["intervals"])
    if len(intervals) != 3 or any(type(n) is not int or n < 64 or n > 4096 for n in intervals):
        raise ValueError("first-response intervals require exactly three bounded integer levels")
    if any(a >= b for a, b in zip(intervals, intervals[1:])):
        raise ValueError("first-response intervals must increase")
    result["intervals"] = list(intervals)
    for key in ("radius", "origin_radius", "value_tolerance", "component_tolerance", "support_tolerance"):
        result[key] = _finite(result[key], key)
        if result[key] <= 0:
            raise ValueError(key+" must be positive")
    if result["radius"] != 4 or result["origin_radius"] != .5:
        raise ValueError("This declared first-response pilot uses R4 and fixed patch .5")
    for key in ("origin_degree", "n_mu", "n_phi", "angular_n_mu", "angular_n_phi",
                "scalar_ell_max", "vector_ell_max", "h_lmax", "guard_lmax", "provider_precision",
                "provider_max_products", "provider_max_terms", "provider_max_moment_steps",
                "provider_max_series_terms", "provider_max_point_evaluations"):
        if type(result[key]) is not int or result[key] <= 0:
            raise ValueError(key+" must be a positive integer")
    if (result["scalar_ell_max"], result["vector_ell_max"], result["h_lmax"], result["guard_lmax"]) != (4, 5, 6, 7):
        raise ValueError("The pilot retains scalar l<=4, vector l<=5, h l<=6 and guard l=7")
    for nmu, nphi in ((result["n_mu"], result["n_phi"]), (result["angular_n_mu"], result["angular_n_phi"])):
        if nmu < 8 or nphi < 15 or nmu*nphi > 1024:
            raise ValueError("sphere grids must cover the complete guard band within the point cap")
    if result["angular_n_mu"] <= result["n_mu"] or result["angular_n_phi"] <= result["n_phi"]:
        raise ValueError("The independent angular slice must refine both sphere orders")
    if not 16 <= result["origin_degree"] <= 32 or not 64 <= result["provider_precision"] <= 128:
        raise ValueError("bounded source Taylor degree and Decimal precision required")
    return result


def _points(radius):
    return (((0., 0., 0.), None), ((.2, -.15, .1), None),
            ((.5, 0., 0.), "left"), ((.5, 0., 0.), "right"),
            ((.7, -.3, .2), None), ((1.2, .6, -.4), None), ((2., -1., .5), None),
            ((radius, 0., 0.), "left"), ((radius, 0., 0.), "right"),
            ((4.8, -1.2, .6), None), ((7., -3., 2.), None))


def _source_record(poly):
    return [{"powers": list(powers), "coefficient": str(c)} for powers, c in sorted(poly.items())]


def _heat_targets(polynomials, points):
    # This import is intentionally confined to the comparison side.
    from gaussian_green_fields import scalar_field, york_field, evaluate_fields
    scalar = {name: scalar_field(polynomials[name], F(1)) for name in ("u", "v")}
    vector = york_field(tuple(polynomials["b1_"+str(i)] for i in range(3)), F(1))
    targets = {}
    field_groups = [(name, (field,)) for name, field in scalar.items()]+[("W", vector)]
    for name, fields in field_groups:
        components = len(fields)
        compiled = tuple(fields)+tuple(field.derivative(a) for a in range(3) for field in fields)
        compiled += tuple(field.derivative(a).derivative(b) for a in range(3) for b in range(3) for field in fields)
        rows = {}
        for point, _ in points:
            if point in rows:
                continue
            values = evaluate_fields(compiled, point)
            if components == 1:
                row = {"value": values[0], "gradient": list(values[1:4]),
                       "hessian": [list(values[4+3*a:7+3*a]) for a in range(3)]}
            else:
                row = {"value": list(values[:3]), "gradient": [list(values[3+3*a:6+3*a]) for a in range(3)],
                       "hessian": [[list(values[12+9*a+3*b:15+9*a+3*b]) for b in range(3)] for a in range(3)]}
            for value in row.values():
                _flat(value)
            rows[point] = row
        targets[name] = rows
    return targets


def run_first_response_checks(config, budget, progress=None):
    """Three radial levels plus one isolated angular slice, with all-m sources."""
    current = "settings"
    result = {"passed": False, "levels": [], "rows": [], "refinements": [],
              "direction": "monopole_mixed", "full_H6_validation": False,
              "independent_mixed_functional_validation": False, "accuracy_certified": False}

    def emit(stage, **payload):
        if progress is not None:
            progress({"stage": stage, **payload})

    def record(row):
        result["rows"].append(row)
        emit("first_response_row", row=row)

    try:
        settings = _settings(config)
        result["settings"], result["work_before"] = settings, budget.snapshot()
        shell_points = settings["n_mu"]*settings["n_phi"]*sum(n+2 for n in settings["intervals"])
        shell_points += settings["angular_n_mu"]*settings["angular_n_phi"]*(settings["intervals"][-1]+2)
        target_points = 4*len(_points(settings["radius"]))
        planned = {"projection_products": 2*(settings["guard_lmax"]+1)**2*shell_points,
                   "point_evaluations": shell_points+3*target_points,
                   "cartesian_calls": 3*shell_points+6*target_points}
        result["planned_adapter_work"] = planned
        emit("first_response_work_preflight", planned=dict(planned), remaining=budget.snapshot())
        for key, count in planned.items():
            if count > budget.remaining(key):
                raise ValueError("Insufficient remaining first-response work cap: "+key)
        current = "actual_geometry"
        zero = [[F(0)]*3 for _ in range(3)]
        axis = [[-F(1, 2), F(0), F(0)], [F(0), -F(1, 2), F(0)], [F(0), F(0), F(1)]]
        xz = [[F(0), F(0), F(1)], [F(0), F(0), F(0)], [F(1), F(0), F(0)]]
        q, s = tensor_from_labels(axis, zero), tensor_from_labels(zero, xz)
        geometry = local_geometry(q, s, local_R_order=2, density_order=4, york_order=1,
                                  progress=lambda event: emit("first_response_geometry_progress", geometry_event=event))
        legacy = compare_legacy(q, s, geometry)
        s2, r2 = tensor_square(s), geometry["R"][2]
        polynomials = {"u": add(s2, r2, -1), "v": add(r2, s2, 7),
                       **{"b1_"+str(i): poly for i, poly in enumerate(geometry["b"][1])}}
        normalized = lambda poly: sum((c*gaussian_moment(powers) for powers, c in poly.items()), F(0))
        masses = {"R2": normalized(r2), "s2": normalized(s2),
                  **{"b1_"+str(i): normalized(poly) for i, poly in enumerate(geometry["b"][1])}}
        mass_targets = {"R2": -F(6237, 128), "s2": F(2079, 8),
                        "b1_0": F(0), "b1_1": F(2079, 32), "b1_2": F(0)}
        rates_passed = geometry["gaussian_rates"]["R"][2] == 2 and geometry["gaussian_rates"]["b"][1] == 2
        anchor = {"name": "actual_source_anchors", "legacy": legacy, "physical_beta": "1",
                  "geometry_rate_is_twice_beta": rates_passed,
                  "moments_without_pi_3_over_2": {name: str(value) for name, value in masses.items()},
                  "targets_without_pi_3_over_2": {name: str(value) for name, value in mass_targets.items()},
                  "passed": legacy["passed"] and rates_passed and masses == mass_targets,
                  "acceptance_required": True, "role": "structural"}
        record(anchor)
        result["sources"] = {name: _source_record(poly) for name, poly in polynomials.items()}
        emit("first_response_source_polynomials", sources=result["sources"], anchors=anchor)
        if not anchor["passed"]:
            raise AssertionError("Actual physical source anchors failed")

        current = "source_fischer_and_targets"
        source_options = {"precision": settings["provider_precision"],
            **{key: settings["provider_"+key] for key in ("max_products", "max_terms", "max_moment_steps", "max_series_terms", "max_point_evaluations")}}
        sources = {name: GaussianRadialSource(poly, F(1), ell_max=settings["scalar_ell_max"] if name in ("u", "v") else settings["vector_ell_max"],
                                            **source_options) for name, poly in polynomials.items()}
        origin_records = {name: {mode: source.origin_data(*mode, degree=settings["origin_degree"],
            patch_radius=F(1, 2), derivative_order=4) for mode in source.mode_keys}
                          for name, source in sources.items()}
        result["source_provenance"] = {name: source.provenance() for name, source in sources.items()}
        emit("first_response_source_preparation", sources=result["source_provenance"])
        mode_count = sum(len(source.mode_keys) for source in sources.values())+(settings["h_lmax"]+1)**2
        fd_count = mode_count*(sum(n+1 for n in settings["intervals"])+settings["intervals"][-1]+1)
        result["planned_adapter_work"]["fd_unknowns"] = fd_count
        emit("first_response_fd_preflight", mode_count=mode_count, planned_fd_unknowns=fd_count,
             remaining_fd_unknowns=budget.remaining("fd_unknowns"), work=budget.snapshot())
        if fd_count > budget.remaining("fd_unknowns"):
            raise ValueError("Insufficient remaining first-response FD unknown cap")
        points = _points(settings["radius"])
        targets = _heat_targets(polynomials, points)
        result["comparison_targets"] = [{"field": name, "point": point, "jets": jets}
                                         for name, table in targets.items() for point, jets in table.items()]
        emit("first_response_targets", comparison_targets=result["comparison_targets"])

        def specification(name, grid):
            source = sources[name]
            tail = source.exterior(grid.radius)
            origins = {mode: OriginSeries(row["patch_radius"], row["coefficients"], row["error_bounds_float_s"])
                       for mode, row in origin_records[name].items()}
            return ScalarSourceSpec(source.sample(grid.radii), origins,
                ExteriorSource(grid.radius, {}, tail.weighted_outer_bounds()),
                residual_policy="exact_correction", residual_provider=tail,
                provenance={"direction": "monopole_mixed", "field": name, "beta": "1",
                            "source": source.provenance(), "origin_degree": settings["origin_degree"],
                            "origin_bounds_include_coefficient_conversion": True})

        def level(intervals, nmu, nphi, kind, required):
            nonlocal current
            current = kind+"_"+str(intervals)
            grid = RadialGrid(settings["radius"], intervals, (settings["origin_radius"],))
            angular = AngularGrid(settings["guard_lmax"], nmu, nphi)
            def source_progress(event):
                emit("first_response_adapter", intervals=intervals, slice=kind, adapter_event=event)
            scalar = {name: solve_scalar_source(grid, specification(name, grid), inverse="C", budget=budget,
                                                progress=source_progress) for name in ("u", "v")}
            york = solve_york_source(grid, tuple(specification("b1_"+str(i), grid) for i in range(3)), angular,
                h_lmax=settings["h_lmax"], guard_lmax=settings["guard_lmax"], support_tolerance=settings["support_tolerance"],
                budget=budget, progress=source_progress)
            fields = {**scalar, "W": york}
            output = {"intervals": intervals, "radius": grid.radius, "sphere": [nmu, nphi], "slice": kind,
                      "points": {}, "diagnostics": {name: field.diagnostics() for name, field in fields.items()}}
            emit("first_response_fields_ready", intervals=intervals, slice=kind, diagnostics=output["diagnostics"])
            checks = []
            for point, side in points:
                label = str((point, side))
                output["points"][label] = {}
                for name, field in fields.items():
                    actual = field.jet(point, order=2, side=side)
                    output["points"][label][name] = {key: actual[key] for key in ("value", "gradient", "hessian")}
                    for key in ("value", "gradient", "hessian"):
                        tolerance = settings["value_tolerance"] if key == "value" else settings["component_tolerance"]
                        comparison = _comparison(actual[key], targets[name][point][key], tolerance)
                        row = {"name": name+"_"+key, "point": point, "side": side,
                               "intervals": intervals, "slice": kind, **comparison,
                               "acceptance_required": required, "role": "physical_target"}
                        record(row)
                        checks.append(row)
            output["raw_passed"] = all(row["passed"] for row in checks)
            output["passed"] = output["raw_passed"] if required else True
            output["work"] = budget.snapshot()
            output["provider_statistics"] = {name: source.statistics() for name, source in sources.items()}
            result["levels"].append(output)
            emit("first_response_level_result", level=output)
            return output

        def compare_levels(old, new, kind, required):
            rows = []
            for name in ("u", "v", "W"):
                for key in ("value", "gradient", "hessian"):
                    tolerance = settings["value_tolerance"] if key == "value" else settings["component_tolerance"]
                    comparisons = {label: _comparison(new["points"][label][name][key], old["points"][label][name][key], tolerance)
                                   for label in old["points"]}
                    maximum = max(item["maximum"] for item in comparisons.values())
                    row = {"name": kind+"_"+name+"_"+key, "slice": kind,
                           "intervals": [old["intervals"], new["intervals"]], "maximum": maximum,
                           "tolerance": tolerance, "passed": maximum <= tolerance,
                           "point_comparisons": comparisons, "acceptance_required": required, "role": "refinement"}
                    record(row)
                    rows.append(row)
            result["refinements"].append({"kind": kind, "rows": rows,
                                          "passed": all(row["passed"] for row in rows) if required else True})

        previous = None
        for index, intervals in enumerate(settings["intervals"]):
            latest = index == len(settings["intervals"])-1
            new = level(intervals, settings["n_mu"], settings["n_phi"], "radial", latest)
            if previous is not None:
                compare_levels(previous, new, "radial", latest)
            previous = new
        changed = level(settings["intervals"][-1], settings["angular_n_mu"], settings["angular_n_phi"], "angular", True)
        compare_levels(previous, changed, "angular", True)
        required = [row for row in result["rows"] if row["acceptance_required"]]
        result.update(passed=bool(required) and all(row["passed"] for row in required),
            required_failures=sum(not row["passed"] for row in required),
            diagnostic_failures=sum(not row["passed"] for row in result["rows"] if not row["acceptance_required"]),
            work_after=budget.snapshot(), provider_statistics={name: source.statistics() for name, source in sources.items()},
            limitations=["Physical source geometry and raw basis are shared with independently checked local modules.",
                         "Only first responses are compared; W2,u3,u4, canonical integrals and mixed H6 remain pending.",
                         "Finite samples and refinement are numerical diagnostics, not field error certificates.",
                         "This pilot isolates angle and radius-grid spacing at R4; a separate exterior-radius slice remains pending.",
                         "HEAT comparison fields never supply boundary/origin/FD data to the independent route."])
        return result
    except Exception as exc:
        failure = {"name": current, "passed": False, "acceptance_required": True,
                   "exception_type": type(exc).__name__, "exception": str(exc)}
        emit("first_response_exception", control_result=failure)
        raise
