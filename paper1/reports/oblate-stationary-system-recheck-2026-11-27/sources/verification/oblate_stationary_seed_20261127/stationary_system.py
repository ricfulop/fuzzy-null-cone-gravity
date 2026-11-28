"""Prepared stationary radial equations; no solution or execution is asserted.

This uses the original Paper I full linearized Einstein map. No radial
projection, endpoint condition, seed, or differential-algebraic reduction is
chosen here. Angular moments are exact for r>0, epsilon>0.
"""
from functools import lru_cache
from contextlib import nullcontext
from pathlib import Path
import sys

import sympy as s

MAP_DIRECTORY = Path(__file__).resolve().parents[1] / "oblate_nonlinear_null_20261127"
sys.path.insert(0, str(MAP_DIRECTORY))
import null_compatibility as n


COLUMNS = (
    ("lambda", 0), ("lambda", 2), ("X", 0), ("X", 2),
    ("V_polar", 2), ("C_polar", 2), ("a", 0), ("a", 2),
)
RESPONSE_TEST = {("X", 4): s.Integer(6), ("C_polar", 4): -s.S.One}


class AngularMoments:
    """Average rational densities with denominator proportional to Sigma**k.

    The recurrence is for integral_0^1 x**(2n)/(r**2+epsilon*x**2)**k dx.
    Denominators outside this class fail explicitly rather than being dropped.
    """
    def __init__(self, r, epsilon, x):
        self.r, self.epsilon, self.x = r, epsilon, x
        self.requested = set()

    @lru_cache(maxsize=None)
    def moment(self, degree, power):
        r, epsilon = self.r, self.epsilon
        if degree < 0 or power < 0:
            raise ValueError("Moment indices must be nonnegative")
        if power == 0:
            return s.Rational(1, 2*degree+1)
        if degree:
            return (self.moment(degree-1, power-1)
                    - r*r*self.moment(degree-1, power))/epsilon
        if power == 1:
            return s.atan(s.sqrt(epsilon)/r)/(r*s.sqrt(epsilon))
        return (1 / (2*r*r*(power-1)*(r*r+epsilon)**(power-1))
                + s.Rational(2*power-3, 2*(power-1))
                * self.moment(0, power-1)/(r*r))

    def average(self, expression):
        x, r, epsilon = self.x, self.r, self.epsilon
        numerator, denominator = s.fraction(s.cancel(expression))
        denominator_degree = s.Poly(denominator, x).degree()
        if denominator_degree % 2:
            raise ValueError("Unsupported angular denominator")
        power = denominator_degree // 2
        scale = s.cancel(denominator/(r*r+epsilon*x*x)**power)
        if scale.has(x):
            raise ValueError("Angular denominator is not proportional to Sigma**k")
        answer = s.S.Zero
        for (degree,), coefficient in s.Poly(numerator, x).terms():
            if degree % 2 == 0:
                self.requested.add((degree//2, power))
                answer += coefficient*self.moment(degree//2, power)/scale
        return s.cancel(answer)


def radial_coefficients(expression, profile, r):
    """Extract coefficients of F,F',F'' without freezing background radius."""
    slots = s.symbols("radial_value radial_first radial_second")
    replacements = {s.diff(profile, r, order): slots[order] for order in range(3)}
    expanded = s.cancel(expression).xreplace(replacements)
    if expanded.has(profile):
        raise ValueError("Unexpected radial derivative order")
    polynomial = s.Poly(expanded, *slots)
    if polynomial.total_degree() > 1 or polynomial.coeff_monomial(1) != 0:
        raise ValueError("Expected a homogeneous linear radial differential expression")
    return [s.cancel(polynomial.coeff_monomial(slot)) for slot in slots]


def static_split_residual(delta_G, fields, coords, epsilon):
    """4D static split cross-check; reuses the same curvature implementation."""
    _, r, x, _ = coords
    g0 = n.background(coords, epsilon)
    sigma = g0[1:, 1:]
    h = n.linear_metric(fields, coords, epsilon)
    G3 = n.linearized_einstein(sigma, h[1:, 1:], coords[1:])
    inverse = sigma.inv()
    Gamma = n.christoffel(sigma, coords[1:])
    Hess = s.Matrix(3, 3, lambda i, j:
        s.diff(fields.a, coords[i+1], coords[j+1])
        - sum(Gamma[k][i][j]*s.diff(fields.a, coords[k+1]) for k in range(3)))
    laplace = sum(inverse[i, j]*Hess[i, j] for i in range(3) for j in range(3))
    expected = s.zeros(4)
    expected[0, 0] = -sum(inverse[i, j]*G3[i, j] for i in range(3) for j in range(3))
    expected[1:, 1:] = G3 + (Hess-sigma*laplace)/2
    return (delta_G-expected).applyfunc(s.cancel)


def assemble(coords, epsilon, recorder=None):
    """Return all 22 retained equations and one response row, before DAE solving."""
    _, r, x, _ = coords
    q, z = n.coefficient_labels()
    rows = q+z
    selected = [rows.index(label) for label in COLUMNS]
    automatic = [index for index in range(len(rows)) if index not in selected]
    all_tests = [n.coefficient_fields({label: s.S.One}, x) for label in rows]
    all_tests.append(n.coefficient_fields(RESPONSE_TEST, x))
    moments = AngularMoments(r, epsilon, x)
    # Last row is D=6 E_X4-E_Cpolar4, not an additional retained equation.
    operators = [s.zeros(len(rows)+1, len(COLUMNS)) for _ in range(3)]
    densities = [s.zeros(len(rows)+1, len(COLUMNS)) for _ in range(3)]
    profile = s.Function("F")(r)
    stage = recorder.stage if recorder is not None else lambda *args, **kwargs: nullcontext()
    for column, label in enumerate(COLUMNS):
        with stage("coefficient_reconstruction", column=column, label=label):
            fields = n.coefficient_fields({label: profile}, x)
        with stage("four_dimensional_curvature", column=column, label=label):
            G = n.linearized_einstein(
                n.background(coords, epsilon), n.linear_metric(fields, coords, epsilon), coords)
        with stage("static_split_control", column=column, label=label):
            if any(value != 0 for value in static_split_residual(G, fields, coords, epsilon)):
                raise AssertionError(("4D/3D static split", label))
        for row, test in enumerate(all_tests):
            with stage("euler_density", column=column, row=row):
                density = n.euler_pair_density(G, test, coords, epsilon)
            with stage("radial_coefficient_extraction", column=column, row=row):
                parts = radial_coefficients(density, profile, r)
            for order, part in enumerate(parts):
                densities[order][row, column] = part
                with stage("angular_average", column=column, row=row, order=order):
                    operators[order][row, column] = moments.average(part)
        if recorder is not None:
            with stage("column_checkpoint", column=column, label=label):
                if any(matrix[row, column].has(s.nan, s.zoo, s.oo, -s.oo)
                       for matrix in operators+densities for row in range(len(rows)+1)):
                    raise ValueError("Nonfinite exact expression in column checkpoint")
                recorder.save_column(column, label, {
                    "retained_row_labels": rows,
                    "selected_row_indices": selected,
                    "automatic_row_indices": automatic,
                    "last_row": "D=6 E_X4-E_Cpolar4; discarded response, not a retained equation",
                    "column_labels": COLUMNS,
                    "derivative_orders": [0, 1, 2],
                    "operators_by_derivative_order": [
                        [str(matrix[row, column]) for row in range(len(rows)+1)]
                        for matrix in operators],
                    "densities_by_derivative_order": [
                        [str(matrix[row, column]) for row in range(len(rows)+1)]
                        for matrix in densities],
                    "requested_moment_indices_through_column": sorted(moments.requested),
                    "static_split_control_passed_for_column": True,
                    "all_global_controls_passed": False,
                    "operator_convention": "sum_j M_j(r,epsilon) u^(j)(r); no derivative on M_j",
                    "angular_measure": "dx/2; Sigma occurs once in the Euler density",
                })
        print(f"ASSEMBLED stationary column {label}", flush=True)
    with stage("all_omitted_retained_rows_control"):
        for order in range(3):
            if any(operators[order][row, col] != 0
                   for row in automatic for col in range(len(COLUMNS))):
                raise AssertionError("An omitted retained row is nonzero")
    return {"rows": rows, "columns": COLUMNS, "selected_row_indices": selected,
            "automatic_row_indices": automatic, "operators": operators,
            "densities": densities, "moments": moments}


def apply_operator(operators, profiles, r):
    """Apply differential expressions to functions; no admissibility inference."""
    profiles = s.Matrix(profiles)
    return sum((matrix*profiles.diff(r, order)
                for order, matrix in enumerate(operators)),
               s.zeros(operators[0].rows, 1)).applyfunc(s.cancel)


def first_order_family_residual(operators, embedding, evolution, r):
    """Future certificate for u=M z, z'=K z; no M,K are supplied here.

    All retained output rows must vanish identically on an interval. The last
    row then gives the response on that actual ODE family, not on free jets.
    This function does not construct or certify the smoothness of M or K.
    """
    M, K = s.Matrix(embedding), s.Matrix(evolution)
    if M.rows != operators[0].cols or K.rows != K.cols or M.cols != K.rows:
        raise ValueError("Incompatible first-order family dimensions")
    first = M.diff(r)+M*K
    second = M.diff(r, 2)+2*M.diff(r)*K+M*K.diff(r)+M*K*K
    return (operators[0]*M+operators[1]*first+operators[2]*second).applyfunc(s.cancel)


def radial_gauge_profiles(r, epsilon, f):
    rho = r*r+epsilon
    lambda0 = (r*r+epsilon/3)/rho
    lambda2 = 2*epsilon/(3*rho)
    return [s.diff(lambda0, r)*f+2*lambda0*s.diff(f, r),
            s.diff(lambda2, r)*f+2*lambda2*s.diff(f, r),
            2*r*f, 0, 0, 0, 0, 0]


def oblate_parameter_profiles(r, epsilon):
    rho = r*r+epsilon
    return [-2*r*r/(3*rho*rho), 2*r*r/(3*rho*rho),
            s.Rational(2, 3), s.Rational(1, 3), 0, -s.Rational(1, 3), 0, 0]


def odd_coordinate_leg(r, epsilon, f):
    return {("C_axial", 2): -s.Rational(2, 3)*(r*r+epsilon)*f,
            ("V_axial", 2): -s.diff(f, r)/3}


def response_prefactor(r, epsilon, f, w):
    return s.Rational(4, 105)*(r*r+epsilon)*f*w
