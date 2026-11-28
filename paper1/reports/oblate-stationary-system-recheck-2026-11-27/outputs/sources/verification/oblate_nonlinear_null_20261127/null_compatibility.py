#!/usr/bin/env python3
"""Prepared, unexecuted maps for Paper I's second-order null compatibility.

All fields here are axisymmetric. Harmonic conventions match all_modes.py at
commit 157ec0f04866ea36509b325ec65f0a68890f9b17. The scalar/tensor coefficient
projection is the fixed round coefficient projection, not a deformed DeWitt
projection. A candidate quadratic source is NOT an obstruction unless its
input satisfies the complete retained linearized equations and constraints.
"""
from dataclasses import dataclass

import sympy as s


@dataclass
class Fields:
    a: object
    b: object
    lam: object
    gamma: object
    V: object
    W: object


def angular_basis(ell, x):
    """Unnormalized real m=0 bases, with the pinned polar/axial signs."""
    d = 1-x*x
    Y = s.legendre(ell, x)
    polar = s.Matrix([d*s.diff(Y, x), 0])
    axial = s.Matrix([0, -s.diff(Y, x)])
    E = s.diag(s.diff(Y, x, 2)/2, -d*d*s.diff(Y, x, 2)/2)
    B = s.Matrix([[0, -d*s.diff(Y, x, 2)/2],
                  [-d*s.diff(Y, x, 2)/2, 0]])
    return Y, polar, axial, E, B


def coefficient_labels(L=2):
    """Named columns; this is an explicit map, not a change of canonical order."""
    q = []
    shifts = []
    for ell in range(L+1):
        q.extend([("lambda", ell), ("X", ell)])
        shifts.append(("b", ell))
        if ell >= 1:
            q.extend([("V_polar", ell), ("V_axial", ell)])
            shifts.extend([("W_polar", ell), ("W_axial", ell)])
        if ell >= 2:
            q.extend([("C_polar", ell), ("C_axial", ell)])
    z = [("a", ell) for ell in range(L+1)] + shifts
    return q, z


def coefficient_fields(coefficients, x):
    d = 1-x*x
    values = {"a": s.S.Zero, "b": s.S.Zero,
              "lambda": s.S.Zero, "X": s.S.Zero}
    V, W, C = s.zeros(2, 1), s.zeros(2, 1), s.zeros(2)
    for (kind, ell), amplitude in coefficients.items():
        Y, polar, axial, E, B = angular_basis(ell, x)
        if kind in values:
            values[kind] += amplitude*Y
        elif kind == "V_polar":
            V += amplitude*polar
        elif kind == "V_axial":
            V += amplitude*axial
        elif kind == "W_polar":
            W += amplitude*polar
        elif kind == "W_axial":
            W += amplitude*axial
        elif kind == "C_polar":
            C += amplitude*E
        elif kind == "C_axial":
            C += amplitude*B
        else:
            raise ValueError(f"Unknown coefficient {kind}")
    gamma = values["X"]*s.diag(1/d, d)+C
    return Fields(values["a"], values["b"], values["lambda"], gamma, V, W)


def background(coords, epsilon):
    _, r, x, _ = coords
    d, rho, Sigma = 1-x*x, r*r+epsilon, r*r+epsilon*x*x
    return s.diag(-1, Sigma/rho, Sigma/d, rho*d)


def reconstruct(fields):
    """Full nonlinear coefficient-to-metric map, coordinate order t,r,x,phi."""
    A = s.Matrix.hstack(fields.W, fields.V)
    h = s.Matrix([[fields.a, fields.b], [fields.b, fields.lam]])
    g = s.zeros(4)
    g[:2, :2] = h+A.T*fields.gamma*A
    g[:2, 2:] = A.T*fields.gamma
    g[2:, :2] = fields.gamma*A
    g[2:, 2:] = fields.gamma
    return g


def linear_metric(fields, coords, epsilon):
    """Derivative of reconstruct at the static oblate coefficient point."""
    g0 = background(coords, epsilon)
    gamma0 = g0[2:, 2:]
    h = s.zeros(4)
    h[0, 0], h[0, 1], h[1, 0], h[1, 1] = (
        fields.a, fields.b, fields.b, fields.lam)
    h[0, 2:] = (gamma0*fields.W).T
    h[2:, 0] = gamma0*fields.W
    h[1, 2:] = (gamma0*fields.V).T
    h[2:, 1] = gamma0*fields.V
    h[2:, 2:] = fields.gamma
    return h


def linear_angular_generator(fields, x, w, angular_profile=None):
    """DF_xi[c1] for xi=w(r)*P(x)*partial_phi and m=0 c1.

    The affine partial_r xi term is in F_xi[c0], not in its derivative.
    This formula keeps the W (temporal) variation.
    """
    P = x if angular_profile is None else angular_profile
    factor = w*s.diff(P, x)
    gamma = s.Matrix([[2*factor*fields.gamma[0, 1],
                       factor*fields.gamma[1, 1]],
                      [factor*fields.gamma[1, 1], 0]])
    V = s.Matrix([0, -factor*fields.V[0]])
    W = s.Matrix([0, -factor*fields.W[0]])
    return Fields(s.S.Zero, s.S.Zero, s.S.Zero, gamma, V, W)


def average(expr, x):
    return s.integrate(s.cancel(expr), (x, -1, 1))/2


def project_coefficients(fields, x, max_degree=4):
    """Exact round positive tensor/vector projection in coefficient space."""
    d = 1-x*x
    X = (d*fields.gamma[0, 0]+fields.gamma[1, 1]/d)/2
    C = fields.gamma-X*s.diag(1/d, d)
    result = {}
    for ell in range(max_degree+1):
        Y, P, A, E, B = angular_basis(ell, x)
        scalar_norm = s.Rational(1, 2*ell+1)
        for kind, value in [("a", fields.a), ("b", fields.b),
                            ("lambda", fields.lam), ("X", X)]:
            result[kind, ell] = s.cancel(average(value*Y, x)/scalar_norm)
        if ell >= 1:
            vector_norm = ell*(ell+1)*scalar_norm
            for kind, value in [("V", fields.V), ("W", fields.W)]:
                for parity, basis in [("polar", P), ("axial", A)]:
                    pairing = value[0]*basis[0]/d+d*value[1]*basis[1]
                    result[kind+"_"+parity, ell] = s.cancel(
                        average(pairing, x)/vector_norm)
        if ell >= 2:
            tensor_norm = (ell-1)*ell*(ell+1)*(ell+2)*scalar_norm/2
            for parity, basis in [("polar", E), ("axial", B)]:
                pairing = (d*d*C[0, 0]*basis[0, 0]
                           + 2*C[0, 1]*basis[0, 1]
                           + C[1, 1]*basis[1, 1]/(d*d))
                result["C_"+parity, ell] = s.cancel(
                    average(pairing, x)/tensor_norm)
    return result


def expected_discarded_generator(coefficients, w):
    """Hand-derived eight-row matrix, pending the prepared independent check."""
    value = lambda kind, ell: coefficients.get((kind, ell), s.S.Zero)
    odd = value("C_axial", 2)
    return {
        ("X", 4): -s.Rational(12, 35)*w*odd,
        ("C_polar", 4): s.Rational(2, 35)*w*odd,
        ("C_axial", 3): -s.Rational(2, 15)*w*value("X", 1),
        ("C_axial", 4): -s.Rational(2, 35)*w*(
            value("X", 2)+value("C_polar", 2)),
        ("V_axial", 3): -s.Rational(2, 15)*w*value("V_polar", 1),
        ("V_axial", 4): -s.Rational(6, 35)*w*value("V_polar", 2),
        ("W_axial", 3): -s.Rational(2, 15)*w*value("W_polar", 1),
        ("W_axial", 4): -s.Rational(6, 35)*w*value("W_polar", 2),
    }


def christoffel(metric, coords):
    inv = metric.inv()
    n = len(coords)
    return [[[s.cancel(sum(inv[a, k]*(
        s.diff(metric[k, c], coords[b])
        + s.diff(metric[k, b], coords[c])
        - s.diff(metric[b, c], coords[k])) for k in range(n))/2)
        for c in range(n)] for b in range(n)] for a in range(n)]


def linearized_einstein(metric, perturbation, coords):
    """FULL covariant delta G on the Ricci-flat background; no projection.

    This implementation is prepared but has not been scientifically executed.
    It retains coordinate derivatives until all variations are assembled.
    """
    n = len(coords)
    inv = metric.inv()
    dinv = -inv*perturbation*inv
    Gamma = christoffel(metric, coords)
    dGamma = [[[s.cancel(sum(
        dinv[a, k]*(s.diff(metric[k, c], coords[b])
                    + s.diff(metric[k, b], coords[c])
                    - s.diff(metric[b, c], coords[k]))
        + inv[a, k]*(s.diff(perturbation[k, c], coords[b])
                     + s.diff(perturbation[k, b], coords[c])
                     - s.diff(perturbation[b, c], coords[k]))
        for k in range(n))/2)
        for c in range(n)] for b in range(n)] for a in range(n)]
    dRicci = s.zeros(n)
    for i in range(n):
        for j in range(n):
            derivative = sum(s.diff(dGamma[k][i][j], coords[k])
                             - s.diff(dGamma[k][i][k], coords[j])
                             for k in range(n))
            products = sum(
                dGamma[k][i][j]*Gamma[l][k][l]
                + Gamma[k][i][j]*dGamma[l][k][l]
                - dGamma[k][i][l]*Gamma[l][k][j]
                - Gamma[k][i][l]*dGamma[l][k][j]
                for k in range(n) for l in range(n))
            dRicci[i, j] = s.cancel(derivative+products)
    dR = sum(inv[i, j]*dRicci[i, j] for i in range(n) for j in range(n))
    return (dRicci-metric*dR/2).applyfunc(s.cancel)


def euler_pair_density(delta_G, test_fields, coords, epsilon):
    """Full coefficient Euler density, to integrate with dr dx/2 for m=0.

    delta_G has covariant indices. The normalized action is (16piG)/(4pi) I.
    The factor Sigma is the background four-volume density; do not add r^2.
    """
    _, r, x, _ = coords
    g0 = background(coords, epsilon)
    inv = g0.inv()
    raised = inv*delta_G*inv
    test = linear_metric(test_fields, coords, epsilon)
    return s.cancel(-(r*r+epsilon*x*x)*sum(
        raised[i, j]*test[i, j] for i in range(4) for j in range(4)))


def candidate_source_density(coefficients, delta_G, coords, epsilon, w):
    """Return the integrand of B_xi[c1], WITHOUT asserting on-shell status.

    The caller must separately certify all retained linearized equations and
    constraints, the radial domain, and consistent derivative data. Off-shell
    nonzero output is not a nonlinear compatibility obstruction.
    """
    defect = expected_discarded_generator(coefficients, w)
    return euler_pair_density(delta_G, coefficient_fields(defect, coords[2]),
                              coords, epsilon)
