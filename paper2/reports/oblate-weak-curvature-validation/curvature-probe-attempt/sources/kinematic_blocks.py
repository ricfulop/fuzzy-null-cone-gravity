"""Prepared, unexecuted kinematic building blocks for the original Paper I action.

The formulas retain general r and the actual coefficient lapse/shift map.
There is no default integration, matrix-rank inference, or scientific main.
Use only after the separately proposed verification has been approved.

Harmonic conventions follow Paper I all_modes.py at commit
157ec0f04866ea36509b325ec65f0a68890f9b17 (source SHA256
55ba3d4ba6c44104668362258aa96465bd485fef8f407f65283e12aa0c643ebe).
General-r maps and the canonical assembly are new preparation, not that
archived verifier's tested output.
"""

import sympy as sp


def harmonic_amplitudes(ell, m, x):
    """Associated-Legendre amplitudes; exp(i*m*phi) is implicit."""
    d = 1 - x**2
    y = sp.assoc_legendre(ell, m, x)
    polar = sp.Matrix([d * sp.diff(y, x), sp.I * m * y / d])
    axial = sp.Matrix([sp.I * m * y, -sp.diff(y, x)])
    electric = sp.zeros(2)
    electric[0, 0] = (
        sp.diff(y, x, 2) - x * sp.diff(y, x) / d
        + ell * (ell + 1) * y / (2 * d)
    )
    electric[0, 1] = electric[1, 0] = sp.I * m * (sp.diff(y, x) + x * y / d)
    electric[1, 1] = -d**2 * electric[0, 0]
    magnetic = sp.Matrix([
        [electric[0, 1] / d, -d * electric[0, 0]],
        [-d * electric[0, 0], -d * electric[0, 1]],
    ])
    return y, polar, axial, electric, magnetic


def sector(L, m):
    """Return unintegrated tensors and labels for one complex angular sector.

    Domain: r>0, -1<x<1, rho>0, Sigma>0; epsilon is real. The radial
    profile is not frozen: shift columns are U0*s + T*J1*s_prime.
    """
    if not isinstance(L, int) or not isinstance(m, int) or not (L >= 2 and 0 <= m <= L):
        raise ValueError("Require integer L>=2 and 0<=m<=L")
    r = sp.Symbol("r", positive=True)
    x = sp.Symbol("x", real=True)
    epsilon = sp.Symbol("epsilon", real=True)
    d = 1 - x**2
    rho, Sigma = r**2 + epsilon, r**2 + epsilon * x**2
    lam, gx, gp = Sigma / rho, Sigma / d, rho * d
    metric = sp.diag(lam, gx, gp)
    tangents, spatial_labels, shifts0, shift_labels = [], [], [], []

    for ell in range(m, L + 1):
        y, polar, axial, electric, magnetic = harmonic_amplitudes(ell, m, x)
        tangents.extend([sp.diag(y, 0, 0), sp.diag(0, y / d, d * y)])
        spatial_labels.extend([("lambda", ell, m), ("X", ell, m)])

        radial = sp.zeros(3)
        radial[0, 0] = -y * sp.diff(lam, r) / lam
        radial[0, 1] = radial[1, 0] = sp.diff(y, x) - y * sp.diff(lam, x) / lam
        radial[0, 2] = radial[2, 0] = sp.I * m * y
        radial[1, 1] = y * sp.diff(gx, r) / lam
        radial[2, 2] = y * sp.diff(gp, r) / lam
        shifts0.append(radial)
        shift_labels.append(("b", ell, m))

        if ell >= 1:
            for parity, vector in (("polar", polar), ("axial", axial)):
                tangent = sp.zeros(3)
                tangent[0, 1] = tangent[1, 0] = gx * vector[0]
                tangent[0, 2] = tangent[2, 0] = gp * vector[1]
                tangents.append(tangent)
                spatial_labels.append(("V_" + parity, ell, m))

                shift = sp.zeros(3)
                shift[0, 0] = vector[0] * sp.diff(lam, x)
                shift[1, 1] = vector[0] * sp.diff(gx, x) + 2 * gx * sp.diff(vector[0], x)
                shift[2, 2] = vector[0] * sp.diff(gp, x) + 2 * gp * sp.I * m * vector[1]
                shift[1, 2] = shift[2, 1] = gx * sp.I * m * vector[0] + gp * sp.diff(vector[1], x)
                shifts0.append(shift)
                shift_labels.append(("W_" + parity, ell, m))

        if ell >= 2:
            for parity, tensor in (("polar", electric), ("axial", magnetic)):
                tangent = sp.zeros(3)
                tangent[1:, 1:] = tensor
                tangents.append(tangent)
                spatial_labels.append(("C_" + parity, ell, m))

    spatial_index = {label: i for i, label in enumerate(spatial_labels)}
    J1 = sp.zeros(len(tangents), len(shifts0))
    for j, (kind, ell, order) in enumerate(shift_labels):
        target = "lambda" if kind == "b" else kind.replace("W_", "V_", 1)
        J1[spatial_index[(target, ell, order)], j] = 2 if kind == "b" else 1

    return {
        "L": L, "m": m, "r": r, "x": x, "epsilon": epsilon,
        "metric": metric, "inverse_diagonal": [rho / Sigma, d / Sigma, 1 / (rho * d)],
        "volume_density": Sigma,
        "spatial_tangents": tangents, "spatial_labels": spatial_labels,
        "shift_columns_without_profile_derivatives": shifts0,
        "shift_labels": shift_labels, "shift_derivative_injection": J1,
        "lapse_labels": [("a", ell, m) for ell in range(m, L + 1)],
        "angular_integration": "one half integral dx from -1 to 1, after equal-m phi average",
        "implementation_status": "prepared_not_scientifically_executed",
    }


def dewitt_density(left, right, data):
    """Hermitian equal-m density, including sqrt(det sigma) but not dx/2.

    Inputs must be the stated real-chart amplitudes with explicit factors
    of I; changing I's sign conjugates them on -1<x<1. This is not a generic
    conjugation convention for arbitrary complex symbolic parameters.
    """
    inverse = data["inverse_diagonal"]
    conjugate_left = left.applyfunc(lambda value: value.subs(sp.I, -sp.I))
    contraction = sum(
        conjugate_left[i, j] * right[i, j] * inverse[i] * inverse[j]
        for i in range(3) for j in range(3)
    )
    left_trace = sum(conjugate_left[i, i] * inverse[i] for i in range(3))
    right_trace = sum(right[i, i] * inverse[i] for i in range(3))
    return data["volume_density"] * (contraction - left_trace * right_trace)


def angular_density_blocks(data):
    """Densities to integrate completely; no intermediate angular cutoff."""
    T = data["spatial_tangents"]
    U = data["shift_columns_without_profile_derivatives"]
    return {
        "S_density": sp.Matrix([[dewitt_density(a, b, data) for b in T] for a in T]),
        "C0_density": sp.Matrix([[dewitt_density(a, b, data) for b in U] for a in T]),
        "U0_norm_density": sp.Matrix([[dewitt_density(a, b, data) for b in U] for a in U]),
    }


def canonical_angular_blocks(S, C0, U0_norm, J1):
    """Assemble after complete angular integration, before radial projection.

    J(r,d_r)=J0+J1*d_r. S must be invertible and Hermitian in the caller's
    declared representation. These conditions and all integration error
    bounds must be verified separately; no rank or positivity is inferred.
    """
    S, C0, U0_norm, J1 = map(sp.Matrix, (S, C0, U0_norm, J1))
    nq, ns = S.rows, C0.cols
    if S.shape != (nq, nq) or C0.rows != nq or U0_norm.shape != (ns, ns) or J1.shape != (nq, ns):
        raise ValueError("Incompatible integrated block shapes")
    Sinv = S.inv()
    return {
        "K_pp": 2 * Sinv,
        "J0": Sinv * C0,
        "J1": J1,
        "A": (C0.conjugate().T * Sinv * C0 - U0_norm) / 2,
        "missing_blocks": ["V_second_variation", "lapse_spatial_operator"],
        "radial_projection_performed": False,
    }


def realify_block(matrix, m):
    """Positive-sine real basis: sqrt(2)*Re(F), sqrt(2)*Im(F) for m>0."""
    matrix = sp.Matrix(matrix)
    real = matrix.applyfunc(sp.re)
    imag = matrix.applyfunc(sp.im)
    if m == 0:
        if any(sp.simplify(value) != 0 for value in imag):
            raise ValueError("m=0 block must be real; cannot discard an imaginary part")
        return real
    if m < 0:
        raise ValueError("Use nonnegative-m sectors")
    return real.row_join(imag).col_join((-imag).row_join(real))
