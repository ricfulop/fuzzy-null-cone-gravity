"""Unexecuted weak-form density assembly for the pinned oblate background.

All functions are preparation. No integrals, rank tests, domain choice, or
scientific checks run at import. Angular integration is dx/2 after equal-m
azimuthal averaging. The metric density Sigma is included exactly once.
"""

import sympy as sp

from kinematic_blocks import harmonic_amplitudes
from potential_hessian import CHRISTOFFEL


def chart_conjugate(value):
    """Known real-chart amplitudes with explicit I; see kinematic_blocks."""
    return value.subs(sp.I, -sp.I)


def connection(data):
    r, x, epsilon = data["r"], data["x"], data["epsilon"]
    names = {"r": r, "x": x, "epsilon": epsilon,
             "rho": r**2 + epsilon, "Sigma": r**2 + epsilon*x**2,
             "d": 1-x**2}
    gamma = [[[sp.S.Zero for _ in range(3)] for _ in range(3)] for _ in range(3)]
    for (upper, lower1, lower2), expression in CHRISTOFFEL.items():
        gamma[upper][lower1][lower2] = sp.sympify(expression, locals=names)
    return gamma


def tensor_jet(tensor, data, value_slot, derivative_slot, gamma=None):
    """Return h and nabla(h) for h=F(r)*tensor with supplied F,F' slots.

    Slots represent a radial jet, not functions to be differentiated again.
    Derivatives of the r-dependent basis tensor are included explicitly.
    """
    if gamma is None:
        gamma = connection(data)
    r, x, m = data["r"], data["x"], data["m"]
    h = value_slot * tensor
    nabla = []
    for a in range(3):
        partial = tensor.diff(r) if a == 0 else tensor.diff(x) if a == 1 else sp.I*m*tensor
        derivative = value_slot * partial
        if a == 0:
            derivative += derivative_slot * tensor
        derivative -= sp.Matrix(3, 3, lambda i, j: sum(
            gamma[b][a][i] * h[b, j] + gamma[b][a][j] * h[i, b]
            for b in range(3)
        ))
        nabla.append(derivative)
    inverse = data["inverse_diagonal"]
    trace = sum(inverse[i] * h[i, i] for i in range(3))
    trace_gradient = [sum(inverse[i] * nabla[a][i, i] for i in range(3)) for a in range(3)]
    divergence_lower = [sum(inverse[i] * nabla[i][i, j] for i in range(3)) for j in range(3)]
    return {"h": h, "nabla": nabla, "trace": trace,
            "trace_gradient": trace_gradient, "divergence_lower": divergence_lower}


def potential_density(left, right, data):
    """The Hermitian first-derivative V'' density; not half of V''."""
    inverse = data["inverse_diagonal"]
    derivative_pair = sum(
        inverse[a] * inverse[i] * inverse[j]
        * chart_conjugate(left["nabla"][a][i, j]) * right["nabla"][a][i, j]
        for a in range(3) for i in range(3) for j in range(3)
    )
    divergence_pair = sum(inverse[j] * chart_conjugate(left["divergence_lower"][j])
                          * right["divergence_lower"][j] for j in range(3))
    cross1 = sum(inverse[j] * chart_conjugate(left["divergence_lower"][j])
                 * right["trace_gradient"][j] for j in range(3))
    cross2 = sum(inverse[j] * chart_conjugate(left["trace_gradient"][j])
                 * right["divergence_lower"][j] for j in range(3))
    trace_pair = sum(inverse[j] * chart_conjugate(left["trace_gradient"][j])
                     * right["trace_gradient"][j] for j in range(3))
    return data["volume_density"] * (derivative_pair/2 - divergence_pair + cross1/2 + cross2/2 - trace_pair/2)


def potential_boundary_flux(left, right, data):
    """Sigma*Theta^r; evaluate its angular integral at the two endpoints."""
    inverse = data["inverse_diagonal"]
    h = left["h"]
    flux = -inverse[0] * sum(
        inverse[i] * inverse[j] * chart_conjugate(h[i, j]) * right["nabla"][0][i, j]
        for i in range(3) for j in range(3)
    ) / 2
    flux += sum(inverse[0] * inverse[j] * chart_conjugate(h[0, j])
                * (right["divergence_lower"][j] - right["trace_gradient"][j]/2)
                for j in range(3))
    flux += chart_conjugate(left["trace"]) * inverse[0] * (
        -right["divergence_lower"][0] + right["trace_gradient"][0]
    ) / 2
    return data["volume_density"] * flux


def lapse_density(y, f_slot, f_prime_slot, right, data):
    """Weak <f,Lq> density; the actual coefficient is delta(a), not delta(alpha)."""
    x, m = data["x"], data["m"]
    gradient_f = [f_prime_slot*y, f_slot*sp.diff(y, x), sp.I*m*f_slot*y]
    return -data["volume_density"] * sum(
        data["inverse_diagonal"][j] * chart_conjugate(gradient_f[j])
        * (right["divergence_lower"][j] - right["trace_gradient"][j])
        for j in range(3)
    ) / 2


def lapse_boundary_flux(y, f_slot, right, data):
    return data["volume_density"] * chart_conjugate(f_slot*y) * data["inverse_diagonal"][0] * (
        right["divergence_lower"][0] - right["trace_gradient"][0]
    ) / 2


def radial_density_blocks(data):
    """Return four derivative-slot matrices for V'' and for L, before integration.

    Key 'ab' means derivative order a on the conjugated left profile and b
    on the right profile. Endpoint fluxes remain separate functions above.
    No symmetry is imposed by replacing an independently assembled block.
    """
    gamma = connection(data)
    tensors = data["spatial_tangents"]
    jets = {
        0: [tensor_jet(tensor, data, 1, 0, gamma) for tensor in tensors],
        1: [tensor_jet(tensor, data, 0, 1, gamma) for tensor in tensors],
    }
    lapse_harmonics = [harmonic_amplitudes(ell, data["m"], data["x"])[0]
                       for _, ell, _ in data["lapse_labels"]]
    potential, lapse = {}, {}
    for a in (0, 1):
        for b in (0, 1):
            key = f"{a}{b}"
            potential[key] = sp.Matrix([[potential_density(left, right, data)
                                         for right in jets[b]] for left in jets[a]])
            lapse[key] = sp.Matrix([[lapse_density(y, int(a == 0), int(a == 1), right, data)
                                     for right in jets[b]] for y in lapse_harmonics])
    return {"potential": potential, "lapse": lapse,
            "angular_integration_performed": False,
            "radial_projection_performed": False,
            "scientific_validation": "not_run"}
