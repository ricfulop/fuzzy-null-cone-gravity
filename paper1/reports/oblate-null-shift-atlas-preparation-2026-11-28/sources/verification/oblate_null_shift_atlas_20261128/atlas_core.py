"""Exact signed-Fourier amplitudes in the chart x=cos(theta), phi.

One amplitude a in mode m represents (sqrt(1-x*x))**(abs(m)%2) * a * exp(i*m*phi).
Only rational functions enter derivatives, products, round projections and equality checks.
"""
from dataclasses import dataclass
from functools import lru_cache
from typing import NamedTuple

import sympy as s


x, r, epsilon = s.symbols("x r epsilon", real=True)
w, wp = s.symbols("w w_prime", real=True)
d = 1 - x*x
round_metric = s.diag(1/d, d)
inv_round = s.diag(d, 1/d)
ZERO2 = s.ImmutableMatrix.zeros(2, 2)
SCALAR = ("a", "b", "lambda", "X")
VECTOR = ("V_polar", "V_axial", "W_polar", "W_axial")
TENSOR = ("C_polar", "C_axial")


class Label(NamedTuple):
    field: str
    ell: int
    m: int


def labels(lower, upper):
    return [Label(kind, ell, m) for ell in range(lower, upper+1)
            for m in range(-ell, ell+1)
            for kind in SCALAR + (VECTOR if ell >= 1 else ()) + (TENSOR if ell >= 2 else ())]


@dataclass(frozen=True)
class Fields:
    base: s.MatrixBase     # symmetric (t,r) block of independent coefficients a,b,lambda
    shifts: s.MatrixBase   # columns W,V in that order
    gamma: s.MatrixBase    # full variation X*q_round+C


def zero_fields():
    return Fields(ZERO2, ZERO2, ZERO2)


def scaled(f, scale):
    return Fields(scale*f.base, scale*f.shifts, scale*f.gamma)


def add(a, b):
    return Fields(a.base+b.base, a.shifts+b.shifts, a.gamma+b.gamma)


def insist_zero(mat, description):
    for value in mat:
        if value != 0 and s.cancel(value) != 0:
            raise AssertionError(description)


def dx(value, m):
    return value.diff(x) - ((abs(m) % 2)*x/d)*value


def product_factor(k, m):
    return d if k % 2 and m % 2 else s.S.One


@lru_cache(None)
def harmonics(ell, m):
    # |m| for both signs: negative modes are conjugates, with no factorial reweighting.
    n = abs(m)
    if n > ell:
        raise ValueError("invalid harmonic label")
    y = (-1)**n * d**(n//2) * s.diff(s.legendre(ell, x), x, n)
    yp = dx(y, m)
    im = s.I*m
    polar = s.ImmutableMatrix([d*yp, im*y/d])
    axial = s.ImmutableMatrix([im*y, -yp])
    exx = dx(yp, m)-x*yp/d+ell*(ell+1)*y/(2*d)
    exp = im*(yp+x*y/d)
    E = s.ImmutableMatrix([[exx, exp], [exp, -d*d*exx]])
    B = s.ImmutableMatrix([[exp/d, -d*exx], [-d*exx, -d*exp]])
    return y, polar.applyfunc(s.cancel), axial.applyfunc(s.cancel), E.applyfunc(s.cancel), B.applyfunc(s.cancel)


@lru_cache(None)
def unit(label):
    kind, ell, m = label
    y, polar, axial, E, B = harmonics(ell, m)
    base, shifts, gamma = s.zeros(2), s.zeros(2), s.zeros(2)
    if kind == "a":
        base[0, 0] = y
    elif kind == "b":
        base[0, 1] = base[1, 0] = y
    elif kind == "lambda":
        base[1, 1] = y
    elif kind == "X":
        gamma = y*round_metric
    elif kind in VECTOR:
        shifts[:, 0 if kind.startswith("W") else 1] = polar if kind.endswith("polar") else axial
    elif kind in TENSOR:
        gamma = E if kind == "C_polar" else B
    else:
        raise ValueError(kind)
    return Fields(base.as_immutable(), shifts.as_immutable(), gamma.as_immutable())


def branches(axis, rotation):
    # n_x=h cos(phi), n_y=h sin(phi), n_z=x; h is implicit for |k|=1.
    profiles = {"x": ((1, s.S.Half), (-1, s.S.Half)),
                "y": ((1, 1/(2*s.I)), (-1, -1/(2*s.I))), "z": ((0, x),)}
    for k, n in profiles[axis]:
        vector = [-s.I*k*n, dx(n, k)] if rotation else [0, n]
        yield k, s.ImmutableMatrix(vector)


def derivative_matrix(vector, k):
    return s.Matrix.hstack(dx(vector, k), s.I*k*vector)


def advect(vector, mat, m):
    return vector[0]*dx(mat, m) + vector[1]*s.I*m*mat


def generator(f, m, vector, k):
    """Return one Fourier branch of DF for unit profile (the product carries m+k)."""
    scale = product_factor(k, m)
    jac = derivative_matrix(vector, k)
    return Fields(scale*advect(vector, f.base, m),
                  scale*(advect(vector, f.shifts, m)-jac*f.shifts),
                  scale*(advect(vector, f.gamma, m)+jac.T*f.gamma+f.gamma*jac))


@lru_cache(None)
def background():
    rho = r*r+epsilon
    sigma = r*r+epsilon*x*x
    return Fields(s.diag(-1, sigma/rho), ZERO2, s.diag(sigma/d, rho*d))


@lru_cache(None)
def affine(vector, k):
    f = scaled(generator(background(), 0, vector, k), w)
    shift = s.Matrix.hstack(s.zeros(2, 1), wp*vector)
    return Fields(f.base, f.shifts+shift, f.gamma)


def average_polynomial(value):
    # Same-mode Hermitian contractions should cancel to polynomials before integration.
    # Poly deliberately raises on a surviving x-dependent denominator.
    poly = s.Poly(s.cancel(value), x)
    return s.cancel(sum(a/(j+1) for (j,), a in poly.terms() if not j % 2))


def projection(raw, m, ell_max=4):
    scalar = {"a": raw.base[0, 0], "b": raw.base[0, 1], "lambda": raw.base[1, 1],
              "X": s.cancel((inv_round*raw.gamma).trace()/2)}
    c = raw.gamma-scalar["X"]*round_metric
    result = {}
    weight = d**(abs(m) % 2)
    for ell in range(abs(m), ell_max+1):
        y, polar, axial, E, B = harmonics(ell, m)
        norm = s.factorial(ell+abs(m))/(s.factorial(ell-abs(m))*(2*ell+1))

        def put(kind, pairing, denom):
            if pairing != 0:
                coeff = s.cancel(average_polynomial(weight*pairing)/denom)
                if coeff != 0:
                    result[Label(kind, ell, m)] = coeff

        for kind, value in scalar.items():
            put(kind, y*value, norm)
        if ell >= 1:
            for name, col in (("W", 0), ("V", 1)):
                v = raw.shifts[:, col]
                for parity, basis in (("polar", polar), ("axial", axial)):
                    put(name+"_"+parity, (basis.conjugate().T*round_metric*v)[0], ell*(ell+1)*norm)
        if ell >= 2:
            tensor_norm = (ell-1)*ell*(ell+1)*(ell+2)*norm/2
            for parity, basis in (("polar", E), ("axial", B)):
                pairing = (basis.conjugate()*inv_round*c*inv_round).trace()
                put("C_"+parity, pairing, tensor_norm)
    return result


def reconstruct_coefficients(raw, projected):
    fitted = zero_fields()
    for label, amplitude in projected.items():
        fitted = add(fitted, scaled(unit(label), amplitude))
    for attr in ("base", "shifts", "gamma"):
        insist_zero(getattr(raw, attr)-getattr(fitted, attr), "round harmonic reconstruction: "+attr)


def blocks(base, mixed, gamma):
    return base.row_join(mixed.T).col_join(mixed.row_join(gamma))


def metric_tangent(f):
    return blocks(f.base, background().gamma*f.shifts, f.gamma)


def metric_check(f, m, vector, k, df):
    """Direct Lie derivative of h versus DR[df]+D²R[f,F(background)]."""
    h = metric_tangent(f)
    # Independent four-coordinate Jacobian; its radial column uses the profile slope.
    xi = s.Matrix([0, 0, *vector])
    jac4 = s.Matrix.hstack(s.zeros(4, 1), wp*xi, w*dx(xi, k), w*s.I*k*xi)
    product = product_factor(k, m)
    direct = product*(w*vector[0]*dx(h, m)+w*vector[1]*s.I*m*h+jac4.T*h+h*jac4)
    F = affine(vector, k)
    mixed = f.gamma*F.shifts+F.gamma*f.shifts
    cross = f.shifts.T*background().gamma*F.shifts
    second = blocks(cross+cross.T, mixed, ZERO2)
    insist_zero(direct-w*metric_tangent(df)-product*second, "differentiated metric reconstruction")


def column_map(label, axis, rotation=False):
    all_projected = {}
    for k, vector in branches(axis, rotation):
        df = generator(unit(label), label.m, vector, k)
        projected = projection(df, label.m+k)
        reconstruct_coefficients(df, projected)
        metric_check(unit(label), label.m, vector, k, df)
        all_projected.update(projected)  # distinct k yield distinct m+k
    return {row: amplitude for row, amplitude in all_projected.items() if row.ell > 2}


def affine_check(axis, rotation):
    for k, vector in branches(axis, rotation):
        f = affine(vector, k)
        projected = projection(f, k)
        reconstruct_coefficients(f, projected)
        if any(row.ell > 2 for row in projected):
            raise AssertionError("background generator outside retained band")
        insist_zero(f.shifts[:, 1]-wp*vector, "background radial shift slope")
        insist_zero(f.shifts[:, 0], "background temporal shift")
