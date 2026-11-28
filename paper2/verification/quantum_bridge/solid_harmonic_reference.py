"""Independent exact real solid harmonics and regular Cartesian shell jets.

Internal coefficients multiply UNNORMALIZED real H_lm(n). For m>0 use
the cosine (real) Condon--Shortley solid; m<0 is its sine (imaginary)
partner with |m|. H_00=1, H_10=z, H_11=-x, H_1,-1=-y. Every m is retained.
norm_over_4pi=h means int H_lm(n)^2 dOmega=4*pi*h; Y=H/sqrt(4*pi*h).

fischer_project(P,L) -> (immutable {(l,m,k):Fraction}, statistics) represents
P(x)=sum c_lmk*r^(2k)*H_lm(x), with mandatory exact reconstruction.
GaussianProfile supplies physical H-shell coefficients and regular amplitudes
G_lm(s), s=r^2, without projecting/dividing by r^l near the origin.
cartesian_jet consumes G,G_s,G_ss,G_sss in the same unnormalized convention.

Default declarations: Cartesian degree<=32, ell<=16, exact products<=5m,
generated polynomial terms<=100000; numerical angular tables have separate
point/work caps. Exceeded support/work limits raise, never prune.
The small bounded basis/monomial caches are immutable mathematical tables;
reported projection work counts exclude cache construction and CPU overhead.

No production HEAT, Fourier, pairing, tail or inverse module is imported.
No basis, field, projection, quadrature or scientific check runs on import.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from math import cos, exp, factorial, fsum, isfinite, pi, sin, sqrt
from numbers import Rational
from types import MappingProxyType


MAX_ELL = 16
MAX_DEGREE = 32


def _integer(value, name, minimum=0, maximum=None):
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        raise ValueError(name+" is outside its declared integer range")
    return value


def _rational(value, name="coefficient"):
    if isinstance(value, bool) or not isinstance(value, Rational):
        raise TypeError(name+" must be an exact real rational")
    return F(value)


def _mode(ell, m):
    _integer(ell, "ell", maximum=MAX_ELL)
    if type(m) is not int or abs(m) > ell:
        raise ValueError("m must be an integer between -ell and ell")
    return ell, m


def _powers(value):
    value = tuple(value)
    if len(value) != 3 or any(type(p) is not int or p < 0 for p in value):
        raise ValueError("Cartesian powers must be three nonnegative integers")
    return value


def _accumulate(out, key, value):
    value = out.get(key, F(0))+value
    if value:
        out[key] = value
    else:
        out.pop(key, None)


def _poly(poly, max_degree=MAX_DEGREE, max_terms=100_000):
    if not isinstance(poly, Mapping):
        raise TypeError("A polynomial must be a mapping")
    if len(poly) > max_terms:
        raise ValueError("Polynomial exceeds max_terms")
    result = {}
    for powers, coefficient in poly.items():
        powers = _powers(powers)
        if sum(powers) > max_degree:
            raise ValueError("Polynomial exceeds max_degree")
        coefficient = _rational(coefficient)
        if coefficient:
            result[powers] = coefficient
    return result


class _ExactWork:
    def __init__(self, max_products=5_000_000, max_terms=100_000):
        self.cap = _integer(max_products, "max_products", 1)
        self.term_cap = _integer(max_terms, "max_terms", 1)
        self.total = self.peak_terms = 0

    def charge(self, count=1):
        if self.total+count > self.cap:
            raise ValueError("Exact harmonic algebra exceeds max_products")
        self.total += count

    def check(self, mapping):
        self.peak_terms = max(self.peak_terms, len(mapping))
        if len(mapping) > self.term_cap:
            raise ValueError("Exact harmonic algebra exceeds max_terms")

    def statistics(self):
        return {"products_and_scatters": self.total, "peak_map_terms": self.peak_terms,
                "max_products": self.cap, "max_terms": self.term_cap}


def _add(a, b, factor=1, work=None):
    result = dict(a)
    for key, value in b.items():
        if work is not None:
            work.charge()
        _accumulate(result, key, factor*value)
        if work is not None:
            work.check(result)
    return result


def _multiply(a, b, work=None):
    if work is not None:
        work.charge(len(a)*len(b))
    result = {}
    for p, c in a.items():
        for q, d in b.items():
            _accumulate(result, tuple(x+y for x, y in zip(p, q)), c*d)
            if work is not None:
                work.check(result)
    return result


def _derivative(poly, axis):
    result = {}
    for powers, coefficient in poly.items():
        if powers[axis]:
            lower = list(powers)
            lower[axis] -= 1
            result[tuple(lower)] = powers[axis]*coefficient
    return result


def _shift(poly, axis, factor=1):
    result = {}
    for powers, coefficient in poly.items():
        raised = list(powers)
        raised[axis] += 1
        result[tuple(raised)] = factor*coefficient
    return result


@lru_cache(maxsize=33)
def _r2_power(k):
    _integer(k, "r2 power", maximum=MAX_DEGREE//2)
    numerator = factorial(k)
    return MappingProxyType({(2*a, 2*b, 2*(k-a-b)): F(numerator, factorial(a)*factorial(b)*factorial(k-a-b))
                             for a in range(k+1) for b in range(k-a+1)})


@lru_cache(maxsize=512)
def _complex_solid(ell, m):
    if ell == m == 0:
        return MappingProxyType({(0, 0, 0): F(1)}), MappingProxyType({})
    if ell == m:
        re, im = _complex_solid(m-1, m-1)
        factor = -(2*m-1)
        real = _add(_shift(re, 0, factor), _shift(im, 1, -factor))
        imag = _add(_shift(im, 0, factor), _shift(re, 1, factor))
    elif ell == m+1:
        re, im = _complex_solid(m, m)
        real, imag = _shift(re, 2, 2*m+1), _shift(im, 2, 2*m+1)
    else:
        first, second = _complex_solid(ell-1, m), _complex_solid(ell-2, m)
        output = []
        for a, b in zip(first, second):
            numerator = _add(_shift(a, 2, 2*ell-1), _multiply(_r2_power(1), b), -(ell+m-1))
            output.append({p: c/F(ell-m) for p, c in numerator.items()})
        real, imag = output
    return MappingProxyType(real), MappingProxyType(imag)


@dataclass(frozen=True)
class SolidHarmonic:
    ell: int
    m: int
    polynomial: Mapping
    norm_over_4pi: F

    @property
    def poly(self):
        return self.polynomial

    def derivative(self, orders):
        orders = _powers(orders)
        if sum(orders) > 3:
            raise ValueError("Solid-harmonic jet order is capped at three")
        return _solid_derivative(self.ell, self.m, orders)


@lru_cache(maxsize=512)
def solid_harmonic(ell, m):
    _mode(ell, m)
    re, im = _complex_solid(ell, abs(m))
    polynomial = im if m < 0 else re
    h = F(1, 2*ell+1) if m == 0 else F(factorial(ell+abs(m)), 2*(2*ell+1)*factorial(ell-abs(m)))
    return SolidHarmonic(ell, m, polynomial, h)


@lru_cache(maxsize=8192)
def _solid_derivative(ell, m, orders):
    result = dict(solid_harmonic(ell, m).polynomial)
    for axis, count in enumerate(orders):
        for _ in range(count):
            result = _derivative(result, axis)
    return MappingProxyType(result)


def _odd_factorial(n):
    result = 1
    for value in range(n, 0, -2):
        result *= value
    return result


@lru_cache(maxsize=32768)
def _sphere_average(powers):
    if any(power % 2 for power in powers):
        return F(0)
    numerator = 1
    for power in powers:
        numerator *= _odd_factorial(power-1)
    return F(numerator, _odd_factorial(sum(powers)+1))


class AngularSupportError(ValueError):
    def __init__(self, remainder, statistics):
        super().__init__("Nonzero exact Fischer remainder outside declared ell_max; terms="+str(len(remainder)))
        self.remainder = MappingProxyType(dict(remainder))
        self.statistics = statistics


def _fischer(poly, ell_max, work):
    homogeneous = {}
    for powers, coefficient in poly.items():
        homogeneous.setdefault(sum(powers), {})[powers] = coefficient
    coefficients, reconstructed = {}, {}
    for degree, source in sorted(homogeneous.items()):
        for ell in range(degree % 2, min(degree, ell_max)+1, 2):
            k = (degree-ell)//2
            for m in range(-ell, ell+1):
                harmonic = solid_harmonic(ell, m)
                work.charge(len(source)*len(harmonic.polynomial))
                projection = F(0)
                for p, c in source.items():
                    for q, d in harmonic.polynomial.items():
                        projection += c*d*_sphere_average(tuple(a+b for a, b in zip(p, q)))
                projection /= harmonic.norm_over_4pi
                if projection:
                    coefficients[ell, m, k] = projection
                    term = _multiply(_r2_power(k), harmonic.polynomial, work)
                    reconstructed = _add(reconstructed, term, projection, work)
                    work.check(coefficients)
    remainder = _add(poly, reconstructed, -1, work)
    return coefficients, remainder


def fischer_project(poly, ell_max, *, max_products=5_000_000, max_terms=100_000, max_degree=MAX_DEGREE):
    """Project a full declared band and require exact reconstruction, all m.

    Raises AngularSupportError with immutable .remainder and .statistics if
    any polynomial component is missing; never silently truncates that part.
    """
    _integer(ell_max, "ell_max", maximum=MAX_ELL)
    _integer(max_degree, "max_degree", maximum=MAX_DEGREE)
    work = _ExactWork(max_products, max_terms)
    poly = _poly(poly, max_degree, max_terms)
    coefficients, remainder = _fischer(poly, ell_max, work)
    statistics = {"ell_max": ell_max, "source_degree": max((sum(p) for p in poly), default=0),
                  "source_terms": len(poly), "fischer_terms": len(coefficients),
                  "remainder_terms": len(remainder), "all_m": True, "work": work.statistics()}
    if remainder:
        raise AngularSupportError(remainder, statistics)
    return MappingProxyType(coefficients), statistics


def _profile_derivative(poly, beta):
    result = {k: -beta*c for k, c in poly.items()}
    for k, c in poly.items():
        if k:
            _accumulate(result, k-1, k*c)
    return {k: c for k, c in result.items() if c}


def _eval_s_poly(poly, s):
    if not poly:
        return 0.0
    value = 0.0
    for k in range(max(poly), -1, -1):
        value = value*s+float(poly.get(k, 0))
    return value


class GaussianProfile:
    """P(r*n)e^-beta*r² in H modes, with explicit regular-origin source data."""
    def __init__(self, coefficients, beta):
        beta = _rational(beta, "beta")
        if beta <= 0:
            raise ValueError("beta must be strictly positive")
        profiles = {}
        if len(coefficients) > 100_000:
            raise ValueError("Gaussian profile coefficient cap exceeded")
        for (ell, m, k), coefficient in coefficients.items():
            _mode(ell, m)
            _integer(k, "radial trace level")
            if ell+2*k > MAX_DEGREE:
                raise ValueError("Gaussian profile exceeds Cartesian degree cap")
            _accumulate(profiles.setdefault((ell, m), {}), k, _rational(coefficient))
        self.beta = beta
        self.profiles = MappingProxyType({mode: MappingProxyType(poly) for mode, poly in profiles.items() if poly})
        self.modes = tuple(sorted(self.profiles))
        derivatives = {}
        for mode, poly in self.profiles.items():
            rows = [dict(poly)]
            for _ in range(3):
                rows.append(_profile_derivative(rows[-1], beta))
            derivatives[mode] = tuple(MappingProxyType(row) for row in rows)
        self.derivatives = MappingProxyType(derivatives)

    def regular_jet(self, ell, m, radius, order=3):
        _mode(ell, m)
        _integer(order, "order", maximum=3)
        radius = float(radius)
        if not isfinite(radius) or radius < 0 or not isfinite(radius*radius):
            raise ValueError("radius must be nonnegative with finite square")
        if (ell, m) not in self.derivatives:
            return (0.0,)*(order+1)
        s = radius*radius
        envelope = exp(-float(self.beta)*s)
        return tuple(envelope*_eval_s_poly(poly, s) for poly in self.derivatives[ell, m][:order+1])

    def sample(self, radii, max_values=2_000_000):
        radii = tuple(float(radius) for radius in radii)
        _integer(max_values, "max_values", 1)
        if len(radii)*len(self.modes) > max_values:
            raise ValueError("Gaussian shell sampling exceeds max_values")
        if any(not isfinite(r) or r < 0 or not isfinite(r*r) for r in radii):
            raise ValueError("Radii must be nonnegative with finite squares")
        return {mode: [r**mode[0]*self.regular_jet(*mode, r, order=0)[0] for r in radii] for mode in self.modes}

    def origin_data(self, ell, m, degree, patch_radius, order=3):
        """Exact G(s) Taylor coefficients and uniform s-derivative error bounds.

        Uses the Lagrange remainder for exp(-beta*s), whose derivatives on
        s>=0 have magnitude <=beta**n. No source fit or inferred error is used.
        Bounds describe G-G_Taylor and its first `order` s derivatives.
        """
        _mode(ell, m)
        _integer(degree, "Taylor degree", maximum=MAX_DEGREE)
        _integer(order, "order", maximum=3)
        patch_radius = _rational(patch_radius, "patch_radius")
        if patch_radius <= 0:
            raise ValueError("patch_radius must be strictly positive")
        poly = self.profiles.get((ell, m), {})
        coefficients = tuple(sum((c*(-self.beta)**(n-k)/factorial(n-k)
                                  for k, c in poly.items() if k <= n), F(0)) for n in range(degree+1))
        s0, bounds = patch_radius*patch_radius, []
        current = dict(poly)
        for j in range(order+1):
            remainder_degree = degree-j
            bound = F(0)
            for k, c in current.items():
                if k <= remainder_degree:
                    n = remainder_degree-k+1
                    bound += abs(c)*self.beta**n*s0**(remainder_degree+1)/factorial(n)
                else:
                    bound += abs(c)*s0**k
            bounds.append(bound)
            current = _profile_derivative(current, self.beta)
        return {"coefficients": coefficients, "patch_radius": patch_radius,
                "error_bounds_s": tuple(bounds), "convention": "unnormalized H, G(s), s=r²",
                "bound_method": "Exact polynomial times exponential Taylor remainder; no fitted source"}


def normalization_factor(ell, m):
    """H coefficient -> Y coefficient multiplier sqrt(4*pi*h)."""
    return sqrt(4*pi*float(solid_harmonic(ell, m).norm_over_4pi))


def convert_coefficients(coefficients, to="Y"):
    if to not in ("H", "Y"):
        raise ValueError("Target convention must be H or Y")
    return {mode: float(value)*normalization_factor(*mode) if to == "Y"
            else float(value)/normalization_factor(*mode) for mode, value in coefficients.items()}


def _evaluate_poly(poly, point):
    return fsum(float(c)*point[0]**p[0]*point[1]**p[1]*point[2]**p[2] for p, c in poly.items())


def _orders(indices):
    return tuple(indices.count(axis) for axis in range(3))


def cartesian_jet(regular_jets, point, order=3, *, max_work=1_000_000):
    """Return value/gradient/hessian/third; unused higher entries are None.

    The supplied tuple for each (l,m) is (G,G_s,...), with respect to s=r².
    There is no division by r or r^l, including at the exact Cartesian origin.
    """
    _integer(order, "order", maximum=3)
    _integer(max_work, "max_work", 1)
    point = tuple(float(x) for x in point)
    if len(point) != 3 or any(not isfinite(x) for x in point):
        raise ValueError("point must contain three finite coordinates")
    indices = {(): []}
    if order >= 1:
        indices.update({(i,): [] for i in range(3)})
    if order >= 2:
        indices.update({(i, j): [] for i in range(3) for j in range(i, 3)})
    if order >= 3:
        indices.update({(i, j, k): [] for i in range(3) for j in range(i, 3) for k in range(j, 3)})
    work = 0
    for (ell, m), jet in regular_jets.items():
        _mode(ell, m)
        jet = tuple(float(value) for value in jet)
        if len(jet) < order+1 or any(not isfinite(value) for value in jet):
            raise ValueError("Each regular jet needs finite G derivatives through the requested order")
        harmonic = solid_harmonic(ell, m)
        derivatives = {}
        for index in indices:
            poly = harmonic.derivative(_orders(index))
            work += len(poly)
            if work > max_work:
                raise ValueError("Cartesian reconstruction exceeds max_work")
            derivatives[index] = _evaluate_poly(poly, point)
        h = derivatives[()]
        indices[()].append(jet[0]*h)
        if order >= 1:
            for i in range(3):
                indices[i,].append(jet[0]*derivatives[i,]+2*point[i]*jet[1]*h)
        if order >= 2:
            for i in range(3):
                for j in range(i, 3):
                    indices[i, j].append(jet[0]*derivatives[i, j]
                        +2*jet[1]*(point[i]*derivatives[j,]+point[j]*derivatives[i,]+(i == j)*h)
                        +4*point[i]*point[j]*jet[2]*h)
        if order >= 3:
            for i in range(3):
                for j in range(i, 3):
                    for k in range(j, 3):
                        indices[i, j, k].append(
                            jet[0]*derivatives[i, j, k]
                            +2*jet[1]*(point[i]*derivatives[j, k]+point[j]*derivatives[i, k]+point[k]*derivatives[i, j]
                                       +(i == j)*derivatives[k,]+(i == k)*derivatives[j,]+(j == k)*derivatives[i,])
                            +4*jet[2]*(point[i]*point[j]*derivatives[k,]+point[i]*point[k]*derivatives[j,]
                                       +point[j]*point[k]*derivatives[i,]
                                       +((i == j)*point[k]+(i == k)*point[j]+(j == k)*point[i])*h)
                            +8*point[i]*point[j]*point[k]*jet[3]*h)
    values = {key: fsum(rows) for key, rows in indices.items()}
    if any(not isfinite(value) for value in values.values()):
        raise ArithmeticError("Nonfinite Cartesian harmonic jet")
    return {"value": values[()],
            "gradient": tuple(values[i,] for i in range(3)) if order >= 1 else None,
            "hessian": tuple(tuple(values[tuple(sorted((i, j)))] for j in range(3)) for i in range(3)) if order >= 2 else None,
            "third": tuple(tuple(tuple(values[tuple(sorted((i, j, k)))] for k in range(3)) for j in range(3)) for i in range(3)) if order >= 3 else None}


def _legendre(n, x):
    previous, value = 1.0, x
    if n == 1:
        return value, 1.0
    for k in range(2, n+1):
        previous, value = value, ((2*k-1)*x*value-(k-1)*previous)/k
    return value, n*(x*value-previous)/(x*x-1)


@lru_cache(maxsize=64)
def _gauss_legendre(n):
    """Independent binary64 angular nodes; ordinary weights sum to TWO."""
    pairs = []
    for i in range(1, n//2+1):
        root = cos(pi*(i-0.25)/(n+0.5))
        for _ in range(40):
            value, derivative = _legendre(n, root)
            step = value/derivative
            root -= step
            if abs(step) <= 4e-16:
                break
        else:
            raise ArithmeticError("Angular Gauss root did not converge")
        _, derivative = _legendre(n, root)
        weight = 2/((1-root*root)*derivative*derivative)
        pairs.extend(((-root, weight), (root, weight)))
    if n % 2:
        _, derivative = _legendre(n, 0.0)
        pairs.append((0.0, 2/(derivative*derivative)))
    return tuple(sorted(pairs))


class AngularGrid:
    """Full-m H projection/reconstruction on an independent full sphere grid."""
    def __init__(self, ell_max, n_mu, n_phi, *, max_points=10_000, max_work=5_000_000):
        _integer(ell_max, "ell_max", maximum=MAX_ELL)
        _integer(n_mu, "n_mu", ell_max+1, 256)
        _integer(n_phi, "n_phi", 2*ell_max+1, 512)
        _integer(max_points, "max_points", 1)
        _integer(max_work, "max_work", 1)
        if n_mu*n_phi > max_points:
            raise ValueError("Angular grid exceeds max_points")
        self.ell_max, self.n_mu, self.n_phi = ell_max, n_mu, n_phi
        self.modes = tuple((ell, m) for ell in range(ell_max+1) for m in range(-ell, ell+1))
        self.factors = tuple(normalization_factor(*mode) for mode in self.modes)
        harmonics = tuple(solid_harmonic(*mode).polynomial for mode in self.modes)
        self.table_work = n_mu*n_phi*sum(len(poly) for poly in harmonics)
        if self.table_work > max_work:
            raise ValueError("Angular harmonic tables exceed max_work")
        directions, weights = [], []
        for mu, weight in _gauss_legendre(n_mu):
            transverse = sqrt(max(0.0, 1-mu*mu))
            for j in range(n_phi):
                angle = 2*pi*j/n_phi
                directions.append((transverse*cos(angle), transverse*sin(angle), mu))
                weights.append(2*pi*weight/n_phi)
        self.directions, self.weights = tuple(directions), tuple(weights)
        self.normalized_table = tuple(tuple(_evaluate_poly(poly, point)/factor
                                           for poly, factor in zip(harmonics, self.factors)) for point in self.directions)

    def project(self, values, *, max_work=5_000_000):
        values = tuple(float(value) for value in values)
        _integer(max_work, "max_work", 1)
        if len(values) != len(self.directions) or any(not isfinite(value) for value in values):
            raise ValueError("One finite field value per angular point is required")
        if len(values)*len(self.modes) > max_work:
            raise ValueError("Angular projection exceeds max_work")
        return {mode: fsum(weight*value*row[index] for weight, value, row in zip(self.weights, values, self.normalized_table))/self.factors[index]
                for index, mode in enumerate(self.modes)}

    def reconstruct(self, coefficients, *, max_work=5_000_000):
        _integer(max_work, "max_work", 1)
        if any(mode not in self.modes for mode in coefficients):
            raise ValueError("Reconstruction includes a mode outside the declared full band")
        if len(self.directions)*len(self.modes) > max_work:
            raise ValueError("Angular reconstruction exceeds max_work")
        values = tuple(float(coefficients.get(mode, 0))*factor for mode, factor in zip(self.modes, self.factors))
        if any(not isfinite(value) for value in values):
            raise ValueError("Harmonic coefficients must be finite")
        return tuple(fsum(coefficient*y for coefficient, y in zip(values, row)) for row in self.normalized_table)

    def statistics(self):
        return {"ell_max": self.ell_max, "modes": len(self.modes), "all_m": True,
                "n_mu": self.n_mu, "n_phi": self.n_phi, "points": len(self.directions),
                "table_monomial_products": self.table_work,
                "coefficient_convention": "unnormalized H", "sphere_weight_sum": fsum(self.weights),
                "normalization": "Y=H/sqrt(4*pi*norm_over_4pi)", "accuracy_certified": False}
