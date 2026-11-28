"""Full fifth-order reduced Hamiltonian as exact angular Green pairings.

Inputs are raw rational TT tensors Q,S with Gaussian exp(-r^2/2), where S
is the dual-Gram momentum direction, not the canonical coefficient vector.
The expressions retain q^5, q^3 p^2 and q p^4. No field construction,
coefficient evaluation or quadrature is performed on import.

Only first scalar and York inverses enter this adjoint form of H5. The
same prepared fields supply independent Green evaluation of H3/H4 against
the earlier exact Fourier results. H6 and four-point scattering are outside
this module's scope.
"""

from fractions import Fraction as F

from reduced_gaunt_exact import add, scale, tensor_square
from reduced_gaunt_vertices import (
    matrix_add, matrix_multiply, trace, pg,
)
from higher_local_geometry import local_geometry, york_data
from gaussian_green_fields import scalar_field, york_field
from gaussian_heat_pairings import PairingAccumulator


def prepare_fields(q, s, geometry=None, progress=None):
    """Construct full local polynomials and first inverses for one raw direction."""
    if geometry is None:
        geometry = local_geometry(q, local_R_order=3, density_order=5,
                                  york_order=0, progress=progress)
    if geometry.get("q") != q or geometry.get("local_R_order", 0) < 3 or geometry.get("density_order", 0) < 5:
        raise ValueError("H5 requires the same q geometry through pointwise R3 and density L5")
    york = york_data(q, s, max_order=2, geometry=geometry, progress=progress)
    s2 = tensor_square(s)
    r2 = geometry["R"][2]
    u = scalar_field(add(s2, r2, -1), F(1))
    v = scalar_field(add(r2, s2, 7), F(1))
    w = york_field(tuple(york["b"][1]), F(1))
    dw = [[w[k].derivative(a) for k in range(3)] for a in range(3)]
    zero = scalar_field({}, F(1))
    divergence = zero
    for k in range(3):
        divergence = divergence.add(dw[k][k])
    t0 = [[dw[i][j].add(dw[j][i]).add(divergence, -F(2, 3) if i == j else F(0))
           for j in range(3)] for i in range(3)]
    return {"q": q, "s": s, "geometry": geometry, "york": york,
            "s2": s2, "u": u, "v": v, "w": w, "dw": dw,
            "divergence_w": divergence, "t0": t0,
            "grad_u": u.gradient(), "grad_v": v.gradient()}


def _add(accumulator, polynomial, beta, fields=(), factor=F(1)):
    """Skip exact zero factors before any angular polynomial convolution."""
    if polynomial and factor and all(field.terms for field in fields):
        accumulator.add(polynomial, beta, fields=fields, factor=F(factor))


def third_fourth_families(data, max_pair_products=20_000_000):
    """Independent spatial Green representation of the previously known H3/H4."""
    def accumulator():
        return PairingAccumulator(max_pair_products=max_pair_products)
    q, s, geometry = data["q"], data["s"], data["geometry"]
    out = {name: accumulator() for name in
           ("H3_curvature", "H4_local", "H4_york", "H4_scalar")}
    _add(out["H3_curvature"], geometry["L"][3], F(3, 2), factor=-1)
    _add(out["H4_local"], geometry["L"][4], F(2), factor=-1)
    commutator = matrix_add(matrix_multiply(q, s), matrix_multiply(s, q), -1)
    raw = trace(matrix_multiply(commutator, commutator))
    _add(out["H4_local"], raw, F(2), factor=F(1, 12))
    for i in range(3):
        _add(out["H4_york"], data["york"]["b"][1][i], F(1), (data["w"][i],), -2)
    b2_scalar = add(geometry["R"][2], data["s2"], 7)
    _add(out["H4_scalar"], b2_scalar, F(1), (data["u"],), -1)
    return out


def fifth_families(data, max_pair_products=20_000_000, progress=None):
    """Assemble all six signed H5 families with exact angular contractions.

    I5=N3=-4<b2,W>+2<W,A1W>. The weak form gives
    2<W,A1W>=-2 int T0:T1+2 int Tr(Q T0^2).
    The scalar part is int R3(v-u)+2 int (S:T0)(7u+v)
    -8 int grad(v) Q grad(u), with u=C(s2-R2),v=C(R2+7s2).
    """
    def accumulator():
        return PairingAccumulator(max_pair_products=max_pair_products)
    names = ("curvature", "york_source", "york_weak", "scalar_curvature",
             "scalar_momentum", "scalar_gradient")
    out = {name: accumulator() for name in names}
    q, s, geometry = data["q"], data["s"], data["geometry"]
    w, dw, t0 = data["w"], data["dw"], data["t0"]
    def done(name):
        if progress:
            progress({"family": name, "completed": names.index(name)+1, "of": len(names)})
    _add(out["curvature"], geometry["L"][5], F(5, 2), factor=-1)
    done("curvature")
    for i in range(3):
        _add(out["york_source"], data["york"]["b"][2][i], F(3, 2), (w[i],), -4)
    done("york_source")
    weak = out["york_weak"]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                # -2 T0_ij (L1 W)_ij, keeping W contravariant.
                _add(weak, pg(q[i][j], k, 1), F(1, 2), (t0[i][j], w[k]), -2)
                _add(weak, q[i][k], F(1, 2), (t0[i][j], dw[j][k]), -2)
                _add(weak, q[j][k], F(1, 2), (t0[i][j], dw[i][k]), -2)
                # +2 Tr(Q T0^2); no extra Gaussian from a Green field.
                _add(weak, q[i][k], F(1, 2), (t0[i][j], t0[k][j]), 2)
            _add(weak, q[i][j], F(1, 2), (t0[i][j], data["divergence_w"]), F(4, 3))
    done("york_weak")
    _add(out["scalar_curvature"], geometry["R"][3], F(3, 2),
         (data["v"].add(data["u"], -1),))
    done("scalar_curvature")
    scalar_momentum = data["u"].scale(7).add(data["v"])
    for i in range(3):
        for j in range(3):
            _add(out["scalar_momentum"], s[i][j], F(1, 2), (t0[i][j], scalar_momentum), 2)
    done("scalar_momentum")
    for i in range(3):
        for j in range(3):
            _add(out["scalar_gradient"], q[i][j], F(1, 2),
                 (data["grad_v"][i], data["grad_u"][j]), -8)
    done("scalar_gradient")
    return out


def merge_families(families, max_pair_products=20_000_000):
    """Merge already collapsed exact bins, without repeating Cartesian products."""
    total = PairingAccumulator(max_pair_products=max_pair_products)
    for family in families.values():
        total.merge(family)
    return total
