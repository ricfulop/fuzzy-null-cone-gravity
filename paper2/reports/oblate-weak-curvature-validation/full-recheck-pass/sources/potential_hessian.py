"""Definition-only scaffold for the oblate spatial-curvature Hessian.

This module records the tensor formulas and coefficient-export shape for the
unit-lapse potential

    V[sigma] = -int sqrt(sigma) R[sigma] dr dOmega / (4*pi).

It is deliberately import-safe: no finite-difference check, quadrature scan,
or gravitational operator assembly is run at import time. The full assembly
must supply the retained basis tensors, their r-dependence, the radial domain,
and boundary convention before evaluating the planned checks below.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Tuple


IndexTriple = Tuple[int, int, int]


@dataclass(frozen=True)
class DiagonalOblateMetric:
    """Symbolic names for the exact flat oblate spatial background."""

    g_rr: str = "Sigma/rho"
    g_xx: str = "Sigma/d"
    g_pp: str = "rho*d"
    rho: str = "r**2 + epsilon"
    sigma: str = "r**2 + epsilon*x**2"
    d: str = "1 - x**2"


@dataclass(frozen=True)
class WeakBilinearExport:
    """A first-derivative radial bilinear for one angular m block.

    The coefficient arrays are indexed by retained spatial coefficient labels.
    They represent

        int dr A_IJ f_I' g_J'
             + B_IJ f_I' g_J
             + conj(B_JI) f_I g_J'
             + C_IJ f_I g_J

    plus a separately exported radial boundary form. Entries should be exact
    symbolic expressions when used as rank evidence; numerical quadrature
    exports must carry their precision and convergence metadata elsewhere.
    """

    cutoff_l: int
    m: int
    labels: Tuple[str, ...]
    a_rr: Tuple[Tuple[str, ...], ...]
    b_r0: Tuple[Tuple[str, ...], ...]
    c_00: Tuple[Tuple[str, ...], ...]
    radial_boundary: Tuple[Tuple[str, ...], ...]


CHRISTOFFEL: Mapping[IndexTriple, str] = {
    (0, 0, 0): "r/Sigma - r/rho",
    (0, 0, 1): "epsilon*x/Sigma",
    (0, 1, 0): "epsilon*x/Sigma",
    (0, 1, 1): "-rho*r/(d*Sigma)",
    (0, 2, 2): "-rho*r*d/Sigma",
    (1, 0, 0): "-epsilon*x*d/(rho*Sigma)",
    (1, 0, 1): "r/Sigma",
    (1, 1, 0): "r/Sigma",
    (1, 1, 1): "epsilon*x/Sigma + x/d",
    (1, 2, 2): "rho*x*d/Sigma",
    (2, 0, 2): "r/rho",
    (2, 2, 0): "r/rho",
    (2, 1, 2): "-x/d",
    (2, 2, 1): "-x/d",
}


LINEARIZED_EINSTEIN_TERMS: Tuple[str, ...] = (
    "-1/2 * nabla^2 h_ij",
    "-1/2 * nabla_i nabla_j tr(h)",
    " 1/2 * nabla_i nabla^a h_aj",
    " 1/2 * nabla_j nabla^a h_ai",
    "-1/2 * sigma_ij * nabla^a nabla^b h_ab",
    " 1/2 * sigma_ij * nabla^2 tr(h)",
)


WEAK_BILINEAR_TERMS: Tuple[str, ...] = (
    " 1/2 * nabla_a h_ij * nabla^a k^ij",
    "-nabla_i h^ij * nabla^a k_aj",
    " 1/2 * nabla_i h^ij * nabla_j tr(k)",
    " 1/2 * nabla_i k^ij * nabla_j tr(h)",
    "-1/2 * nabla_a tr(h) * nabla^a tr(k)",
)


PLANNED_CONTROLS: Tuple[str, ...] = (
    "weak Hermiticity up to the exported radial boundary remainder",
    "zero bulk Hessian on compactly supported continuum pure gradients",
    "V'' J_kernel = 0 on exactly retained oblate shift-kernel sections",
    "positive Cartesian TT gradient sign and 1/2 normalization",
    "finite-difference agreement with -int sqrt(sigma) R",
    "separate angular-quadrature and radial-Galerkin convergence",
)


def covariant_derivative_recipe() -> str:
    """Return the coefficient-basis derivative rule used by the export."""

    return (
        "nabla_a h_ij = "
        "(partial_a T_I,ij - Gamma^b_ai T_I,bj - Gamma^b_aj T_I,ib) u^I "
        "+ T_I,ij partial_a u^I"
    )


def validate_export_shape(block: WeakBilinearExport) -> None:
    """Validate dimensions of an already-built block.

    This catches interface errors only. It does not validate any scientific
    matrix element and intentionally performs no assembly.
    """

    n = len(block.labels)
    for name, matrix in (
        ("a_rr", block.a_rr),
        ("b_r0", block.b_r0),
        ("c_00", block.c_00),
        ("radial_boundary", block.radial_boundary),
    ):
        if len(matrix) != n or any(len(row) != n for row in matrix):
            raise ValueError(f"{name} must be a square {n}x{n} matrix")


__all__ = [
    "CHRISTOFFEL",
    "LINEARIZED_EINSTEIN_TERMS",
    "PLANNED_CONTROLS",
    "WEAK_BILINEAR_TERMS",
    "DiagonalOblateMetric",
    "WeakBilinearExport",
    "covariant_derivative_recipe",
    "validate_export_shape",
]
