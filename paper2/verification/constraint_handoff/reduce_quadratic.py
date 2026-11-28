"""Exact algebra for a finite quadratic constrained Hamiltonian.

This is a downstream adapter, not a gravitational calculation. Input is the
full quadratic Hamiltonian Hessian at an exact stationary solution, in
canonical order w=(q,p), with nondynamical coordinates z:
    H2 = (w.T K w + 2 z.T E w + z.T D z)/2.
Radial differential operators must first be given an explicit domain and a
compatible finite representation. A frozen-radius matrix is not that export.
Only exact rational matrices are accepted; floating-point rank inference is
deliberately outside this verifier. Run --self-check for manufactured cases.
"""

import argparse
import json
import re
from pathlib import Path

import sympy as s


def columns(vectors, rows):
    return s.Matrix.hstack(*vectors) if vectors else s.zeros(rows, 0)


def null_basis(matrix):
    return columns(matrix.nullspace(), matrix.cols)


def row_basis(matrix):
    return columns(matrix.T.columnspace(), matrix.cols).T


def canonical_form(dimension):
    if dimension % 2:
        raise ValueError("w must contain equally many q and p coordinates")
    n = dimension // 2
    return s.zeros(n).row_join(s.eye(n)).col_join(
        (-s.eye(n)).row_join(s.zeros(n))
    )


def reduce_quadratic(K, E, D):
    dim, nz = K.rows, D.rows
    if K.shape != (dim, dim) or D.shape != (nz, nz) or E.shape != (nz, dim):
        raise ValueError("Expected K=(2d,2d), E=(nz,2d), D=(nz,nz)")
    if K != K.T or D != D.T:
        raise ValueError("K and D must be symmetric")
    if any(not entry.is_Rational for matrix in (K, E, D) for entry in matrix):
        raise ValueError("Use exact rational entries; do not rationalize numerical ranks")
    J = canonical_form(dim)
    Q = columns(D.columnspace(), nz)
    N = null_basis(D)
    A = Q.T * D * Q
    if Q.cols:
        Keff = K - E.T * Q * A.inv() * Q.T * E
    else:
        Keff = K.copy()
    B = N.T * E
    M = E * J * E.T
    bracket = s.zeros(nz).row_join(-D).col_join(D.row_join(M))
    decomposed_rank = 2 * D.rank() + (N.T * M * N).rank()
    assert bracket.rank() == decomposed_rank

    # Keep only the original multiplier columns B.T. Later constraints do not
    # automatically receive new, independent multipliers.
    R = row_basis(B)
    constraint_counts = [R.rows]
    stabilization_steps = []
    for _ in range(dim + 1):
        multiplier_matrix = R * J * B.T
        left_kernel = null_basis(multiplier_matrix.T)
        new_rows = left_kernel.T * R * J * Keff
        updated = row_basis(R.col_join(new_rows))
        stabilization_steps.append({
            "current_rows": R,
            "multiplier_matrix": multiplier_matrix,
            "left_kernel_basis": left_kernel,
            "preservation_rows": new_rows,
            "updated_rows": updated,
            "rank_before": R.rows,
            "rank_after": updated.rows,
        })
        if updated.rows == R.rows:
            break
        R = updated
        constraint_counts.append(R.rows)
    else:
        raise AssertionError("Constraint row rank did not stabilize")

    V = null_basis(R)
    assert left_kernel.T * R * J * Keff * V == s.zeros(left_kernel.cols, V.cols)
    constraint_bracket = R * J * R.T
    gauge = J * R.T * null_basis(constraint_bracket)
    assert R * gauge == s.zeros(R.rows, gauge.cols)
    assert gauge.T * Keff * V == s.zeros(gauge.cols, V.cols)

    # Choose a complement to the null directions of the pulled-back symplectic
    # form. This describes the stabilized linear model only.
    span = gauge
    complement = []
    for vector in V.columnspace():
        candidate = span.row_join(vector)
        if candidate.rank() > span.cols:
            span = candidate
            complement.append(vector)
    T = columns(complement, dim)
    omega = T.T * (-J) * T
    assert omega.rank() == T.cols
    assert T.cols == dim - 2 * R.rows + constraint_bracket.rank()
    Kphysical = T.T * Keff * T
    evolution = omega.inv() * Kphysical if T.cols else s.zeros(0)
    lam = s.Symbol("lambda")
    return {
        "rank_D": D.rank(),
        "rank_initial_primary_secondary_matrix": decomposed_rank,
        "null_secondary_bracket_rank": (N.T * M * N).rank(),
        "constraint_row_counts": constraint_counts,
        "stabilization_steps": stabilization_steps,
        "final_constraint_rows": R,
        "final_constraint_bracket_rank": constraint_bracket.rank(),
        "null_multiplier_source": left_kernel.T * R * J * Keff,
        "gauge_direction_count": gauge.cols,
        "linearized_symplectic_dimension": T.cols,
        "reduced_hamiltonian_hessian": Kphysical,
        "reduced_symplectic_form": omega,
        "evolution_characteristic_polynomial": s.factor(evolution.charpoly(lam).as_expr()),
        "scope": "Exact finite linear model only; not nonlinear Dirac closure, a gravitational mode count, or a quantum measure",
    }


def self_check():
    lam = s.Symbol("lambda")
    passed = []
    # Ordinary oscillator; verifies signs and canonical ordering.
    r = reduce_quadratic(s.diag(9, 1), s.zeros(0, 2), s.zeros(0))
    assert r["evolution_characteristic_polynomial"] == lam**2 + 9
    passed.append("unconstrained oscillator has frequencies +/-3i")

    # Genuine algebraic elimination changes the frequency; dropping it fails.
    r = reduce_quadratic(s.diag(5, 1), s.Matrix([[2, 0]]), s.Matrix([[2]]))
    assert r["evolution_characteristic_polynomial"] == lam**2 + 3
    assert r["evolution_characteristic_polynomial"] != lam**2 + 5
    passed.append("auxiliary Schur correction changes the physical frequency")

    # Indefinite D can have zero diagonals without being singular.
    D = s.Matrix([[0, 1], [1, 0]])
    r = reduce_quadratic(s.diag(5, 1), s.zeros(2, 2), D)
    assert r["rank_D"] == 2
    passed.append("off-diagonal indefinite auxiliary block is invertible")

    # A nonorthogonal image basis requires coordinate congruence, not an
    # orthogonal-projector shortcut. Hand elimination gives diag(5,2).
    S = s.Matrix([[1, 1, 0], [0, 1, 1], [0, 0, 1]])
    D = S.T * s.diag(2, 3, 0) * S
    E = S.T * s.Matrix([[2, 0], [0, 3], [0, 0]])
    Q = D.columnspace()
    assert Q[0].dot(Q[1]) != 0
    r = reduce_quadratic(s.diag(7, 5), E, D)
    assert r["evolution_characteristic_polynomial"] == lam**2 + 10
    passed.append("nonorthogonal image basis agrees with hand Schur elimination")

    # q1=0 generates p1=0 on preservation: a Hessian-only count is wrong.
    K = s.diag(1, 9, 1, 1)
    r = reduce_quadratic(K, s.Matrix([[1, 0, 0, 0]]), s.zeros(1))
    assert r["constraint_row_counts"] == [1, 2]
    assert r["linearized_symplectic_dimension"] == 2
    assert r["evolution_characteristic_polynomial"] == lam**2 + 9
    passed.append("null-secondary preservation finds the missing momentum constraint")

    # Stopping at the initial secondary surface is invalid even if a naive
    # dimension count happens to agree with the final quotient dimension.
    first_step = r["stabilization_steps"][0]
    residual = first_step["preservation_rows"] * null_basis(first_step["current_rows"])
    assert residual != s.zeros(*residual.shape)
    assert r["stabilization_steps"][-1]["rank_before"] == 2
    assert r["stabilization_steps"][-1]["rank_after"] == 2
    passed.append("premature chain termination leaves a nonzero compatibility source")

    # p1 generates p2; both are first class, leaving only q3,p3.
    K = s.zeros(6)
    K[0, 4] = K[4, 0] = 1
    K[2, 2], K[5, 5] = 4, 1
    r = reduce_quadratic(K, s.Matrix([[0, 0, 0, 1, 0, 0]]), s.zeros(1))
    assert r["constraint_row_counts"] == [1, 2]
    assert r["gauge_direction_count"] == 2
    assert r["evolution_characteristic_polynomial"] == lam**2 + 4
    passed.append("first-class chain removes two gauge directions")

    # D=0 does not mean all constraints are first class.
    K = s.diag(1, 4, 1, 1)
    E = s.Matrix([[1, 0, 0, 0], [0, 0, 1, 0]])
    r = reduce_quadratic(K, E, s.zeros(2))
    assert r["null_secondary_bracket_rank"] == 2
    assert r["rank_initial_primary_secondary_matrix"] == 2
    assert r["gauge_direction_count"] == 0
    assert r["evolution_characteristic_polynomial"] == lam**2 + 4
    passed.append("nonzero null-secondary bracket survives a zero temporal Hessian")

    # Recombine auxiliary and null temporal directions, including a nonzero
    # off-null secondary bracket. Physical eigenvalues must be unchanged.
    K = s.diag(5, 9, 1, 1)
    E = s.Matrix([[0, 0, 1, 1], [1, 0, 0, 0]])
    D = s.diag(2, 0)
    base = reduce_quadratic(K, E, D)
    S = s.Matrix([[1, 2], [3, 1]])
    moved = reduce_quadratic(K, S.T * E, S.T * D * S)
    for key in ("rank_initial_primary_secondary_matrix", "linearized_symplectic_dimension", "evolution_characteristic_polynomial"):
        assert base[key] == moved[key], key
    passed.append("mixed temporal-coordinate changes preserve the stabilized spectrum")
    return {"status": "passed", "checks": passed, "sympy": s.__version__,
            "scope": "Manufactured exact algebra checks; no gravitational experiment or physical rank measurement"}


def json_value(value):
    if isinstance(value, s.MatrixBase):
        return {"shape": list(value.shape), "entries": [[str(x) for x in row] for row in value.tolist()]}
    if isinstance(value, s.Basic):
        return str(value)
    raise TypeError(type(value).__name__)


def rational_matrix(data):
    rows, cols = data["shape"]
    entries = data["entries"]
    if len(entries) != rows or any(len(row) != cols for row in entries):
        raise ValueError("Matrix shape does not match entries")
    values = [x for row in entries for x in row]
    if any(isinstance(x, bool) or not (
        isinstance(x, int) or isinstance(x, str) and
        re.fullmatch(r"[+-]?\d+(?:/[1-9]\d*)?", x)
    ) for x in values):
        raise ValueError("Use integers or explicit fractions; decimals must not be promoted to exact rank evidence")
    return s.Matrix(rows, cols, [s.Rational(x) for x in values])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--self-check", action="store_true")
    modes.add_argument("--input", type=Path, help="Exact finite quadratic-model JSON with K,E,D; no numeric rank inference")
    args = parser.parse_args()
    if args.self_check:
        result = self_check()
    else:
        data = json.loads(args.input.read_text())
        result = reduce_quadratic(*(rational_matrix(data[key]) for key in ("K", "E", "D")))
    print(json.dumps(result, default=json_value, indent=2))
