"""Diagnostic right kernel of the pinned step-02 constraint matrix over QQ(r,T)."""

import argparse
import ast
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import tempfile
import time

import sympy as sp
from sympy.polys.matrices import DomainMatrix


R, T = sp.symbols("r T", positive=True)
FIELD = sp.QQ.frac_field(R, T)
SHAPE = (6, 14)
CHECKPOINT_STATE = "intermediate DAE elimination stage; not a seed certificate"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def report(**fields):
    print("::orx:: " + json.dumps(fields, allow_nan=False), flush=True)


@contextmanager
def stage(name, durations):
    started = time.monotonic()
    report(phase=name, step=0, of=1, unit="check")
    try:
        yield
    except BaseException as error:
        print("PROGRESS " + json.dumps({"stage": name, "event": "error",
              "seconds": time.monotonic() - started, "error": repr(error)}), flush=True)
        raise
    durations[name] = time.monotonic() - started
    report(phase=name, step=1, of=1, unit="check", m={name + "_seconds": durations[name]})


def read_checkpoint(path, expected_digest):
    require(re.fullmatch("[0-9a-f]{64}", expected_digest), "Expected a lowercase checkpoint SHA-256")
    data = path.read_bytes()
    actual = hashlib.sha256(data).hexdigest()
    require(actual == expected_digest, "Checkpoint SHA-256 mismatch")
    obj = json.loads(data)
    require(obj["state"] == CHECKPOINT_STATE and obj["step"] == 2,
            "Expected the intermediate step-02 elimination checkpoint")
    require(obj["total_constraint_count"] == SHAPE[0], "Unexpected checkpoint constraint count")
    rows = obj["constraints"]
    require(isinstance(rows, list) and len(rows) == SHAPE[0], "Wrong constraint row count")
    require(all(isinstance(row, list) and len(row) == SHAPE[1] for row in rows),
            "Wrong constraint column count")
    return rows, actual


def parse_entry(text):
    # These are authenticated reducer expressions. Check the integer syntax
    # before SymPy evaluation can simplify an inexact literal out of existence.
    require(isinstance(text, str), "Constraint entries must be strings")
    tree = ast.parse(text, mode="eval")
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant):
            require(type(node.value) is int, "Expected an integer literal")
        elif isinstance(node, ast.Name):
            require(node.id in ("r", "T"), "Unexpected field symbol")
        else:
            require(isinstance(node, (ast.Expression, ast.Load, ast.UnaryOp, ast.UAdd, ast.USub,
                                      ast.BinOp, ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow)),
                    "Unexpected expression syntax")
    return sp.sympify(text, locals={"r": R, "T": T})


def field_matrix(rows, shape):
    entries = {}
    for i, row in enumerate(rows):
        for j, expr in enumerate(row):
            require(not expr.has(sp.Float), "Inexact coefficient in constraint matrix")
            value = FIELD.from_sympy(expr)
            if value:
                entries.setdefault(i, {})[j] = value
    return DomainMatrix.from_dod(entries, shape, FIELD)


def rref_with_denominator(original):
    if 0 in original.shape:
        return original, FIELD.one, ()
    return original.rref_den(method="CD", keep_domain=False)


def normalize(numerator, denominator, pivots):
    scale = FIELD.convert(denominator, numerator.domain)
    require(bool(scale), "RREF denominator is zero")
    reduced = numerator.convert_to(FIELD) / scale
    rows, columns = reduced.shape
    rank = len(pivots)
    require(reduced.extract(range(rows), pivots) == DomainMatrix.eye((rows, rank), FIELD),
            "Normalized pivot columns do not form the expected identity block")
    require(reduced.extract(range(rank, rows), range(columns)).is_zero_matrix,
            "RREF has a nonzero trailing row")
    return reduced


def kernel_embedding(reduced, pivots):
    columns = reduced.shape[1]
    free = [column for column in range(columns) if column not in pivots]
    entries = {column: {j: FIELD.one} for j, column in enumerate(free)}
    for i, pivot in enumerate(pivots):
        for j, column in enumerate(free):
            value = reduced[i, column].element
            if value:
                entries.setdefault(pivot, {})[j] = -value
    embedding = DomainMatrix.from_dod(entries, (columns, len(free)), FIELD)
    require(embedding.extract(free, range(len(free))) == DomainMatrix.eye(len(free), FIELD),
            "Free-coordinate rows are not the identity")
    return embedding, free


def check_original_product(original, embedding):
    require(original.matmul(embedding).is_zero_matrix, "Original C times candidate M is nonzero")


def manufactured_controls():
    # The middle column is free; its first entry in the kernel is a proper rational function.
    expressions = [["r+1", "(T+1)/(r+2)", "0"], ["0", "0", "1/(T+2)"],
                   ["2*(r+1)", "2*(T+1)/(r+2)", "0"]]
    original = field_matrix([[parse_entry(v) for v in row] for row in expressions], (3, 3))
    numerator, denominator, pivots = rref_with_denominator(original)
    reduced = normalize(numerator, denominator, pivots)
    embedding, free = kernel_embedding(reduced, pivots)
    first = FIELD.from_sympy(-(T+1)/((R+1)*(R+2)))
    expected = DomainMatrix.from_dod({0: {0: first}, 1: {0: FIELD.one}}, (3, 1), FIELD)
    require(tuple(pivots) == (0, 2) and free == [1] and embedding == expected, "Wrong known kernel")
    check_original_product(original, embedding)
    wrong = expected + DomainMatrix.from_dod({0: {0: FIELD.one}}, (3, 1), FIELD)
    try:
        check_original_product(original, wrong)
    except ValueError:
        pass
    else:
        raise AssertionError("The original-matrix check accepted a nonzero rational product")

    for shape in ((0, 3), (2, 0), (0, 0), (2, 3)):
        original = DomainMatrix.zeros(shape, FIELD)
        numerator, denominator, pivots = rref_with_denominator(original)
        embedding, free = kernel_embedding(normalize(numerator, denominator, pivots), pivots)
        require(embedding.shape == (shape[1], shape[1]) and free == list(range(shape[1]))
                and embedding == DomainMatrix.eye(shape[1], FIELD),
                "Wrong zero-matrix kernel or shape")
        check_original_product(original, embedding)
    original = DomainMatrix.eye(2, FIELD)
    numerator, denominator, pivots = rref_with_denominator(original)
    embedding, free = kernel_embedding(normalize(numerator, denominator, pivots), pivots)
    require(not free and embedding.shape == (2, 0), "Wrong full-rank kernel shape")
    check_original_product(original, embedding)
    for text in ("0.5", "0*1.0"):
        try:
            parse_entry(text)
        except ValueError:
            pass
        else:
            raise AssertionError("An inexact literal was accepted")
    return ["known_rational_kernel", "nonzero_original_product_rejected", "empty_dimensions",
            "zero_matrix", "full_column_rank", "inexact_literals_rejected"]


def publish_new(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent, delete=False) as stream:
            temporary = Path(stream.name)
            json.dump(payload, stream, indent=2, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def probe(checkpoint, digest, output):
    require(not output.exists(), "Output already exists")
    durations = {}
    with stage("input_pin_verification", durations):
        texts, actual = read_checkpoint(checkpoint, digest)
    with stage("manufactured_controls", durations):
        controls = manufactured_controls()
    report(m={"manufactured_control_count": len(controls)})
    with stage("checkpoint_expression_parsing", durations):
        exprs = [[parse_entry(value) for value in row] for row in texts]
    with stage("field_conversion", durations):
        original = field_matrix(exprs, SHAPE)
    with stage("constraint_rref", durations):
        numerator, denominator, pivots = rref_with_denominator(original)
    with stage("denominator_normalization", durations):
        reduced = normalize(numerator, denominator, pivots)
    with stage("kernel_embedding", durations):
        embedding, free = kernel_embedding(reduced, pivots)
    with stage("original_constraint_product", durations):
        check_original_product(original, embedding)
    with stage("serialization", durations):
        encoded = [[str(FIELD.to_sympy(value)) for value in row] for row in embedding.to_list()]
        publish_new(output, {
            "state": "diagnostic terminal constraint kernel verified", "checkpoint_sha256": actual,
            "checkpoint_step": 2, "epsilon0": "1/100", "field": "QQ(r,T)", "sympy": sp.__version__,
            "python": sys.version, "executable": sys.executable, "method": "rref_den CD",
            "constraint_shape": list(original.shape), "formal_constraint_rank": len(pivots),
            "pivots": list(pivots), "free_columns": free, "embedding_shape": list(embedding.shape),
            "embedding": encoded, "original_constraint_product_zero": True, "free_rows_identity": True,
            "manufactured_controls_passed": controls, "stage_seconds_before_serialization": dict(durations),
            "terminal_constraint_preservation_checked": False, "original_22_retained_rows_checked": False,
            "reduced_dynamics_reconstructed": False, "discarded_response_checked": False,
            "dae_closure_certified": False, "smooth_solution_germ_certified": False,
        })
    report(m={"diagnostic_kernel_verified": 1, "formal_constraint_rank": len(pivots),
              "free_column_count": len(free)})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--checkpoint-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    probe(args.checkpoint, args.checkpoint_sha256, args.output)
