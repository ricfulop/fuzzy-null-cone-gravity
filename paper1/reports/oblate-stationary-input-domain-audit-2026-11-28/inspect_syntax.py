"""Evidence only. Inspect the pinned JSON syntax without evaluating any coefficient."""
import argparse
import ast
import hashlib
import json
from collections import Counter
from pathlib import Path

OPERATORS = "reports/oblate-stationary-system-recheck-2026-11-27/outputs/stationary-operators.json"
STATE = "reports/oblate-stationary-field-dynamics-timeout-2026-11-28/outputs/terminal-state-before-preservation.json"
OPERATOR_PIN = "07ebdcd4e9fe08fa21e1e72ad5729da05ea7d2fe7f0b04904bb1d489cb985c60"
STATE_PIN = "bc7d3b8a0ffccc771538da1c24cf3bbc66daccdab278a01fcfa4d1ae70d24d08"


def read_pinned(path, pin):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != pin:
        raise ValueError(f"Input hash mismatch: {path.name}")
    return json.loads(raw)


def positive_integer(n):
    return isinstance(n, ast.Constant) and type(n.value) is int and n.value > 0


def positive_exponent(n):
    return positive_integer(n) or (
        isinstance(n, ast.BinOp) and isinstance(n.op, ast.Div)
        and positive_integer(n.left) and positive_integer(n.right)
    )


def positive_on_shell(n):
    if positive_integer(n):
        return True
    if isinstance(n, ast.Name):
        return n.id in ("r", "epsilon")
    if isinstance(n, ast.BinOp):
        if isinstance(n.op, (ast.Add, ast.Mult, ast.Div)):
            return positive_on_shell(n.left) and positive_on_shell(n.right)
        if isinstance(n.op, ast.Pow):
            return positive_on_shell(n.left) and positive_exponent(n.right)
    if isinstance(n, ast.Call):
        return (isinstance(n.func, ast.Name) and n.func.id == "sqrt"
                and len(n.args) == 1 and not n.keywords and positive_on_shell(n.args[0]))
    return False


def polynomial(n):
    """Polynomial syntax in r,T with rational coefficients (integer divisors allowed)."""
    if isinstance(n, ast.Constant):
        return type(n.value) is int
    if isinstance(n, ast.Name):
        return n.id in ("r", "T")
    if isinstance(n, ast.UnaryOp):
        return isinstance(n.op, (ast.UAdd, ast.USub)) and polynomial(n.operand)
    if isinstance(n, ast.BinOp):
        if isinstance(n.op, (ast.Add, ast.Sub, ast.Mult)):
            return polynomial(n.left) and polynomial(n.right)
        if isinstance(n.op, ast.Div):
            return polynomial(n.left) and positive_integer(n.right)
        if isinstance(n.op, ast.Pow):
            return polynomial(n.left) and positive_integer(n.right)
    return False


def inspect(paper):
    ops = read_pinned(paper / OPERATORS, OPERATOR_PIN)
    state = read_pinned(paper / STATE, STATE_PIN)
    rhs = set()
    entries = 0
    for matrix in ops["matrices_by_derivative_order"]:
        for row in matrix:
            for expression in row:
                entries += 1
                for node in ast.walk(ast.parse(expression, mode="eval")):
                    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
                        if not positive_on_shell(node.right):
                            raise ValueError("Unclassified operator divisor")
                        rhs.add(ast.dump(node.right, include_attributes=False))
                    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Pow):
                        if not positive_exponent(node.right):
                            raise ValueError("Unclassified operator exponent")
                    if isinstance(node, ast.Call):
                        if not isinstance(node.func, ast.Name) or node.func.id not in ("sqrt", "atan"):
                            raise ValueError("Unclassified operator function")
                        if node.func.id == "sqrt" and not positive_on_shell(node.args[0]):
                            raise ValueError("Unclassified square-root argument")

    kinds = {}
    for name in ("C", "K", "B", "M"):
        counts = Counter()
        for row in state["matrices"][name]["entries"]:
            for expression in row:
                root = ast.parse(expression, mode="eval").body
                if polynomial(root):
                    counts["polynomial"] += 1
                elif (isinstance(root, ast.BinOp) and isinstance(root.op, ast.Div)
                      and polynomial(root.left) and polynomial(root.right)
                      and any(isinstance(n, ast.Name) for n in ast.walk(root.right))):
                    counts["outer_quotient_with_variable_polynomial_denominator"] += 1
                else:
                    raise ValueError(f"Unclassified terminal expression in {name}")
        kinds[name] = dict(counts)
    return {"audit_kind": "static JSON syntax inspection", "operator_sha256": OPERATOR_PIN,
            "terminal_state_sha256": STATE_PIN, "operator_entries": entries,
            "distinct_positive_operator_division_rhs_ast": len(rhs),
            "all_operator_exponents_positive_rational_literals": True,
            "all_operator_sqrt_arguments_positive_on_shell": True,
            "terminal_entry_shapes": kinds, "candidate_evaluated": False,
            "scientific_modules_imported": False, "certificate_accepted_by_audit": False}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper", type=Path, default=Path("."))
    print(json.dumps(inspect(parser.parse_args().paper), indent=2))
