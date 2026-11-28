"""Complete 20-variable H3/H4 of the selected reduced Gaunt model.

The invariant coefficients are calculated from the verified spatial
functionals, not tuned to loop outputs. N=3 with ten selected TT modes.
This module contains no higher-order gravitational vertices.
"""

from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from math import factorial, pi, prod, sqrt
from pathlib import Path
import json
import random

from reduced_gaunt_exact import five_polarizations
from reduced_gaunt_vertices import (
    tensor_from_labels, curvature_data, york_source_and_raw_gram, contract,
    integrate, coulomb, york_norm, invariants, solve_rational,
    label_matrix, label_add, label_multiply as lm, label_trace as lt,
)


def matrix_key(a):
    return tuple(tuple(F(x) for x in row) for row in a)


ZERO = matrix_key([[0]*3 for _ in range(3)])


def cross_invariants(a, b, c, d):
    return [lt(lm(a,b))*lt(lm(c,d)), lt(lm(a,c))*lt(lm(b,d)),
            lt(lm(a,d))*lt(lm(b,c)), lt(lm(lm(lm(a,b),c),d)),
            lt(lm(lm(lm(a,c),b),d))]


def paired_invariants(a, b):
    return invariants(a,b)["quartic"][2:]


def matrix_rank(rows):
    if not rows:
        return 0
    a = [[F(x) for x in row] for row in rows]
    pivot = 0
    for col in range(len(a[0])):
        target = next((i for i in range(pivot,len(a)) if a[i][col]), None)
        if target is None:
            continue
        a[pivot], a[target] = a[target], a[pivot]
        divisor = a[pivot][col]
        a[pivot] = [x/divisor for x in a[pivot]]
        for i in range(pivot+1,len(a)):
            factor = a[i][col]
            a[i] = [x-factor*y for x,y in zip(a[i],a[pivot])]
        pivot += 1
        if pivot == len(a):
            break
    return pivot


def cross_candidates(seed):
    rng = random.Random(seed)
    basis = five_polarizations()
    for _ in range(1000):
        group = []
        for _ in range(4):
            value = ZERO
            for i in rng.sample(range(5),2):
                value = label_add(value,basis[i],F(rng.choice((-2,-1,1,2)),2))
            group.append(matrix_key(value))
        yield tuple(group)


def cross_directions():
    rows, labels = [], []
    for group in cross_candidates(103):
        row = cross_invariants(*group)
        if matrix_rank(rows+[row]) > len(rows):
            rows.append(row)
            labels.append(tuple(group))
        if len(labels) == 5:
            return labels, rows
    raise ArithmeticError("Failed to select the preregistered invariant basis")


def cross_validation_directions(old, fitted_rows):
    rows = [cross_invariants(*(matrix_key(a) for a in case["labels_Eq_Bq_Es_Bs"]))
            for case in old["exact_directions"]]
    initial_rank = matrix_rank(rows)
    selected = []
    excluded = {tuple(row) for row in fitted_rows}
    for group in cross_candidates(271):
        row = cross_invariants(*group)
        if tuple(row) in excluded:
            continue
        if matrix_rank(rows+[row]) > matrix_rank(rows):
            rows.append(row)
            selected.append(group)
        if matrix_rank(rows) == 5:
            return selected, initial_rank
    raise ArithmeticError("Held-out cross validation did not span all five invariants")


@lru_cache(None)
def profile(e, b):
    return tensor_from_labels(e,b)


@lru_cache(None)
def geometry(e, b):
    return curvature_data(profile(e,b))


def mixed_coefficient(eq,bq,es,bs):
    q, s = profile(eq,bq), profile(es,bs)
    g = geometry(eq,bq)
    vector, commutator, _ = york_source_and_raw_gram(q,s,g)
    return (integrate(commutator,4)/24 + york_norm(vector)
            +6*coulomb(contract(s,s),g["r2"]))


def momentum_coefficient(es,bs):
    s = profile(es,bs)
    square = contract(s,s)
    return -7*coulomb(square,square)


def exact_dot(a,b):
    return sum((F(x)*F(y) for x,y in zip(a,b)),F(0))


def evaluate_mixed(coefficients, labels):
    eq,bq,es,bs = labels
    value = F(0)
    for qi,q in enumerate((eq,bq)):
        for si,s in enumerate((es,bs)):
            value += exact_dot(coefficients["paired"][f"{qi}{si}"],paired_invariants(q,s))
    return value + exact_dot(coefficients["cross"], cross_invariants(*labels))


def calculate_coefficients(progress=None):
    old_path = Path(__file__).resolve().parents[2]/"reports/quantum-bridge/data/reduced-vertices-results.json"
    old = json.loads(old_path.read_text())
    assert old["source_commit"] == "fb41c4bb0b423696ea12eb90cf75afd6214bc64c"
    base = [matrix_key(x) for x in five_polarizations()]
    axial, d, xz = base[1],base[0],base[3]
    pair_rows = [paired_invariants(axial,x) for x in (axial,d,xz)]
    paired, raw = {}, []
    for qi in range(2):
        for si in range(2):
            values = []
            for name, s in zip(("aligned","orthogonal","noncommuting"),(axial,d,xz)):
                labels = [ZERO]*4
                labels[qi], labels[2+si] = axial,s
                value = mixed_coefficient(*labels)
                values.append(value)
                row={"kind":"paired","species":f"{qi}{si}","direction":name,"coefficient":str(value)}
                raw.append(row)
                if progress:
                    progress(row)
            paired[f"{qi}{si}"] = solve_rational(pair_rows,values)
    cross_labels, cross_rows = cross_directions()
    cross_values=[]
    for index, labels in enumerate(cross_labels):
        total = mixed_coefficient(*labels)
        pure = evaluate_mixed({"paired":paired,"cross":[0]*5},labels)
        cross_values.append(total-pure)
        row={"kind":"cross","direction":index,"coefficient":str(total),
             "labels":[[[str(x) for x in row] for row in a] for a in labels]}
        raw.append(row)
        if progress:
            progress(row)
    cross = solve_rational(cross_rows,cross_values)
    p_labels=[(axial,ZERO),(ZERO,axial),(axial,axial),(axial,d),(axial,xz)]
    p_values=[]
    for index,labels in enumerate(p_labels):
        value=momentum_coefficient(*labels)
        p_values.append(value)
        row={"kind":"momentum","direction":index,"coefficient":str(value)}
        raw.append(row)
        if progress:
            progress(row)
    p_coefficients=solve_rational([invariants(*x)["quartic"] for x in p_labels],p_values)
    result={"cubic":[F(x) for x in old["potential_invariants"]["cubic_coefficients"]],
            "potential":[F(x) for x in old["potential_invariants"]["quartic_coefficients"]],
            "momentum":p_coefficients,"mixed":{"paired":paired,"cross":cross}}
    validations=[]
    for case in old["exact_directions"]:
        labels=tuple(matrix_key(a) for a in case["labels_Eq_Bq_Es_Bs"])
        expected=case["coefficients"]
        assert evaluate_mixed(result["mixed"],labels)==F(expected["H4_mixed"]),case["name"]
        assert exact_dot(p_coefficients,invariants(*labels[2:])["quartic"])==F(expected["H4_momentum"])
        validations.append({"name":case["name"],"mixed_and_momentum_residual":"0",
                            "source_run":old["run"]})
    held_out, initial_rank = cross_validation_directions(old,cross_rows)
    # These totals are preregistered in run-view.json; changing the source
    # cases requires updating that plan instead of silently changing its work.
    assert initial_rank == 3 and len(held_out) == 2
    for index, labels in enumerate(held_out):
        direct = mixed_coefficient(*labels)
        reconstructed = evaluate_mixed(result["mixed"],labels)
        assert reconstructed == direct,(index,reconstructed,direct)
        row={"kind":"held_out_cross","direction":index,"coefficient":str(direct),
             "reconstructed_coefficient":str(reconstructed),"residual":"0",
             "labels":[[[str(x) for x in row] for row in a] for a in labels]}
        raw.append(row)
        if progress:
            progress(row)
    validation_rows = [cross_invariants(*(matrix_key(a) for a in case["labels_Eq_Bq_Es_Bs"]))
                       for case in old["exact_directions"]]
    validation_rows += [cross_invariants(*labels) for labels in held_out]
    assert matrix_rank(validation_rows) == 5
    validations.append({"name":"held_out_cross_coverage","historical_rank":initial_rank,
                        "additional_directions":len(held_out),"final_rank":5,
                        "direct_functional_residual":"0"})
    return result,raw,validations


# Sparse polynomials use sorted tuples of variable indices, allowing repeated
# indices. Their short degree-three/four keys avoid a dense 20-index tensor.
def padd(a,b,factor=F(1)):
    out=dict(a)
    for key,value in b.items():
        out[key]=out.get(key,F(0))+factor*value
    return {k:v for k,v in out.items() if v}


def pscale(a,factor):
    return {k:v*factor for k,v in a.items() if v*factor}


def pmul(a,b):
    out={}
    for ka,va in a.items():
        for kb,vb in b.items():
            key=tuple(sorted(ka+kb))
            out[key]=out.get(key,F(0))+va*vb
    return {k:v for k,v in out.items() if v}


def psum(items):
    out={}
    for a in items:
        out=padd(out,a)
    return out


def pmatmul(a,b):
    return [[psum(pmul(a[i][k],b[k][j]) for k in range(3)) for j in range(3)] for i in range(3)]


def ptrace(a):
    return psum(a[i][i] for i in range(3))


def trace_product(*matrices):
    out=matrices[0]
    for a in matrices[1:]:
        out=pmatmul(out,a)
    return ptrace(out)


def polynomial_invariants(a,b):
    aa,bb,ab=trace_product(a,a),trace_product(b,b),trace_product(a,b)
    return [pmul(aa,aa),pmul(bb,bb),pmul(aa,bb),pmul(ab,ab),trace_product(a,a,b,b)]


def polynomial_cross(a,b,c,d):
    return [pmul(trace_product(a,b),trace_product(c,d)),
            pmul(trace_product(a,c),trace_product(b,d)),
            pmul(trace_product(a,d),trace_product(b,c)),
            trace_product(a,b,c,d),trace_product(a,c,b,d)]


def label_polynomial(parity,momentum):
    basis=five_polarizations()
    return [[{(2*(5*parity+k)+momentum,):a[i][j] for k,a in enumerate(basis) if a[i][j]}
             for j in range(3)] for i in range(3)]


def complete_polynomials(coefficients):
    eq,bq=label_polynomial(0,0),label_polynomial(1,0)
    es,bs=label_polynomial(0,1),label_polynomial(1,1)
    cubic=psum(pscale(a,c) for a,c in zip(
        (trace_product(eq,eq,eq),trace_product(eq,bq,bq)),coefficients["cubic"]))
    quartic=psum(pscale(a,c) for a,c in zip(polynomial_invariants(eq,bq),coefficients["potential"]))
    quartic=padd(quartic,psum(pscale(a,c) for a,c in zip(polynomial_invariants(es,bs),coefficients["momentum"])))
    for qi,q in enumerate((eq,bq)):
        for si,s in enumerate((es,bs)):
            quartic=padd(quartic,psum(pscale(a,c) for a,c in zip(
                polynomial_invariants(q,s)[2:],coefficients["mixed"]["paired"][f"{qi}{si}"])))
    quartic=padd(quartic,psum(pscale(a,c) for a,c in zip(
        polynomial_cross(eq,bq,es,bs),coefficients["mixed"]["cross"])))
    # Physical time reversal and spatial parity must survive coefficient
    # assembly, including the cross-parity mixed sector.
    for key in quartic:
        assert sum(i%2 for i in key)%2==0
        assert sum(i//2>=5 for i in key)%2==0
    assert all(i%2==0 for key in cubic for i in key)
    return cubic,quartic


def canonical_factors():
    basis=five_polarizations()
    gram=[float(F(189,8)*lt(lm(a,a)))*pi**1.5 for a in basis]
    gram += [g*5.5 for g in gram]
    factors=[]
    for g in gram:
        factors.extend((sqrt(2/g),1/sqrt(2*g)))
    return gram,factors,[sqrt(5.5)]*5+[sqrt(6.5)]*5


def derivative_at_zero(poly, indices, unit, factors):
    key=tuple(sorted(indices))
    return (float(poly.get(key,0))*unit*prod(factors[i] for i in key)
            *prod(factorial(n) for n in Counter(key).values()))


def hessian_block(poly, external, unit, factors):
    return [[derivative_at_zero(poly,(*external,i,j),unit,factors)
             for j in range(20)] for i in range(20)]


def json_coefficients(result):
    if isinstance(result,F):
        return str(result)
    if isinstance(result,dict):
        return {k:json_coefficients(v) for k,v in result.items()}
    if isinstance(result,list):
        return [json_coefficients(v) for v in result]
    return result
