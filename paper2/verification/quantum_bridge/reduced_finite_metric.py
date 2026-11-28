"""Independent finite-metric checks of the weak-field vertex ingredients.

This evaluates the full matrix exponential and its spatial/directional jets
with Decimal arithmetic. It does not call the polynomial curvature series.
An independent matrix inverse, Christoffel derivative and Ricci contraction
are used before centered differences extrapolate the weak-field coefficients.
"""

from decimal import Decimal as D, localcontext
from fractions import Fraction as F

from reduced_gaunt_exact import gaussian_derivative


def decimal_fraction(x):
    x = F(x)
    return D(x.numerator) / D(x.denominator)


def zeros():
    return [[D(0) for _ in range(3)] for _ in range(3)]


def madd(a, b, factor=D(1)):
    return [[a[i][j]+factor*b[i][j] for j in range(3)] for i in range(3)]


def mscale(a, factor):
    return [[factor*x for x in row] for row in a]


def mmul(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(3)), D(0))
             for j in range(3)] for i in range(3)]


def trace(a):
    return sum((a[i][i] for i in range(3)), D(0))


def norm(a):
    return max(abs(x) for row in a for x in row)


def polyval(poly, x):
    return sum((decimal_fraction(v) * (x[0]**p[0] if p[0] else D(1))
                * (x[1]**p[1] if p[1] else D(1)) * (x[2]**p[2] if p[2] else D(1))
                for p, v in poly.items()), D(0))


def field_values(tensor, x):
    envelope = (-sum((v*v for v in x), D(0))/2).exp()
    evaluate = lambda p: envelope*polyval(p, x)
    value = [[evaluate(p) for p in row] for row in tensor]
    first_poly = [[[gaussian_derivative(tensor[i][j], a) for j in range(3)]
                   for i in range(3)] for a in range(3)]
    first = [[[evaluate(p) for p in row] for row in t] for t in first_poly]
    second = [[[[evaluate(gaussian_derivative(first_poly[a][i][j], b))
                 for j in range(3)] for i in range(3)] for b in range(3)] for a in range(3)]
    return value, first, second


def exp_jet(value, first, second, precision):
    """Matrix exponential with four first/second jet directions.

    Three directions are spatial. The fourth is an independent variation S,
    whose mixed spatial derivatives give D_i(D exp_Q[S]).
    """
    eye = [[D(int(i == j)) for j in range(3)] for i in range(3)]
    pv, pd, pdd = eye, [zeros() for _ in range(4)], [[zeros() for _ in range(4)] for _ in range(4)]
    ev, ed, edd = eye, [zeros() for _ in range(4)], [[zeros() for _ in range(4)] for _ in range(4)]
    reciprocal_factorial = D(1)
    threshold = D(10) ** (-(precision-10))
    for n in range(1, 81):
        nv = mmul(pv, value)
        nd = [madd(mmul(pd[a], value), mmul(pv, first[a])) for a in range(4)]
        ndd = [[madd(madd(madd(mmul(pdd[a][b], value), mmul(pd[a], first[b])),
                          mmul(pd[b], first[a])), mmul(pv, second[a][b]))
                for b in range(4)] for a in range(4)]
        reciprocal_factorial /= n
        ev = madd(ev, nv, reciprocal_factorial)
        ed = [madd(ed[a], nd[a], reciprocal_factorial) for a in range(4)]
        edd = [[madd(edd[a][b], ndd[a][b], reciprocal_factorial) for b in range(4)] for a in range(4)]
        tail = reciprocal_factorial*max(norm(nv), *(norm(x) for x in nd),
                                         *(norm(x) for row in ndd for x in row))
        pv, pd, pdd = nv, nd, ndd
        if n >= 8 and tail < threshold:
            return ev, ed, edd, n, tail
    raise ArithmeticError("Matrix exponential jet did not reach its declared precision")


def inverse_and_determinant(a):
    work = [[*row, *[D(int(i == j)) for j in range(3)]] for i, row in enumerate(a)]
    det = D(1)
    for i in range(3):
        pivot = max(range(i, 3), key=lambda k: abs(work[k][i]))
        if pivot != i:
            work[i], work[pivot] = work[pivot], work[i]
            det = -det
        factor = work[i][i]
        assert factor
        det *= factor
        work[i] = [x/factor for x in work[i]]
        for k in range(3):
            if k != i:
                factor = work[k][i]
                work[k] = [x-factor*y for x, y in zip(work[k], work[i])]
    return [row[3:] for row in work], det


def geometry_at_amplitude(qvalues, svalues, amplitude, precision):
    q, dq, ddq = qvalues
    s, ds, _ = svalues
    value = mscale(q, amplitude)
    first = [mscale(x, amplitude) for x in dq] + [s]
    second = [[zeros() for _ in range(4)] for _ in range(4)]
    for a in range(3):
        for b in range(3):
            second[a][b] = mscale(ddq[a][b], amplitude)
        second[a][3] = second[3][a] = ds[a]
    g, dg, ddg, order, tail = exp_jet(value, first, second, precision)
    gi, det = inverse_and_determinant(g)
    dgi = [mscale(mmul(mmul(gi, dg[a]), gi), -1) for a in range(3)]
    gamma = [[[sum((gi[k][l]*(dg[i][j][l]+dg[j][i][l]-dg[l][i][j])
                   for l in range(3)), D(0))/2 for j in range(3)] for i in range(3)] for k in range(3)]
    dgamma = [[[[sum((
        dgi[a][k][l]*(dg[i][j][l]+dg[j][i][l]-dg[l][i][j])
        +gi[k][l]*(ddg[a][i][j][l]+ddg[a][j][i][l]-ddg[a][l][i][j])
        for l in range(3)), D(0))/2 for j in range(3)] for i in range(3)] for k in range(3)] for a in range(3)]
    ricci = [[sum((dgamma[k][k][i][j]-dgamma[j][k][i][k]
                  for k in range(3)), D(0))
              +sum((gamma[k][k][l]*gamma[l][i][j]-gamma[k][j][l]*gamma[l][i][k]
                    for k in range(3) for l in range(3)), D(0))
              for j in range(3)] for i in range(3)]
    curvature = sum((gi[i][j]*ricci[i][j] for i in range(3) for j in range(3)), D(0))
    gamma_gamma = sum((gi[i][j]*(gamma[k][i][l]*gamma[l][j][k]-gamma[k][i][j]*gamma[l][k][l])
                       for i in range(3) for j in range(3) for k in range(3) for l in range(3)), D(0))
    tangent = dg[3]
    a = mmul(gi, tangent)
    raw_gram = trace(mmul(a, a))
    divergence = [sum((gi[j][k]*(ddg[k][3][i][j]
        -sum((gamma[l][k][i]*tangent[l][j]+gamma[l][k][j]*tangent[i][l]
              for l in range(3)), D(0))) for j in range(3) for k in range(3)), D(0)) for i in range(3)]
    return {"R": curvature, "GG": gamma_gamma, "raw_gram": raw_gram,
            "divergence": divergence, "determinant_error": abs(det-1),
            "unimodular_connection_error": max(abs(sum((gamma[k][i][k] for k in range(3)), D(0))) for i in range(3)),
            "exponential_terms": order, "last_jet_term": tail}


def independent_point_check(q, s, expected_polynomials, point, steps, precision):
    with localcontext() as context:
        context.prec = precision
        x = [D(a) for a in point]
        hs = [D(a) for a in steps]
        assert len(hs) == 3 and hs[0] == 2*hs[1] == 4*hs[2]
        qvalues, svalues = field_values(q, x), field_values(s, x)
        amplitudes = sorted({sign*h*factor for h in hs for sign in (-1, 1) for factor in (1, 2)})
        finite = {h: geometry_at_amplitude(qvalues, svalues, h, precision) for h in amplitudes}
        f = (-sum((v*v for v in x), D(0))/2).exp()
        expected = {"GG2": polyval(expected_polynomials["l2"], x)*f**2,
                    "GG3": polyval(expected_polynomials["l3"], x)*f**3,
                    "GG4": polyval(expected_polynomials["l4"], x)*f**4,
                    "R2": polyval(expected_polynomials["r2"], x)*f**2,
                    "raw_gram_1": D(0),
                    "raw_gram_2": polyval(expected_polynomials["raw2"], x)*f**4}
        expected.update({"york_source_"+str(i): polyval(p, x)*f**2
                         for i, p in enumerate(expected_polynomials["vector"])})
        raw0 = trace(mmul(svalues[0], svalues[0]))
        estimates = []
        for h in hs:
            positive, negative = finite[h], finite[-h]
            twice_pos, twice_neg = finite[2*h], finite[-2*h]
            row = {
                "GG2": (positive["GG"]+negative["GG"])/(2*h**2),
                "GG3": (positive["GG"]-negative["GG"])/(2*h**3),
                "GG4": (twice_pos["GG"]+twice_neg["GG"]-4*(positive["GG"]+negative["GG"]))/(24*h**4),
                "R2": (positive["R"]+negative["R"])/(2*h**2),
                "raw_gram_1": (positive["raw_gram"]-negative["raw_gram"])/(2*h),
                "raw_gram_2": (positive["raw_gram"]+negative["raw_gram"]-2*raw0)/(2*h**2),
            }
            row.update({"york_source_"+str(i): (positive["divergence"][i]-negative["divergence"][i])/(2*h)
                        for i in range(3)})
            estimates.append(row)
        refined = {k: (4*estimates[2][k]-estimates[1][k])/3 for k in expected}
        errors = {k: abs(refined[k]-v)/(1+abs(v)) for k, v in expected.items()}
        return {
            "point": point, "steps": steps, "precision": precision,
            "expected_coefficients": {k: str(v) for k, v in expected.items()},
            "centered_estimates": [{k: str(v) for k, v in row.items()} for row in estimates],
            "richardson_coefficients": {k: str(v) for k, v in refined.items()},
            "scaled_errors": {k: str(v) for k, v in errors.items()},
            "maximum_scaled_error": str(max(errors.values())),
            "maximum_determinant_error": str(max(row["determinant_error"] for row in finite.values())),
            "maximum_connection_trace_error": str(max(row["unimodular_connection_error"] for row in finite.values())),
            "maximum_exponential_terms": max(row["exponential_terms"] for row in finite.values()),
            "maximum_last_jet_term": str(max(row["last_jet_term"] for row in finite.values())),
            "finite_amplitudes": len(finite),
        }
