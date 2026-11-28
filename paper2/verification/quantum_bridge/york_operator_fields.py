"""Pointwise York operator from physical metric and contravariant-vector jets.

This implements reports/quantum-bridge/york-operator-implementation.md, using
only standard-library scalar arithmetic. It constructs neither polynomial
geometry nor Green functions, and performs no evaluation on import.

Spatial dimension is three. The metric and its supplied jets must describe a
unimodular metric (det g=1 throughout a neighborhood), such as exp(Q) with
trace-free Q. The formulas are NOT the general non-unimodular York operator.
Symmetry, unimodularity and consistency of the supplied jets/series are the
caller's responsibility; this module does not silently symmetrize them.

Index conventions, identical for finite fields and each series coefficient:
  g[i][j]          = g_ij
  dg[a][i][j]      = partial_a g_ij
  ddg[a][b][i][j]  = partial_a partial_b g_ij
  W[k]            = W^k                 (contravariant, held fixed)
  dW[a][k]        = partial_a W^k
  ddW[a][b][k]    = partial_a partial_b W^k
  inverse[i][j]   = g^ij
  gamma[k][i][j]  = Gamma^k_ij
Outputs L[i][j] and dL[a][i][j] are covariant tensor jets. A[i] is a
COVECTOR: g^ja (partial_a L_ij - Gamma^m_ai L_mj - Gamma^m_aj L_im).

Use floats or Decimal throughout an invocation; integers are allowed with
either. Decimal inputs retain the active decimal context. Mixing float and
Decimal is rejected instead of introducing implicit binary rounding. Plain
integer-only input is promoted to float. No NumPy dependency is required.

Public API:
  york_apply(g, dg, ddg, W, dW, ddW)
      -> {L, dL, A, inverse, gamma, determinant}
  longitudinal_jet(g, dg, ddg, W, dW, ddW)
      -> {L, dL}, without matrix inversion or connection construction
  york_series(metric, dmetric, ddmetric, inverse, gamma, W, dW, ddW,
              max_order=2)
      -> {L: [T0,...], dL: [dT0,...], A: [A0W,...], max_order}

Series arguments have an OUTER homogeneous-order index; inverse and gamma
are supplied by the caller. gamma[0]=None is accepted as a zero connection.
Coefficients use g(lambda)=sum lambda^n*g_n, not derivative normalization:
A[2] is A2W, whereas the second lambda derivative is 2*A2W. W and all its
spatial jets remain fixed across orders. Feed physical Gaussian-polynomial
metric jets and physical HeatField vector jets, including all envelopes;
this module does not infer or apply a shared Gaussian envelope.
"""

from decimal import Decimal
from math import isfinite


def _has_decimal(value):
    if isinstance(value, (list, tuple)):
        return any(_has_decimal(item) for item in value)
    return isinstance(value, Decimal)


def _array(value, shape, name, decimal_mode):
    """Copy a small nested array, checking shape and preserving one scalar type."""
    if shape:
        if not isinstance(value, (list, tuple)) or len(value) != shape[0]:
            raise ValueError(name + " must have shape " + str(shape))
        return [_array(item, shape[1:], name + "[" + str(i) + "]", decimal_mode)
                for i, item in enumerate(value)]
    if isinstance(value, bool) or not isinstance(value, (int, float, Decimal)):
        raise TypeError(name + " must be a real int, float, or Decimal")
    if decimal_mode:
        if isinstance(value, float):
            raise TypeError("do not mix float and Decimal jets; convert explicitly before calling")
        result = value if isinstance(value, Decimal) else Decimal(value)
        if not result.is_finite():
            raise ValueError(name + " must be finite")
    else:
        result = float(value)
        if not isfinite(result):
            raise ValueError(name + " must be finite")
    return result


def _finite_inputs(g, dg, ddg, W, dW, ddW):
    values = (g, dg, ddg, W, dW, ddW)
    shapes = ((3, 3), (3, 3, 3), (3, 3, 3, 3), (3,), (3, 3), (3, 3, 3))
    names = ("g", "dg", "ddg", "W", "dW", "ddW")
    decimal_mode = any(_has_decimal(value) for value in values)
    return tuple(_array(value, shape, name, decimal_mode)
                 for value, shape, name in zip(values, shapes, names))


def _inverse(matrix):
    """Three-by-three Gauss-Jordan inverse with partial pivoting and determinant."""
    zero = matrix[0][0] * 0
    one = zero + 1
    rows = [list(row) + [one if i == j else zero for j in range(3)]
            for i, row in enumerate(matrix)]
    determinant = one
    for column in range(3):
        pivot = max(range(column, 3), key=lambda i: abs(rows[i][column]))
        if rows[pivot][column] == zero:
            raise ValueError("singular metric in york_apply")
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            determinant = -determinant
        divisor = rows[column][column]
        determinant *= divisor
        rows[column] = [value / divisor for value in rows[column]]
        for i in range(3):
            if i != column:
                factor = rows[i][column]
                rows[i] = [value - factor * other
                           for value, other in zip(rows[i], rows[column])]
    return [row[3:] for row in rows], determinant


def _christoffel(inverse, dg):
    zero = inverse[0][0] * 0
    return [[[sum((inverse[k][ell] * (dg[i][ell][j] + dg[j][ell][i] - dg[ell][i][j])
                   for ell in range(3)), zero) / 2
              for j in range(3)] for i in range(3)] for k in range(3)]


def _longitudinal(g, dg, ddg, W, dW, ddW):
    zero = g[0][0] * 0
    divergence = sum((dW[k][k] for k in range(3)), zero)
    divergence_gradient = [sum((ddW[a][k][k] for k in range(3)), zero) for a in range(3)]
    L = [[sum((W[k] * dg[k][i][j] + g[i][k] * dW[j][k] + g[j][k] * dW[i][k]
               for k in range(3)), zero) - (2 * g[i][j] * divergence) / 3
          for j in range(3)] for i in range(3)]
    dL = [[[sum((dW[a][k] * dg[k][i][j] + W[k] * ddg[a][k][i][j]
                 + dg[a][i][k] * dW[j][k] + g[i][k] * ddW[a][j][k]
                 + dg[a][j][k] * dW[i][k] + g[j][k] * ddW[a][i][k]
                 for k in range(3)), zero)
            - (2 * (dg[a][i][j] * divergence + g[i][j] * divergence_gradient[a])) / 3
            for j in range(3)] for i in range(3)] for a in range(3)]
    return L, dL


def _connection_action(gamma, tensor):
    """Both lower-index connection terms, with output indices [a][i][j]."""
    zero = tensor[0][0] * 0
    return [[[sum((gamma[m][a][i] * tensor[m][j] + gamma[m][a][j] * tensor[i][m]
                   for m in range(3)), zero)
              for j in range(3)] for i in range(3)] for a in range(3)]


def longitudinal_jet(g, dg, ddg, W, dW, ddW):
    """Return the covariant tensor L_g W and its first spatial derivative.

    Also valid for a single homogeneous metric coefficient g_n, since this
    expression is linear in g and its spatial jets. No inverse is computed.
    """
    fields = _finite_inputs(g, dg, ddg, W, dW, ddW)
    L, dL = _longitudinal(*fields)
    return {"L": L, "dL": dL}


def york_apply(g, dg, ddg, W, dW, ddW):
    """Finite-metric York application with an independently built inverse/Gamma.

    W is contravariant. A is a covector. The caller must supply a consistent
    unimodular metric jet. The returned determinant is metadata, not an
    assertion of unimodularity or a replacement for checking its derivatives.
    """
    g, dg, ddg, W, dW, ddW = _finite_inputs(g, dg, ddg, W, dW, ddW)
    inverse, determinant = _inverse(g)
    gamma = _christoffel(inverse, dg)
    L, dL = _longitudinal(g, dg, ddg, W, dW, ddW)
    connection = _connection_action(gamma, L)
    zero = g[0][0] * 0
    A = [sum((inverse[j][a] * (dL[a][i][j] - connection[a][i][j])
              for j in range(3) for a in range(3)), zero) for i in range(3)]
    return {"L": L, "dL": dL, "A": A, "inverse": inverse,
            "gamma": gamma, "determinant": determinant}


def york_series(metric, dmetric, ddmetric, inverse, gamma, W, dW, ddW, max_order=2):
    """Return T_n, partial T_n and A_n W through max_order (0, 1 or 2).

    Supply coefficient arrays with at least max_order+1 entries. Spatial
    indices follow the module conventions after the outer order index.
    Extra coefficients are ignored. gamma[0]=None means a zero connection;
    an explicit order-zero array is also accepted. No metric inversion or
    Taylor reconstruction is performed. The supplied inverse/Gamma series
    must agree with the metric series. In the intended flat expansion,
    metric[0]=inverse[0]=I and its spatial derivatives and gamma[0] vanish.
    The convolution also includes a supplied nonzero gamma[0], if present.
    """
    if type(max_order) is not int or max_order not in (0, 1, 2):
        raise ValueError("max_order must be 0, 1, or 2")
    count = max_order + 1
    names = ("metric", "dmetric", "ddmetric", "inverse", "gamma")
    series = (metric, dmetric, ddmetric, inverse, gamma)
    for name, value in zip(names, series):
        if not isinstance(value, (list, tuple)) or len(value) < count:
            raise ValueError(name + " must provide coefficients 0 through " + str(max_order))
    # Look only at the requested orders; unrelated higher coefficients do
    # not influence numeric type selection or cause additional work.
    decimal_mode = any(_has_decimal(value[:count]) for value in series)
    decimal_mode = decimal_mode or any(_has_decimal(value) for value in (W, dW, ddW))
    shapes = ((3, 3), (3, 3, 3), (3, 3, 3, 3), (3, 3))
    prepared = [[_array(value[n], shape, name + "[" + str(n) + "]", decimal_mode)
                 for n in range(count)]
                for name, value, shape in zip(names[:4], series[:4], shapes)]
    metric, dmetric, ddmetric, inverse = prepared
    W = _array(W, (3,), "W", decimal_mode)
    dW = _array(dW, (3, 3), "dW", decimal_mode)
    ddW = _array(ddW, (3, 3, 3), "ddW", decimal_mode)
    zero = metric[0][0][0] * 0
    connections = []
    for n in range(count):
        if n == 0 and gamma[n] is None:
            connections.append([[[zero for _ in range(3)] for _ in range(3)] for _ in range(3)])
        else:
            connections.append(_array(gamma[n], (3, 3, 3), "gamma[" + str(n) + "]", decimal_mode))

    L, dL = [], []
    for n in range(count):
        tensor, derivative = _longitudinal(metric[n], dmetric[n], ddmetric[n], W, dW, ddW)
        L.append(tensor)
        dL.append(derivative)
    covariant, A = [], []
    for n in range(count):
        correction = [_connection_action(connections[b], L[n - b]) for b in range(n + 1)]
        covariant.append([[[dL[n][a][i][j] - sum((term[a][i][j] for term in correction), zero)
                            for j in range(3)] for i in range(3)] for a in range(3)])
        A.append([sum((inverse[b][j][a] * covariant[n - b][a][i][j]
                       for b in range(n + 1) for j in range(3) for a in range(3)), zero)
                  for i in range(3)])
    return {"L": L, "dL": dL, "A": A, "max_order": max_order}
