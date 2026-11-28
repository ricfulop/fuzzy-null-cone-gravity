# Independent terminal germ certificate — preparation

This checker is prepared but has not been executed. It requires the completed
original stationary operator export and a completed terminal DAE closure. The
active closure calculation has not yet supplied the latter input. No solution
germ or nonlinear compatibility result is asserted by this directory.

The scientific input is the original Paper I angular Galerkin coefficient
embedding, with the stationary, axisymmetric, reflection-even eight-profile
sector at fixed `epsilon = 1/100`. The seven gauge-fixed profiles omit `X0`, but
the equation associated with `X0` remains among all 22 retained Euler rows.
The separate 23rd row is `D = 6 E_X4 - E_Cpolar4`, not a retained equation.

## Terminal identities

Let the supplied fourteen-row state embedding have upper and lower blocks
`P,S`, so that the intended physical profile and its radial derivative are
`v = P y` and `v' = S y`. Let `y' = K y + B u` be the proposed reduced system.
All coefficients belong to `Q(r,T)`, with
`T = atan(1/(10r))` and `D_r = partial_r - 10/(100r^2+1) partial_T`.

The independent checker requires

```
D_r(P) + P K = S,
P B = 0,
E0 P + E1 S + E2 (D_r(S) + S K) = 0,
E2 S B = 0.
```

The last two equalities use every original retained row, including the `X0`
equation. The response is recomputed with the separate discarded row and must
equal the saved state/control response. The checker also reconstructs the
constant-control augmented ODE and verifies the same identities by directly
differentiating its physical profile map twice.

Only the original operator and terminal matrices enter this verification. The
checker does not import the discovery reducer or the curvature exporter, and
does not rely on their intermediate elimination history.

## From formal coefficients to an actual germ

The exact rational alternating series encloses `atan(1/(10r0))` at a positive
rational radius. Rational interval Horner evaluation must separate every
relevant coefficient denominator from zero. Analyticity then supplies a
nonempty open radial interval around that radius. Constant remaining controls
can be adjoined as coordinates with zero derivative, giving a linear analytic
ODE and an analytic local solution for each initial vector.

A responding witness requires both denominator regularity and separation of a
response numerator at the **same radius**. A coordinate unit vector supplies
the corresponding initial state or constant control. Failure to separate a
formally nonzero response is inconclusive. An identically zero response proves
silence only on the certified forward subfamily; a universal zero would need a
converse parametrization or an explicit differential-consequence certificate.

## Meaning of a responding witness

For a real stationary even solution `e`, the known retained odd coordinate
perturbation `g_f` has zero unrestricted linearized Einstein residual. The
prepared coefficient-map identity gives

```
B_xi[e + g_f] = (4/105) integral (r^2 + epsilon) w f D[e] dr.
```

If `D[e]` is nonzero at the certified radius, continuity permits compact real
profiles `f,w` supported in a constant-sign neighborhood with positive overlap.
The integral is then nonzero. The second-order Noether necessary condition
excludes a family of exact retained bulk solutions on a common spacetime
neighborhood, twice differentiable in amplitude in the topology controlling
the Euler jets, with this specified spacetime tangent `e + g_f`.

That implication allows arbitrary time-dependent second-order corrections. It
does not assert that the oblate background is inadmissible, that every nearby
solution is obstructed, that every linearized gauge representative behaves the
same way, or that a canonical Dirac classification, physical degree-of-freedom
count, or global radial boundary realization has been obtained.

## Launch preparation

A future tracked, bounded local invocation must pin SHA-256 hashes of the
original `stationary-operators.json` and completed `dae-closure.json`, freeze the
reviewed checker source, preserve stdout/stderr and process outcome, and write
to a fresh output location. Its entry point is:

```
python certify_germ.py --config CONFIG --operators ORIGINAL_OPERATORS \
  --closure TERMINAL_CLOSURE --output NEW_CERTIFICATE
```

No input hashes or completed result are invented while the closure input is
absent. Static preparation does not execute the manufactured scientific
controls; those belong to the eventual bounded invocation.
