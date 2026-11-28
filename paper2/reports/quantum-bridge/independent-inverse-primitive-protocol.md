# Independent inverse primitives: initial numerical protocol

This protocol is prospective. No value in it is a measured success. Its
purpose is to test the independent spatial machinery needed for the
[unreduced sixth-order reference](unreduced-sixth-order-reference.md), after
the [twenty-term reduced H6 pilot](full-sixth-order-results.md) obtained a
numerically converged mixed coefficient. It does not evaluate gravitational
W2/u3/u4, independently validate mixed H6, or compute any loop amplitude.

## Numerical and analytic independence

The new implementation uses exact real solid harmonics, full magnetic-mode
projection, radial finite differences, analytic exterior power/log sources,
and an explicit regular-origin patch. It imports none of the production
HEAT, Gaussian Fourier-bilinear, pairing, or tail implementations. The
synthetic solutions and integral targets are given in
[the control catalogue](spherical-green-control-design.md).

Scalar P=(-Delta)^(-1), C=P/8, and York Y=Delta_L^(-1) conventions are fixed.
The vector control actually solves three V=P f fields, derives h=div(V)
from those computed fields including their origin and exterior coefficients,
solves chi=P h, and combines W=-V-grad(chi)/4. Exact targets are held out;
they do not set a solved homogeneous coefficient or a boundary moment.

## Declared work

The committed configuration selects `independent_spherical_green` through
the unchanged run command. Stages are setup, angular controls, scalar
controls, York controls, and export. Python's standard library is sufficient.

- Angular support includes every mode through l=10. Sphere grids are
  12x24 and 16x32. Required checks include exact harmonic/norm identities,
  full normalized Gram and reconstruct/project tests, signed Fischer
  reconstruction, regular Cartesian jets, nontrivial Gaussian origin bounds,
  and the Kelvin traced-second-moment coefficient. Algebra and numerical
  quadrature have distinct recorded work limits, with 50 million counted
  operations across this suite.
- Seven scalar fixtures use R=4 and N=256,512,1024 intervals. They include
  l=0,1,2,8,10; the resonant p=4 logarithm; the bounded p=2 dipole;
  and both C/P normalizations. R=2,4,6 at fixed h=1/256 supplies a separate
  outer-radius slice. The origin patch has fixed physical radius 1/2.
- Two vector fixtures use R=4, N=128,256,512 and sphere grid 8x16. A
  separate 12x24 sphere calculation holds R=4,N=512 fixed. R=2,N=256
  and R=6,N=768 hold h=1/128 and sphere grid 8x16 fixed. The compactified
  exterior radial quadrature has independent order 12 in every invocation.
  Physical joins at r=1 and the origin-bridge edge r=1/2 are explicit,
  with both source traces preserved in the second solve.

The fixed physical origin patch is essential. A fixed number of cells would
leave a nonvanishing high-l homogeneous-coefficient bias in this radial FD
scheme. Tests therefore check its coefficient against an analytic target,
as well as values and derivatives on both sides of the patch.

The scalar value, derivative and origin-coefficient gates are respectively
1e-4 (scaled), 0.002 (scaled), and 0.003 (relative). Independent fourth-order
raw-grid PDE residuals have a 0.001 gate; they do not reuse ODE-derived jets.
The separate exact core-moment quadrature comparison has a 1e-7 gate.
Successive scalar value/origin errors must decrease by at least 2.5 unless
the finer error is already below 1e-10. This exception cannot hide a
regression from a small coarse error to a larger fine error.

Vector value and component gates are 0.0003 and 0.003; the integral gate is
0.003. Unsupported source channels are checked at 1e-10 before any declared
support restriction. The source pairing, full longitudinal norm, exterior
norm, and exterior flux are computed along distinct integration paths.
The localized auxiliary monopole is retained. A nonintegrable algebraic
monopole is a rejection condition, not a term to discard.

Structural comparisons are required at every resolution. Physical target
comparisons and latest radial differences are required at the finest
declared grid; coarse failures remain recorded diagnostics. Single-grid
angular and radius slices must meet all physical gates. Cross-slice checks
compare common physical points; R-dependent exterior norm/flux values are
compared to their exact targets at each R instead of being equated across R.

## Resources and provenance

One standard-library Python CPU process is allowed up to 1800 seconds and
a sampled peak RSS of 2 GiB. The wall limit is also capped by the standing
authorization deadline, 2026-11-28 16:29 UTC. Scalar unknown counts are
preflighted against 200,000; cumulative vector point evaluations across
all slices are limited to five million. Work counters bound their documented
operations, not every interpreter instruction or allocation. RSS is checked
at progress boundaries, not enforced by an operating-system reservation.

Completed rows, complete source diagnostics, and complete level records are
atomically checkpointed before advancing. Exceptions preserve these records
and fail the run. The result records source commit, run identity, effective
configuration, elapsed time and sampled peak RSS. There is no automatic
radius increase, resolution increase, or tolerance relaxation.

Passing would establish only the specified inverse-primitive diagnostics.
It would not certify quadrature error, supply a uniform physical large-N
bound, or close the independent mixed-functional, corner-state, or
Euclidean--Lorentzian requirements of the complete theory.
