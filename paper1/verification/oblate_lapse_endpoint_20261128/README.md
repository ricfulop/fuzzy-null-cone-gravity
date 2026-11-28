# Oblate lapse endpoint checker

**Prepared source only.** `check_lapse_endpoints.py` is a standalone
SymPy 1.14 checker for the two endpoint identities prepared in the
`ae32314d27e20ac682c5ed1b71b4d7f11bee7133` publication of
`reports/oblate-common-domain-preparation-2026-11-28`. The mathematical
definitions in that note were reviewed against source revision
`773ccc4e3e8a1a89e06086163042fbf45e171c76`.

It fixes $\epsilon=1/100$, the radial shell $[1,2]$, and the
$N=3$, $m=0$, $\ell=2$ retained scalar $Y=P_2(x)$. The checker derives
the linearized scalar curvature from varied Christoffels on the actual
oblate spatial metric, checks $\det\sigma=\Sigma^2$ for the positive
density $\Sigma$ on the selected shell, derives the two-Hessian flux
independently, and compares the Dirichlet radial fluxes with the
reviewed endpoint predictions before emitting an atomic JSON payload.

The symbolic identities are stored without the common angular
normalizations; the scalar-curvature Green identity carries the
$1/(8\pi)$ prefactor in the report convention, while the Hessian identity
carries $1/(4\pi)$. The endpoint sign samples are exact controls at
$x=0$, not integrated or global angular sign certificates.

No scientific execution has been run from this preparation. The output,
when later produced by a bounded runner, will certify only these selected
original lapse/adjoint endpoint identities and trace controls. It will
not establish a common canonical operator domain, Hamiltonian evolution
invariance, a physical null boundary, a final Dirac classification, or
a physical degree-of-freedom count.
