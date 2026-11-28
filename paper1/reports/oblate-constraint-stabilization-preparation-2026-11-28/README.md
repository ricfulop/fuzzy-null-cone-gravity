# Oblate constraint stabilization: N=3/4 preparation atlas

**Preparation only.** This note inventories the certified oblate constraint
facts for the original Paper I coefficient-Galerkin action and derives the
next preservation equations in notation. It adds no experiment, no executable
symbolic checker, no manuscript change, no rank upgrade, and no physical
degree-of-freedom count.

The target background is the exact static oblate Minkowski family

$$
\sigma_\epsilon
=\frac{\Sigma}{\rho}dr^2+\frac{\Sigma}{1-x^2}dx^2
  +\rho(1-x^2)d\phi^2,\qquad
\rho=r^2+\epsilon,\qquad
\Sigma=r^2+\epsilon x^2,
$$

with $a=h_{tt}=-1$, $b=h_{tr}=0$, $W=A_t=0$, $p=0$ and
$\Lambda=0$. All coefficient fields of this family have angular degree at
most two, so it belongs to the actual finite band for $N=3$ and $N=4$.

## 1. Certified data at the exact oblate solution

Let $L=N-1$. The temporal coefficients comprise one lapse coefficient $a$,
one radial-shift coefficient $b$, and one angular-shift vector $W$:

$$
\#a=(L+1)^2,\qquad
\#(b,W)=3(L+1)^2-2,\qquad
\#z=4(L+1)^2-2.
$$

The exact small-oblate all-mode calculation proves that the only nonzero
temporal Hessian block at the static point is the shift block

$$
A(s,t)=-\frac12G(n_s,n_t),\qquad
n_s=(1-\Pi)\mathcal L_{v_s}\sigma_\epsilon,
\qquad v_s=(b/\lambda)\partial_r+W^A\partial_A ,
$$

and that radial derivatives in $\partial_r b$ and $\partial_rW$ cancel
because they are exact retained spatial tangents. The DeWitt form is
indefinite; the leading minor proof checks the Gram rank on the discarded
image, not a positive minimization.

Sources: the [all-mode rank report](https://github.com/ricfulop/fuzzy-null-cone-gravity/blob/773ccc4e3e8a1a89e06086163042fbf45e171c76/paper1/reports/oblate-constraint-rank-2026-11-27/README.md),
its [all-m verifier](https://github.com/ricfulop/fuzzy-null-cone-gravity/blob/773ccc4e3e8a1a89e06086163042fbf45e171c76/paper1/verification/oblate_constraint_rank_20261127/all_modes.py),
and the integrated [canonical oblate manuscript section](https://github.com/ricfulop/fuzzy-null-cone-gravity/blob/773ccc4e3e8a1a89e06086163042fbf45e171c76/paper1/sections/canonical_oblate.tex).

| item | $N=3$ | $N=4$ | established status |
|---|---:|---:|---|
| lapse coefficient functions | 9 | 16 | exact coefficient count |
| shift coefficient functions | 25 | 46 | exact coefficient count |
| temporal primaries | 34 | 62 | exact coefficient count |
| shift-Hessian rank $R=8L+2$ | 18 | 26 | exact for every sufficiently small nonzero $\epsilon$ on a compact positive-radius shell |
| shift-kernel dimension $3L^2-2L-1$ | 7 | 20 | exact, equal to the tangent-shift kernel |
| primary-plus-secondary Poisson rank $2R$ | 36 | 52 | exact at the primary/secondary stage only |
| temporal-Hessian nullity | 16 | 36 | exact at this static $p=0$ point; not a gauge count |

The exact kernel is the direct sum of

$$
b/\lambda\in V_{L-2},\qquad W^A\partial_A x\in V_{L-1},
$$

with arbitrary radial profiles. The radial part has dimension
$(L-1)^2$, while the angular part has dimension $2L^2-2$. For $N=3$,
this gives one base-radial shift plus six angular shifts; the angular
six split into the three rotations $K_i$ and the three background-null
directions $n_iK_z$. These product fields are still retained: the
$K_z$ multiples have $W^A\partial_Ax=0$, and a degree-one scalar
$n_i$ times the degree-one vector harmonic $K_z$ decomposes only into
$j\le2$ vector harmonics by rotational equivariance. For $N=4$, the
same formula gives four radial kernel directions and sixteen angular
kernel directions. The rotations are exact ansatz-preserving symmetries.
The remaining background-null directions are not thereby classified as
nonlinear gauge generators.

The mixed lapse/shift secondary bracket factors through the same normal
map:

$$
\{C_a[f],C_s[s]\}
=-\frac12G(\nabla^2 f,n_s).
$$

It annihilates $\ker A$. After eliminating the $18$ or $26$
nondegenerate shift primary-secondary pairs, the certified
primary-plus-secondary rank is therefore $36$ or $52$. No tertiary or
higher-preservation rank is included in that number.

## 2. Linearized lapse prolongation

On the exact static family, the corrected lapse checker eliminates the
nondegenerate shift complement $R$ and defines

$$
K_{\rm eff}=K-J_RA_R^{-1}J_R^*,\qquad
L_fh=\frac12G(\nabla^2 f,Th),\qquad
\mathcal T_f(p)=L_fK_{\rm eff}p.
$$

Here $K=2(T^*GT)^{-1}$ is the spatial kinetic inverse and $J_R$ maps the
nondegenerate shift complement into retained spatial velocity. All
adjoints retain the radial integrations by parts.

Let $P_{\mathcal N}$ project, in the integrated DeWitt form, onto the
nondegenerate normal shift image. The reduced lapse-tertiary bracket is

$$
\{L_fh,\mathcal T_g\}_{\rm lin}
=-\frac12G(Z_f,Z_g),\qquad
Z_f=(1-P_{\mathcal N})(1-\Pi)\nabla^2f.
$$

The checker proves exact gradient identities for the lapse descriptors
whose spatial gradients are retained shifts. These make the associated
tertiaries redundant modulo the remaining shift secondaries:

| item | $N=3$ | $N=4$ | established status |
|---|---:|---:|---|
| redundant lapse-profile families | 1 | 4 | exact retained-gradient identities |
| complementary real lapse channels | 8 | 12 | exact leading forms plus a small-deformation $H^1_0$ weak-form argument |

This is still a **linearized** prolongation result. It proves additional
independent lapse-preservation equations near the oblate solution. It
does not prove nonlinear closure of the remaining null directions, a
constant-rank nonlinear constraint manifold, a common original canonical
domain, or a final Dirac count.

Sources: the [corrected lapse report](https://github.com/ricfulop/fuzzy-null-cone-gravity/blob/773ccc4e3e8a1a89e06086163042fbf45e171c76/paper1/reports/oblate-lapse-consistency-2026-11-27/README.md)
and its [lapse checker](https://github.com/ricfulop/fuzzy-null-cone-gravity/blob/773ccc4e3e8a1a89e06086163042fbf45e171c76/paper1/verification/oblate_constraint_rank_20261127/lapse_prolongation.py).

## 3. Symbolic next preservation equations

Write the canonical variables as

$$
y=(q,p,z),\qquad q=(\lambda,X,V,C),\qquad
z=(a,b,W),
$$

and let $\pi_A=0$ be the primaries conjugate to the temporal
coefficients $z^A$. The coefficient Hamiltonian $H(q,p,z)$ gives
secondaries

$$
C_A=D_{z^A}H .
$$

With total Hamiltonian $H_T=H+\mu^A\pi_A$, their first preservation is

$$
0=\dot C_A=\{C_A,H\}_{q,p}
  +D_{z^B}C_A\,\mu^B
  \equiv F_A+A_{AB}\mu^B .
$$

On the exact oblate solution, split temporal tests into a
nondegenerate shift complement $R$ and null tests $i$. In the
compact-test linearization at this fixed static background, with this
constant-rank complement on the shell, the secondary row first gives

$$
C_{R,1}=J_R^*p+A_{RR}\delta s_R=0,
$$

with zero null cross block, and algebraically solves the corresponding
shift-coordinate perturbation $\delta s_R$. The full block row of
first-preservation equations is formally

$$
A_{RR}\mu^R=-F_R-A_{Ri}\mu^i .
$$

At $\mu_0=0$, its linearization is a row in $\delta F_R$ and
$\delta\mu$. Since $A_{RR}$ is invertible on the complement at the
background, it determines the corresponding multiplier rates:

$$
\mu^R=-A_{RR}^{-1}(F_R+A_{Ri}\mu^i).
$$

Substitution leaves reduced linearized first-preservation equations

$$
\Phi_i(y,\mu^j)=0 .
$$

This is not an unconditional nonlinear implicit-function step. At a
generic point, $D_zC$ can be a radial differential operator; solving the
nonlinear secondaries or their preservation equations in a neighborhood
requires an isomorphism on a specified common function-space domain, or
a separate proof that the needed algebraic block structure persists.

The lapse components of the linearization of $\Phi_i$ are the
$\mathcal T_f(p)=0$ equations above. The null-shift components vanish
at linear order by exact tangency and flat spatial diffeomorphism
invariance, which is why the next nontrivial test for them is quadratic.

The next stabilization equations have the schematic form

$$
0=\dot\Phi_i
=D_q\Phi_i[\dot q]
 +D_p\Phi_i[\dot p]
 +D_z\Phi_i[\mu]
 +D_{\mu^j}\Phi_i[\nu^j],
\qquad
\nu^j=\dot\mu^j .
$$

For a complementary lapse channel this linearizes to the
lapse-tertiary chain

$$
\pi_a\longrightarrow Lh\longrightarrow\mathcal T(p)
       \longrightarrow\mathcal Q(h,\delta a),
$$

where the bracket with $\delta a$ is nondegenerate only on the eight or
twelve complementary real channels. The low retained-gradient lapse
tertiaries reduce to null shift secondaries.

For a non-rotation null angular shift descriptor $k$, the first
possible nonzero term is quadratic. The rotation descriptors are exact
and remain zero. Holding the background descriptor fixed,

$$
C_k=D_zH[k],\qquad
\mathcal F_k=\{C_k,H\}_{q,p}+D_zC_k[\mu].
$$

At the oblate static point, $H_1=0$, $\mu_0=0$, and
$C_{k,1}=J_k^*p$. Together with $V''J_k=L_fJ_k=0$, these identities
make both $\mathcal F_k$ and $D\mathcal F_k$ vanish strongly on the
compact-test phase-space tangent. If $H_j$ and $C_{k,j}$ denote
homogeneous Taylor terms, including factorials, the quadratic
coefficient is

$$
\mathcal O_k[y_1,\mu_1]
=\{C_{k,1},H_3\}_{q,p}
 +\{C_{k,2},H_2\}_{q,p}
 +D_zC_{k,2}[\mu_1].
$$

The coefficient-space Noether identity supplies the matching necessary
condition on admissible first-order histories without building the
cubic Hamiltonian. For the angular vector $\xi$ corresponding to $k$,

$$
B_\xi[c_1]
=E_{\rm full}^{(1)}[c_1]\bigl[Q_cDF_\xi(c_0)c_1\bigr]=0 .
$$

The equality $\mathcal O_k=-B_\xi$ is an on-history comparison: the
first-order fields must solve the retained linearized equations and
supply their consistent time and radial derivatives. The projector
$Q_c$ is the fixed coefficient-band projector, not the DeWitt normal
projector. The unknown second amplitude correction drops out because
$E_{\rm full}(c_0)=0$ and $Q_cF_\xi(c_0)=0$ at the exact oblate
background.

Sources: the [nonlinear null preparation](https://github.com/ricfulop/fuzzy-null-cone-gravity/blob/773ccc4e3e8a1a89e06086163042fbf45e171c76/paper1/reports/oblate-nonlinear-null-preparation-2026-11-27/README.md),
the [Noether compatibility source](https://github.com/ricfulop/fuzzy-null-cone-gravity/blob/773ccc4e3e8a1a89e06086163042fbf45e171c76/paper1/verification/oblate_nonlinear_null_20261127/null_compatibility.py),
and the [original coefficient action section](https://github.com/ricfulop/fuzzy-null-cone-gravity/blob/773ccc4e3e8a1a89e06086163042fbf45e171c76/paper1/sections/nonlinear_action.tex).

## 4. What c925 resolves

The accepted direct-state certificate supplies a real interval-local
stationary polar germ $e$ for the original $N=3$, $m=0$,
reflection-even, $X_0=0$ problem at $\epsilon=1/100$. It verifies the
fourteen-state identities

$$
CM=0,\qquad JM=I_8,\qquad (DC+CK)M=0,
$$

the seven-profile lift, all 22 retained Euler rows including the
restored $X_0$ row, 57 same-radius proof denominators, a nonzero
rank-six minor of the original $6\times14$ constraint matrix, and a
nonzero discarded response numerator at $r_0=1$. The certified
responding initial coordinate is the second standard vector in the
eight-dimensional forward family.

This is a direct Lagrangian forward-germ certificate at a fixed oblate
deformation. It does not require identifying $\epsilon=1/100$ with the
unspecified small-deformation interval used in the canonical rank
theorem.

The coefficient-map check and stationary-seed calculation then add the
known odd continuum coordinate leg

$$
g_f:\quad C_{{\rm a},2}=-\frac{2\rho f}{3},\qquad
V_{{\rm a},2}=-\frac{f'}3,\qquad W=0,
$$

and the axisymmetric background-null descriptor

$$
\xi=w(r)x\partial_\phi .
$$

For the mixed tangent $c_1=e+g_f$,

$$
B_\xi[c_1]
=\frac4{105}\int dr\,\rho\,w f
\left(6\mathscr E_{X,4}^{(1)}[e]
      -\mathscr E_{C_{\rm p},4}^{(1)}[e]\right),
$$

with no extra radial measure factor. The certificate gives a radius at
which the parenthesized response is nonzero; shrinking to a compact
shell and choosing $f,w$ with $fw\ge0$ and positive overlap makes the
integral nonzero.

This proves that the retained linearized Euler equations have at least
one spurious tangent direction: the specified spacetime tangent
$e+g_f$ cannot be tangent to any exact retained bulk solution family
that is $C^2$ in amplitude on a common spacetime neighborhood and in a
topology controlling the Euler jets. It also proves that the
second-order correction, including arbitrary time-dependent temporal
data, cannot cancel this selected Noether test.

This does **not** by itself prove that the intrinsic exact-solution set
is a singular variety. A nonregular defining equation such as
$x^2=0$ can have a smooth zero set with an enlarged linearization. The
c925 result proves failure of linearization stability for a specified
tangent of the retained Euler presentation. It does not determine the
regularized tangent space of the exact-solution set, remove every
accidental null direction, establish an $N=4$ analogue, or add a
tertiary rank to the primary-plus-secondary ranks above.

Sources: the accepted [direct stationary-germ certificate report](https://github.com/ricfulop/fuzzy-null-cone-gravity/blob/773ccc4e3e8a1a89e06086163042fbf45e171c76/paper1/reports/oblate-stationary-state-germ-pass-2026-11-28/README.md),
the [input-domain audit](https://github.com/ricfulop/fuzzy-null-cone-gravity/blob/773ccc4e3e8a1a89e06086163042fbf45e171c76/paper1/reports/oblate-stationary-input-domain-audit-2026-11-28/README.md),
the [stationary-to-nonlinear linkage note](https://github.com/ricfulop/fuzzy-null-cone-gravity/blob/773ccc4e3e8a1a89e06086163042fbf45e171c76/paper1/reports/oblate-stationary-nonlinear-linkage-preparation-2026-11-28/README.md),
and the integrated [local obstruction proposition](https://github.com/ricfulop/fuzzy-null-cone-gravity/blob/773ccc4e3e8a1a89e06086163042fbf45e171c76/paper1/sections/oblate_local_obstruction.tex).

## 5. Smallest calculation that would constrain the nonlinear stratum

The narrowest next calculation is an $N=3$ quadratic null-shift atlas at
the same exact oblate background.

**Operator to derive.** For the three non-rotation angular kernel
descriptors

$$
\xi_i=w_i(r)n_iK_z,\qquad i=x,y,z,
$$

derive the discarded coefficient maps

$$
Q_cDF_{\xi_i}(c_0)c_1
$$

in all magnetic sectors. Rotationally equivariant Lie products make the
discarded coefficient maps land in degrees at most four for $N=3$; this
is not a degree bound on the full continuum Einstein residual. Pair the
resulting discarded coefficient tests with the full linearized
continuum coefficient Euler residual and export the radial density rows

$$
B_i[c_1]
=E_{\rm full}^{(1)}[c_1]
  \bigl[Q_cDF_{\xi_i}(c_0)c_1\bigr].
$$

The $i=z$ row must reduce to the accepted $4/105$ response only on the
special mixed tangent $e+g_f$, where the odd leg is the pure coordinate
mode and $E_{\rm full}^{(1)}[g_f]=0$. A generic axisymmetric odd leg
retains the other cross-pairing and the other discarded tests.

**How it constrains the stratum.** Any first-order spacetime tangent to
a $C^2$ exact retained family must solve

$$
B_i[c_1]=0
$$

as a compact-test distribution for all three $i$. The equations are
quadratic conditions on the retained linearized solution space. They
can cut spurious tangent directions, as the c925 example shows, but
they should not be converted into additive Dirac ranks until their
independence, their preservation, and their common radial domain are
proved.

**Positive controls.**

1. Recover the archived axisymmetric formula
   $X_4=-12wC_{{\rm a},2}/35$ and
   $C_{{\rm p},4}=2wC_{{\rm a},2}/35$, and its resulting
   $4/105$ compact-test normalization for the odd coordinate leg.
2. Pair the atlas with the c925 responding stationary germ and reproduce
   a nonzero $i=z$ response on the same certified local shell.
3. Keep the affine $\partial_r\xi_i$ contribution in $F_\xi V$ and
   verify it is retained at the background; no radial-freezing shortcut
   is allowed.

**Negative and zero controls.**

1. The three exact rotations $K_i$ must have identically vanishing
   discarded maps at the background.
2. Pure continuum-coordinate first-order seeds must give zero
   unrestricted linearized Einstein residual even if their finite
   coordinate orbits leave the coefficient band.
3. Pure reflection parity in the axisymmetric $n_zK_z$ channel must
   give zero by symmetry; a zero pure-parity test is not evidence of a
   nonlinear gauge symmetry.
4. An off-shell conformal perturbation should remain a nonzero density
   control, confirming that the Euler pairing has not silently
   discarded the relevant continuum residual.

**Completion criteria.** A completed checkpoint should list the exact
coefficient rows for $\xi_x,\xi_y,\xi_z$, verify harmonic
projection/reconstruction and the rotation zeros, state which rows are
exact identities on certified first-order solution families, and give
at least one same-domain nonzero response if one is found. A zero result
must be reported as a hidden identity only after proving that the tested
first-order family covers the claimed sector; a forward-family zero
alone would be weaker.

## 6. Relation to the remaining Paper I gate

This atlas would address the first nonlinear null directions after the
primary-plus-secondary oblate rank theorem and the linearized lapse
tertiaries. It would still leave the full Paper I gate open:

- the complementary lapse chains need nonlinear preservation on a
  common domain, not only linearized nondegenerate brackets;
- the null-shift quadratic conditions need their own preservation;
- the $N=4$ angular kernel contains thirteen non-rotation angular
  directions not tested by the $N=3$ axisymmetric obstruction;
- radial boundary conditions must make the canonical generators and
  their adjoints act on the same physical domain;
- only after these pieces are regular and independent can one classify
  first- and second-class constraints and count physical degrees of
  freedom.

The right next output is therefore another preparation or verifier for
the quadratic null-shift atlas, not a manuscript rank update.

## Preparation review and source access

Independent read-only review found no remaining mathematical or scope blocker (msg_0457afae). Source links pin the reviewed repository revision; the local [validation record](validation.json) records source hashes, review scope and packaging. No proposed scientific calculation has been executed, and neither report establishes a common physical domain or a final constraint classification.
