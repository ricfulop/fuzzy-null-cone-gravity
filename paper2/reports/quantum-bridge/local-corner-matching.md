# Local area matching for the reduced canonical Gaunt candidate

## Scope and result

This report uses **the same reduction-first candidate** as
[reduced-canonical-gaunt.md](reduced-canonical-gaunt.md): ten canonical
TT mode pairs, the ordinary matrix exponential, complete Gaunt products,
the dual-Gram momentum map, and the untruncated dependent conformal
solution. It also uses the area--boost polarization of
[corner-quantization-and-transfer.md](corner-quantization-and-transfer.md).
No new microscopic state multiplicity or change to the bulk Hamiltonian
is assumed. The corner remains a prescribed coordinate/observer sphere
in the asymptotically flat construction, not yet a relationally defined
causal-diamond boundary.

The bulk canonical form was obtained by integrating the ADM potential
over the whole asymptotically flat slice, not just its portion inside the
sphere. Adding the independent corner pair below is therefore a test of
the proposed extension, not an already derived subregion symplectic form.

**New derivations below:** the local area expansion, its matching-constraint
algebra, an explicit ten-mode nonzero bracket, the dimension and Dirac
bracket of the constrained system, and a local canonical conversion.
These are analytic calculations, not results of an additional numerical
run. The only external input used for the broader corner interpretation
is the continuum surface-symmetry construction identified in Section 7.

The main conclusions are:

- The original constraints `independent area = reconstructed bulk area`
  are all first class **only if the reconstructed area smearings mutually
  Poisson commute**. They do not do so in the finite-radius ten-mode model.
- At a regular point their nonzero bracket rank gives second-class
  directions; the kernel gives the remaining first-class directions.
  Fixing every boost to zero is therefore not, in general, a gauge fixing
  of those original constraints.
- Imposing both matching and zero boost as second-class conditions has an
  exact, canonical, twenty-dimensional reduction. This is a specified
  elimination of all independent edge variables, not a derivation that
  they were all gauge.
- There is also an explicit local first-class conversion using the existing
  boost pairs, with valid gauge `eta = 0`. It changes the area matching
  off that section. Its identification with gravitational edge charges
  requires boundary information absent from the candidate.

## 1. Define the area density before truncating its harmonics

Use the candidate's dimensionless packet width one and $16\pi G=1$ for
all coefficient calculations below. In these units $1/(8\pi G)=2$.
The bulk brackets are

$$
\Omega_B=\sum_{\alpha=1}^{10}dp_\alpha\wedge dq^\alpha,
\qquad \{q^\alpha,p_\beta\}_B=\delta^\alpha_\beta.
$$

Let $n=\mathbf x/R$, and let $d\Omega$ be the unit-sphere measure. Define
$F(n)=a_{\rm bulk}(n)$ by $F(n)d\Omega=dA/(8\pi G)$, so that
$F[1]=\int Fd\Omega=A_R/(8\pi G)$. In arbitrary units the area formula
is $R^2\psi^4\sqrt{n^T e^{-Q}n}/(8\pi G)$; in the normalization used here,

$$
\boxed{F(n)=2R^2\psi(Rn)^4\sqrt{n^Te^{-Q(Rn)}n}.}                 \tag{1}
$$

Indeed, the tangential determinant of $e^Q$ is its normal cofactor,
$\det(e^Q)n^Te^{-Q}n=n^Te^{-Q}n$, because $\operatorname{tr}Q=0$.
This retains radial--tangential mixing in the induced metric. The density
is positive on the positive conformal-solution branch. The sphere and its
angular labels are held fixed in taking phase-space derivatives.

Write, as in [reduced-vertices.md](reduced-vertices.md),

$$
Q=\sum_\alpha q^\alpha e_\alpha,\qquad
S=\sum_\alpha\frac{p_\alpha}{g_\alpha}e_\alpha,\qquad
u_2=\mathcal C(S:S-R_2[Q]),\qquad
\mathcal C=\frac18(-\Delta)^{-1}.
$$

Here the ten chosen polarizations are Gram-orthogonal,
$\langle e_\alpha,e_\beta\rangle=g_\alpha\delta_{\alpha\beta}$,
and $\psi=1+u_2+O(\Phi^3)$. Thus

$$
\begin{aligned}
F^{(0)}&=2R^2,\\
F^{(1)}&=-R^2 Q_{nn},\\
F^{(2)}&=2R^2\left[
4\mathcal C(S:S-R_2[Q])
+\frac14(Q^2)_{nn}-\frac18 Q_{nn}^2\right]_{r=R}.
\end{aligned}                                                    \tag{2}
$$

The pointwise $R_2$, rather than its integrated replacement by a
$\Gamma\Gamma$ density, is required inside $\mathcal C$. Equation (2)
includes both the direct boundary metric and the nonlocal momentum
contribution. Retaining only the monopole of $u_2$ loses precisely the
angular information needed for local matching.

For any fixed real smearing $f$, define $F_f=\int fF\,d\Omega$.
It is useful to write

$$
F_f=F_f^{(0)}+l_f^Tq+q^TV_fq+p^TK_fp+O(\Phi^3).                 \tag{3}
$$

These matrices are specified, without a fitted angular ansatz, by

$$
\begin{aligned}
(l_f)_\alpha&=-R^2\int d\Omega\,f(e_\alpha)_{nn},\\
(K_f)_{\alpha\beta}
&=\frac{8R^2}{g_\alpha g_\beta}
\int d\Omega\,f\,[\mathcal C(e_\alpha:e_\beta)](Rn),\\
(V_f)_{\alpha\beta}
&=2R^2\int d\Omega\,f\left[
-4\mathcal C R_{\alpha\beta}
+\frac18\bigl(e_\alpha e_\beta+e_\beta e_\alpha\bigr)_{nn}
-\frac18(e_\alpha)_{nn}(e_\beta)_{nn}\right]_{r=R},
\end{aligned}                                                    \tag{4}
$$

where $R_2[Q]=q^\alpha q^\beta R_{\alpha\beta}$ defines the symmetric
polarization of the pointwise second-order curvature. Equations (2)--(4)
are concrete leading local observables for all ten modes.

At this order the scalar harmonics have $\ell\leq4$: a linear spin-two
piece and products of two spin-two pieces. The magnetic linear normal
contraction vanishes. Higher orders generate higher dependent harmonics;
there is no new cutoff of the reconstructed area at $\ell=2$ or $4$.

## 2. The exact constraint algebra and its classification

Extend the bulk phase space by a local pair with the stated orientation,

$$
\Theta_{\rm ext}=\sum_\alpha p_\alpha dq^\alpha+
\int d\Omega\,a(n)d\eta(n),\qquad
\{\eta(n),a(n')\}=\delta_{S^2}(n,n').
$$

The field $a$ is a density relative to the fixed $d\Omega$ convention.
Only its pointwise value must be positive; its nonconstant harmonic
coefficients need not be positive. Define

$$
C(n)=a(n)-F(n),\qquad C_f=a_f-F_f.
$$

Since the bulk and independent corner variables initially commute,

$$
\boxed{\{C_f,C_g\}_{\rm ext}=\{F_f,F_g\}_B=:B_{fg}(q,p).}        \tag{5}
$$

This bracket is not proportional to $C$ merely because $a=F$ is imposed.
The constraint graph contains every bulk point in the candidate's local
existence domain. Consequently, first-classness of all these constraints
requires $B_{fg}=0$ throughout that domain. Even a bulk relation
$\{F_A,F_B\}=f_{AB}{}^D F_D$ would not suffice with independent Abelian
$a_A$: its restriction to the graph is $f_{AB}{}^D a_D$, not zero.

For a finite set of $m$ real angular test modes use dual coefficient
coordinates so that $\{\eta^A,a_B\}=\delta^A_B$. This is an auxiliary
finite description of the constraint algebra, not a new cap on the
reconstructed geometry. Put $B_{AB}=\{F_A,F_B\}_B$.
The constraints are independent because their derivatives with respect
to $a_A$ form the identity. At a constant-rank point, if
$\operatorname{rank}B=r$, then $r$ is even, $r$ combinations are second
class, and $m-r$ combinations are first class. Hence, kinematically,

$$
\dim\mathcal P_{C\text{-only}}=(20+2m)-r-2(m-r)=20+r.           \tag{6}
$$

The full system need not be purely second class. In particular,
$B=DF\,J\,(DF)^T$ has rank at most 20 for this bulk candidate, even if
arbitrarily many angular labels are used. The corresponding continuum
kernel is a finite-rank kernel at each bulk point. Rank can change, and
at the Minkowski origin it vanishes; classification there cannot establish
first-classness on a neighborhood.

The pullback to the matching graph is

$$
\Omega_M=\Omega_B+\sum_A dF_A\wedge d\eta^A.                     \tag{7}
$$

Its null vectors have the form

$$
\delta\eta^A=\lambda^A,\qquad
\delta z=-\{z,F_A\}_B\lambda^A,\qquad B_{AB}\lambda^B=0,        \tag{8}
$$

where $z=(q,p)$. This displays exactly which combined bulk/boost motions
are gauge. For unrestricted $\lambda$, the putative matching generator
instead obeys

$$
\delta_\lambda C_A=\{C_A,C_B\}\lambda^B=B_{AB}\lambda^B,         \tag{9}
$$

so it is not tangent to the constraint surface. The fact that
$\{\eta^A,C_B\}=\delta^A_B$ does not repair (9).

Pure shifts of $\eta$, generated by the independent charges $a_A$, are
a different matter: they preserve $C=0$ but are generally **physical
surface motions**, not null vectors of (7). This distinction prevents
confusing a surface symmetry with a matching gauge transformation.

If $B$ is invertible, reduction by $C$ alone gives, for bulk functions,

$$
\{f,g\}_{C}=
\{f,g\}_B-\{f,F_A\}_B(B^{-1})^{AB}\{F_B,g\}_B.                 \tag{10}
$$

In particular $\{F_A,F_B\}_{C}=0$. The retained edge variables supply
extra physical directions, and the old bulk bracket has changed.
Thus one cannot both keep the original bulk Poisson structure and claim
that $C$ alone removes every extended pair.

## 3. An explicit nonzero bracket in the ten-mode model

The obstruction is not only a generic warning about momentum dependence.
For the electric profile generated by an STF label $E$, the TT differential
operator of the candidate gives, with $t=r^2$,

$$
\begin{aligned}
T^E_{ij}(x)=e^{-t/2}\bigl[&
(t^2-8t+6)E_{ij}
+(6-t)(x_i(Ex)_j+x_j(Ex)_i)\\
&+\tfrac12 x_ix_j(x^TEx)
+\tfrac12\delta_{ij}(t-8)(x^TEx)\bigr].
\end{aligned}                                                    \tag{11}
$$

Its radial contraction simplifies exactly to

$$
(T^E)_{nn}(Rn)=6e^{-R^2/2}n^TEn.                                \tag{12}
$$

The magnetic contraction is zero: normal contraction is a rotation-
equivariant map to scalar spin two, while the magnetic profile has the
opposite spatial parity to scalar $\ell=2$.

Choose the electric axis label
$E=\operatorname{diag}(-1/2,-1/2,1)$, denote its tensor by $T$, and use
an orthogonal polarization basis containing this direction. Set

$$
q=0,\qquad p_E=g\varepsilon,\qquad p_{\alpha\ne E}=0,
\qquad g=\int T:T=\frac{567}{16}\pi^{3/2}.
$$

Let $u_{20}(R)$ be the spherical mean of $\mathcal C(T:T)$.
Using the source monopole already derived in the candidate,

$$
\begin{aligned}
s_0(r)&=e^{-r^2}
\left(54-72r^2+48r^4-\frac{48}{5}r^6+\frac35r^8\right),\\
u_{20}(R)&=\frac18\left[
\frac1R\int_0^R r^2s_0(r)\,dr+
\int_R^\infty rs_0(r)\,dr\right]>0.                             \tag{13}
\end{aligned}
$$

The positivity follows also directly from the positive Newton kernel
and the nonnegative source $T:T$. Here the subscript $0$ denotes the
spherical mean, not a unit-normalized $Y_{00}$ coefficient.

The total charge has no linear term, and its quadratic restriction is

$$
F_1=8\pi R^2+32\pi R^2u_{20}(R)\varepsilon^2+O(\varepsilon^4).
                                                                    \tag{14}
$$

The subscript $1$ on $F_1$ means the smearing
$f=1$, not perturbation order. Rotational invariance and parity imply
that at this point its leading momentum gradient has only the chosen
electric component:

$$
\partial_{p_E}F_1=
\frac{64\pi R^2u_{20}(R)}{g}\varepsilon+O(\varepsilon^3).
$$

Equations (2), (12) and (14) therefore give the local-to-total bracket

$$
\boxed{
\{F(n),F_1\}_B\big|_{q=0,p_E=g\varepsilon}
=-\frac{384\pi R^4u_{20}(R)}{g}
 e^{-R^2/2}\varepsilon\,n^TEn+O(\varepsilon^3).
}                                                               \tag{15}
$$

The remainder is odd in momentum because the reconstructed metric is
exactly even under $p\mapsto-p$. For any fixed $R>0$, (15) is nonzero
at sufficiently small nonzero $\varepsilon$ and angles with $n^TEn\ne0$.
It can be stated entirely with smooth smearings: for $f(n)=n^TEn$,
$\int f^2d\Omega=4\pi/5$, so the leading coefficient of
$\{F_f,F_1\}$ is $4\pi/5$ times the coefficient in (15).
This is an explicit nonzero entry of the constraint matrix on its
constraint surface, and proves that arbitrary local matching in the
original product extension is not first class.

At large radius,

$$
u_{20}(R)=\frac{g}{32\pi R}
+\text{Gaussian-tail corrections},
$$

so (15) becomes $-12R^3e^{-R^2/2}\varepsilon\,n^TEn$ to leading
far-sphere accuracy. It is small but is not zero at finite radius.
This does not invalidate the earlier approximation for the **single**
uniform constraint, whose bracket with itself is identically zero.
It invalidates extending that approximation to an exact local gauge algebra.

### If direct Gaussian boundary terms are deliberately neglected

At fixed finite radius, the leading general bracket from (3) is

$$
B_{fg}=2(l_f^TK_g-l_g^TK_f)p+O(\Phi^2).                         \tag{16}
$$

At the next order one must include $\{F^{(1)},F^{(3)}\}$ as well
as $\{F^{(2)},F^{(2)}\}$. Dropping cubic area data at finite $R$
is not a complete computation of the degree-two bracket.

In the separate far-sphere approximation that discards all direct Gaussian
boundary terms, $l_f$ vanishes, and the leading bracket becomes

$$
B_{fg}=4q^T(V_fK_g-V_gK_f)p+O(\Phi^3),                         \tag{17}
$$

with the nonlocal multipoles in (4) still retained. All local quadratic
charges commute only if

$$
V_fK_g=V_gK_f\quad\hbox{for every pair }f,g.                    \tag{18}
$$

For total area the earlier leading relation gives
$K_1=RG^{-1}$ and $V_1=RG\omega^2/4$, where
$G=\operatorname{diag}g_\alpha$ and
$\omega^2=\operatorname{diag}(11/2,\ldots,11/2,13/2,\ldots,13/2)$.
Already its bracket with a single nonconstant smearing requires

$$
V_fG^{-1}=\frac14G\omega^2K_f.                                 \tag{19}
$$

Equivalently, in the normalized canonical oscillator coordinates, write
$F_f^{(2)}=\mathsf Q^T\widetilde V_f\mathsf Q+
\mathsf P^T\widetilde K_f\mathsf P$. Equation (19) reads

$$
\widetilde V_f=\omega^2\widetilde K_f,
\qquad [\omega^2,\widetilde K_f]=0.                              \tag{20}
$$

The second condition follows from symmetry of both matrices; in particular,
it excludes electric--magnetic entries of $\widetilde K_f$ between the two
different free frequencies. These are explicit tests on the Green-function
matrices (4), not identities implied by TT transversality. No evaluated
claim about every matrix in this tail-discarding approximation is made here.
Equation (15) already settles the exact finite-radius classification.
If all angular multipoles as well as the Gaussian tails are discarded,
every smearing reduces to a constant plus a multiple of $H_2$ and commutes;
that is the uniform approximation again, not local matching.

## 4. What imposing eta = 0 actually does

For the original constraints, introduce $\chi=(C_A,\eta^A)$.
Their full bracket matrix and inverse are

$$
\mathcal D=
\begin{pmatrix}B&-I\\ I&0\end{pmatrix},\qquad
\mathcal D^{-1}=
\begin{pmatrix}0&I\\-I&B\end{pmatrix}.                          \tag{21}
$$

They are invertible for **every** $B$, with determinant one. On
$C=0,\eta=0$,

$$
\Theta_{\rm red}=\sum_\alpha p_\alpha dq^\alpha,
\qquad \Omega_{\rm red}=\Omega_B,
\qquad \{f(q,p),g(q,p)\}_{\chi}=\{f,g\}_B.                     \tag{22}
$$

The last equality follows from the vanishing $CC$ block of
$\mathcal D^{-1}$ and $\{f,\eta\}=0$. It is exact, including all nonlinear
area dependence. In particular, after this reduction

$$
\{a_A,a_B\}_{\chi}=B_{AB},\qquad a_A=F_A(q,p).                  \tag{23}
$$

The areas are now noncommuting bulk observables, not independent commuting
boost momenta. At finite $m$, the second-class phase-space measure has
$\sqrt{\det\mathcal D}=1$; integrating
$\delta(C)\delta(\eta)$ over the canonical edge variables returns the bulk
Liouville measure. This is a finite-dimensional statement, not a derivation
of a continuum functional determinant or of a quantum configuration-space
measure for gravity.

If $B=0$ identically, this is the usual matching quotient and a valid
local gauge section, as in the uniform example. If $\operatorname{rank}B=r>0$,
only $m-r$ of the boost conditions fix gauge. The other $r$ conditions
remove physical directions of the $20+r$ dimensional system (6).
Thus (22) is a justified symplectic elimination **conditional on choosing
no independent physical corner sector**. It is not the gauge reduction
of $C$ alone.

Dynamics must respect this distinction. A bare bulk Hamiltonian obeys
$\{C_A,H_B\}=-\{F_A,H_B\}_B$ and need not preserve matching.
For $C$ alone, stabilization with multipliers of $C$ requires
$B_{AB}\lambda^B=\{F_A,H_B\}$; null directions impose consistency
conditions unless the Hamiltonian is extended. For the complete
second-class set (21), the Dirac evolution is unambiguous and reduces
exactly to $H_B$. For example, on the section an extended total Hamiltonian
can use the term $-\eta^A\{F_A,H_B\}$ to enforce preservation of $C$.
One must not infer preservation from the leading special case
$F_1=\mathrm{constant}+RH_2$.

## 5. An exact local first-class conversion, and its price

There is no mathematical theorem that extra *pairs beyond the already
introduced area--boost pairs* are unavoidable. One can instead change the
off-section constraints. The following explicit construction demonstrates
both the possibility and what has changed.

For finitely many test modes choose reference constants $F_A^\star$ and
form the auxiliary Hamiltonian

$$
W=\eta^A(F_A(z)-F_A^\star),\qquad
D_W h=\{h,W\}_{\rm ext}.
$$

Let $\mathcal T=\exp D_W$ be its local canonical flow acting on functions.
Define

$$
\widetilde C_A=\mathcal T(a_A-F_A^\star).
$$

Because a canonical transformation preserves brackets,

$$
\boxed{\{\widetilde C_A,\widetilde C_B\}=0,
\qquad \{\eta^A,\widetilde C_B\}=\delta^A_B.}                  \tag{24}
$$

The second identity holds because $\mathcal T\eta=\eta$.
Expanding the flow gives

$$
\widetilde C_A=
 a_A-F_A-\frac12\eta^B B_{AB}
 -\frac16\eta^B\eta^D\{B_{AB},F_D\}_B+O(\eta^3).               \tag{25}
$$

For example, the two brackets of $a$ with the linear-in-$\eta$ correction
contribute $-B_{AB}$ and cancel (5). Higher terms are fixed by the flow,
not guessed counterterms. Smooth local flow suffices for the exact
statement (24); convergence of a written Lie series need not be assumed.
The reference constants permit working near a positive-area background;
only local flows staying in the allowed neighborhood are asserted.

The section $\eta=0$ intersects each sufficiently small converted gauge
orbit once, imposes $a=F$, and reduces to exactly (22). Gauge-invariant
extensions of the bulk observables and Hamiltonian are

$$
Z=\mathcal T z,\qquad
\widetilde H=\mathcal T H_B
=H_B-\eta^A\{F_A,H_B\}_B+O(\eta^2),
$$

since they commute with every $\widetilde C_A$. On the section $Z=z$
and $\widetilde H=H_B$.

**The price is substantive:** away from $\eta=0$, (25) no longer imposes
$a_A=F_A(z)$. It identifies $a_A$ with a boost-dependent expression instead.
This is a local parametrization of the same twenty-dimensional bulk
system, not an equivalence to the original $C$-only system (6), and not a
proof that these converted translations are physical normal-frame boosts.
The gauge-invariant geometric areas are $F_A(Z)$; their mutual bracket
is still the dressed version of $B_{AB}$. Independent $a_A$ are not those
gauge-invariant observables.

Equivalently, writing $\widetilde C_A=a_A-\mathcal F_A(\eta,z)$,
first-class Abelian matching demands the flatness equation

$$
\partial_{\eta^A}\mathcal F_B-
\partial_{\eta^B}\mathcal F_A+
\{\mathcal F_A,\mathcal F_B\}_B=0,\qquad
\mathcal F_A(0,z)=F_A(z).                                      \tag{26}
$$

The undressed choice $\mathcal F_A=F_A$ fails exactly by $B_{AB}$.
The local canonical flow above solves (26). The continuum use of an
arbitrary smooth $\eta(n)$ additionally needs functional regularity;
the finite-mode proof does not establish that analytic limit.

## 6. What follows for quantum positivity and transfer

A product of commuting independent area representations does not solve
noncommuting bulk matching. Formally,

$$
[\widehat C_f,\widehat C_g]=
[\widehat F_f,\widehat F_g],                                   \tag{27}
$$

assuming the initially independent sectors commute. For Weyl-ordered
quadratic area approximants the commutator is exactly $i\hbar$ times the
quantized Poisson bracket of those approximants. Thus the leading
obstruction is already present without invoking an ordering anomaly.
Simultaneous strong equations $\widehat C_f\Psi=0$ entail additional
commutator conditions; a common generalized kernel is not guaranteed.
The uniform single-constraint rigging map cannot be multiplied pointwise
into a local positive averaging prescription. There is no joint area
spectral measure unless strong commutation is established.

If one **chooses** either the second-class elimination (21) or the
converted gauge reduction (24), reduction first leaves the original bulk
canonical phase space. The previously specified free quantization then
has its positive ten-oscillator Fock space and its free Hamiltonian
$\widehat H_2$. Local area observables must be quantized on that space as
functions of the same oscillators; no independent edge degeneracy is
created. The exact nonlinear construction has only a local classical
chart, so this statement does not supply its global interacting quantum
domain or an everywhere positive local-area operator.

In particular, classical positivity of (1) does not by itself prove
positivity, self-adjointness, or a unique ordering of a quantum operator
for (1). Nor does the classical conversion automatically yield a
unitary equivalence of quantizations on positive-area half-lines.
The earlier half-line obstruction to an unrestricted self-adjoint boost
coordinate remains relevant to any attempt at quantizing the extension
before reduction.

For any reduction whose quantum completion supplies a positive Hilbert
space and a self-adjoint Hamiltonian bounded below, the same conditional
reconstruction statement remains valid:

$$
T(\tau)=e^{-\tau(\widehat H-E_0)/\hbar},\qquad
U(t)=e^{-it(\widehat H-E_0)/\hbar},
$$

$$
\sum_{i,j}\langle v_i,T(\tau_i+\tau_j)v_j\rangle
=\left\|\sum_jT(\tau_j)v_j\right\|^2\geq0.                      \tag{28}
$$

For the free reduced model these statements are already available.
For the interacting model they remain conditional on the operator
construction; local corner matching has not established them.
A nonconstant boost rate would involve $F_\nu=\int\nu Fd\Omega$,
which generally neither commutes with $H_B$ nor with all other local
areas. Its exponentials cannot be assumed diagonal in the uniform-area
basis or factorized into independent angular transfer kernels. For a
fixed self-adjoint, lower-bounded generator the spectral theorem still
applies, irrespective of that lack of simultaneous diagonalization.

No finite microscopic entropy follows. The unreduced boundary field
sector, the $C$-only system with surviving edge directions, and the
no-independent-edge reduction are different theories of the physical
corner content. Even the last has an infinite-dimensional free Fock
space. Its leading finite total-energy/area-window counts remain counts
of bulk oscillator states under the earlier specified approximation;
they are not local corner multiplicities. Noncommuting local area
operators also preclude assigning simultaneous sharp area bins by
multiplying spectral projectors without an additional measurement or
constraint prescription.

## 7. Relation to genuine gravitational edge charges and remaining input

Donnelly and Freidel derive surface diffeomorphisms and a local
$SL(2,\mathbb R)$ normal-frame algebra in an extended continuum phase
space, with charges obeying their Eqs. (3.66)--(3.70). They distinguish
the area Casimir from the single boost generator obtained after fixing
the normal metric; the corresponding quantum spectra need not agree.
Their concise statement is: “Thus we recover Carlip and Teitelboim's
result as a gauge fixing of ours.”
[Local subsystems in gauge theory and gravity, Section 3.6](https://www.alphaxiv.org/pdf/1601.04744v2)
The formulae and qualifications in Sections 2--5 above are new deductions
for the present finite candidate, not results claimed by that paper.

There is no contradiction between a commuting unreduced metric area
and (15): solving the constraints makes the metric depend on momenta,
and restriction to a symplectic finite-dimensional submanifold is not
in general a Poisson map for arbitrary observables. Likewise, the
continuum area Casimir algebra cannot simply be imported as the bracket
of the fixed-sphere functions (1).

If one insists on genuine local gravitational matching, with geometric
area fixed in every frame rather than only on a chosen section, the
present candidate must acquire the missing boundary symplectic data.
For example, in a product construction with non-Abelian edge charges
$J_A$ and bulk charges $Q_A$, matching their difference requires an
oppositely oriented edge algebra so that schematically

$$
\{J_A,J_B\}=-f_{AB}{}^D J_D,\qquad
\{Q_A,Q_B\}=f_{AB}{}^D Q_D
\quad\Longrightarrow\quad
\{J_A-Q_A,J_B-Q_B\}=-f_{AB}{}^D(J_D-Q_D).                       \tag{29}
$$

This example assumes a common, closed Lie charge algebra and commuting
sectors. For field-dependent structure functions, additional consistency
and cross terms must be supplied. Merely adjoining inert edge momenta
does not cancel $B$: one must derive their charges, cross brackets and
matching conditions, or a modified mixed symplectic potential. The
nonlocal finite-rank kernel $B_{fg}(q,p)$ in (5) has not been identified
with the continuum normal-frame or surface-diffeomorphism algebra.

What is still undetermined is therefore precise: the boundary embedding
and normal-frame variables of this candidate, which surface transformations
are gauge after gluing, their closed charge algebra and Hamiltonian,
whether elimination (21) or a gravitational realization of (24) is the
appropriate physical reduction, and the positive interacting operator
representation of that choice. Equations (15), (21) and (24) settle the
local algebraic alternatives without resolving those gravitational inputs
or replacing them with a guessed finite state count.
