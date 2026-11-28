# Quantizing the uniform corner mode: states, gluing, and Euclidean evolution

## Scope and result

Within the extended area--boost polarization specified in the
[corner-action analysis](corners-and-reconstruction.md), the spatially
constant positive-area mode admits an explicit positive Hilbert space,
self-adjoint area and dilation operators, unitary boost evolution, and
an exactly matching positive Euclidean transfer semigroup. Normal-boost
matching across two sides can also be solved by group averaging.

This derives a quantum realization of that boundary sector. Its area
spectrum is continuous and its Hilbert space remains infinite-dimensional
inside every area interval of nonzero width. Its Euclidean transfer
operator has no finite Gibbs trace. These are results of the specified
quantization, not a microscopic finite-state derivation.

The full gravitational corner algebra contains more than this single
mode. We also prove that its undeformed normal-frame $SL(2,\mathbb R)$
subalgebra has no nontrivial finite-dimensional unitary representation.
An eventual finite-state theory must explain which reduction, deformation,
or limiting representation replaces that requirement.

Neither conclusion supplies the missing bulk constraint reduction,
the one-loop amplitude, or reflection positivity of the archived
simulation measure. The full project objective remains open.

## 1. Boundary data being quantized

Choose the corner orientation, null-normal normalization and polarization
as in the previous report. For the uniform mode,

$$
a=\frac{A}{8\pi G}>0,\qquad
\Theta_C=a\,\delta\eta,\qquad
\Omega_C=\delta a\wedge\delta\eta,\qquad
\{\eta,a\}=1.
$$

The Lorentzian normal boost has parameter $\eta\in\mathbb R$. For a
prescribed constant positive boost rate $\nu$, the boundary charge is
$H_C=\nu a$, giving $\dot a=0$ and $\dot\eta=\nu$. This is a boundary
observer's evolution, not an already derived Hamiltonian of a closed
gravitational diamond.

The source of the area--boost pairing and its boundary conventions is
[Carlip--Teitelboim, PDF pp. 3--5](https://www.alphaxiv.org/pdf/gr-qc/9312002v3?page=4).
The null-joint normalization is fixed by
[Lehner et al., PDF pp. 24--25](https://www.alphaxiv.org/pdf/1609.00207?page=25).

There is a further distinction at the quantum level: in an extended
normal-frame treatment, area is related to the $SL(2,\mathbb R)$ Casimir.
Fixing the normal metric at the classical level leaves a single boost
whose charge is area. Quantizing that reduced pair need not give the same
spectrum as quantizing the entire normal-frame algebra and then reducing.
[Donnelly--Freidel, PDF pp. 31--32](https://www.alphaxiv.org/pdf/1601.04744v2?page=31)
make this distinction explicit. We quantize the reduced uniform pair here;
we do not identify it with their complete surface Hilbert space.

## 2. A positive-area quantum representation

Introduce a reference area charge $a_*>0$ solely to make a logarithm
dimensionless:

$$
q=\log(a/a_*),\qquad d=-a\eta,\qquad
\Omega_C=\delta d\wedge\delta q,\qquad \{q,d\}=1.
$$

Use the regular Schrödinger representation of this canonical pair on
$L^2(\mathbb R,dq)$. The unitary change of representation

$$
\psi(a)=a^{-1/2}\phi(\log(a/a_*))
$$

gives the equivalent positive-area Hilbert space

$$
\mathcal H_C=L^2(\mathbb R_+,da),\qquad
(\hat a\psi)(a)=a\psi(a),\qquad
\hat d=-i\hbar\left(a\frac{d}{da}+\frac12\right).
$$

The multiplication operator $\hat a$ is self-adjoint on
$\{\psi:\int_0^\infty a^2|\psi(a)|^2da<\infty\}$.
The closure of $\hat d$ on $C_c^\infty(\mathbb R_+)$ is self-adjoint
because it is unitarily equivalent to $-i\hbar\,d/dq$ on the real line.
On this common core,

$$
[\hat d,\hat a]=-i\hbar\hat a,\qquad
\widehat A=8\pi G\,\hat a.
$$

Generalized area states obey
$\langle a|a'\rangle=\delta(a-a')$; normalizable states are square-integrable
area wave packets. This minimal representation has multiplicity one.
Additional shape, normal-frame, matter and bulk labels have not been
quantized or removed by this construction.

The formal conjugate $\hat\eta=i\hbar\,d/da$ is symmetric on
$C_c^\infty(\mathbb R_+)$ but has unequal deficiency indices on the
half-line and no self-adjoint extension there. Thus we do not assert
an everywhere self-adjoint angle obeying unrestricted translations of
a positive area spectrum. Boost evolution is nevertheless well defined:
its generator is the self-adjoint $\hat a$, not $\hat\eta$.

This is a specified canonical quantization, not a proof that every possible
quantum theory of gravitational area has a continuous spectrum. Changing
the reduced observable algebra or the quantization assumptions changes
that question.

## 3. Finite dimension and the normal-frame symmetry

Two exact statements constrain a finite-state proposal.

**Affine pair.** If finite-dimensional matrices with positive
$\hat a\ne0$ obey $[\hat d,\hat a]=-i\hbar\hat a$, the trace gives
$0=-i\hbar\operatorname{Tr}\hat a$, a contradiction. This proof does not
require a self-adjoint angle operator.

More quantitatively, if a finite-dimensional area-window model has
$\hat a\geq a_{\min}I$ with $a_{\min}>0$ and

$$
[\hat d,\hat a]=-i\hbar\hat a+\Delta ,
$$

then

$$
\operatorname{Tr}\Delta=i\hbar\operatorname{Tr}\hat a,\qquad
\|\Delta\|\geq
\hbar\,\frac{\operatorname{Tr}\hat a}{\dim\mathcal H}
\geq\hbar a_{\min}.
$$

Consequently the exact affine relation cannot be recovered uniformly in
operator norm by such finite-dimensional models while that positive area
window is held fixed. Matrix elements on a chosen interior set of states
can still converge; the bound does not rule out that weaker limit.

**Full normal-frame algebra.** Rescale the generators of
Donnelly--Freidel's Eq. (3.70) to self-adjoint operators obeying

$$
[J_0,J_1]=i\hbar J_2,\qquad
[J_0,J_2]=-i\hbar J_1,\qquad
[J_1,J_2]=-i\hbar J_0 .
$$

For finite matrices, cyclicity of the trace implies

$$
-i\hbar\operatorname{Tr}J_0^2
=\operatorname{Tr}([J_1,J_2]J_0)
=\operatorname{Tr}(J_1[J_2,J_0])
=i\hbar\operatorname{Tr}J_1^2 .
$$

Positivity of the Hilbert inner product forces $J_0=J_1=0$, and then
$J_2=0$. There is no nontrivial finite-dimensional unitary realization of
this undeformed Lie algebra. Nonzero geometric area cannot be represented
through its nonzero Casimir in that finite realization.

The premise is the normal-frame algebra derived in the extended continuum
phase space, not a claim that a finite regulator must retain it exactly.
An angular cutoff alone does not remove its spatially constant subgroup.
The single boost group $\mathbb R$ by itself *can* act unitarily on a finite
Hilbert space; it is the full affine or $SL(2,\mathbb R)$ relations that
produce these obstructions.

A physical quotient may remove some extended corner observables, and a
microscopic theory may deform their algebra. Those are possible responses,
but the quotient or deformation and its dynamics must be derived. The
proof does not forbid finite physical state sectors after such a reduction.

## 4. Solve the uniform normal-boost matching constraint

Give the two sides of an interface opposite symplectic orientations:

$$
\Omega_{LR}=\delta a_L\wedge\delta\eta_L
           -\delta a_R\wedge\delta\eta_R,\qquad
C=a_L-a_R.
$$

The constraint generates the common shift of $\eta_L$ and $\eta_R$.
On $C=0$, the reduced symplectic form is
$\delta a\wedge\delta(\eta_L-\eta_R)$.
Matching this charge removes one pair and leaves a relative boost pair.

Quantize the two sides using a dual representation for the oppositely
oriented side. On $L^2(\mathbb R_+^2,da_Lda_R)$ the matching operator is
multiplication by $a_L-a_R$. Its kernel contains no nonzero square-integrable
functions: the diagonal has measure zero. Treating matching as an ordinary
orthogonal projector would therefore give the wrong answer.

Instead choose test functions in $C_c^\infty(\mathbb R_+^2)$ and average
the unitary matching action distributionally:

$$
\int_{\mathbb R}\frac{d\lambda}{2\pi\hbar}
e^{-i\lambda(a_L-a_R)/\hbar}=\delta(a_L-a_R).
$$

The resulting physical pairing is

$$
\langle\phi,\psi\rangle_{\rm match}
=\int_0^\infty da\,\overline{\phi(a,a)}\,\psi(a,a)\geq0
\quad\text{when }\phi=\psi.
$$

Quotient by test functions vanishing on the diagonal and complete. The
result is again $L^2(\mathbb R_+,da)$. This is a solved Abelian matching
constraint with a positive physical pairing, not the full gravitational
constraint problem. The general distinction between an ordinary
projector and a rigging map is described by
[Marolf, PDF pp. 2--3](https://www.alphaxiv.org/pdf/gr-qc/0011112v1?page=3);
the calculation above establishes positivity directly for this case.

For oppositely oriented side evolutions
$H_{LR}=\nu_La_L-\nu_Ra_R$, the induced charge is
$(\nu_L-\nu_R)a$. Equal common boost rates are generated by the constraint
and act trivially after matching. A nonzero boundary evolution therefore
requires a specified relative boost rate, not the common gauge motion.

The identity kernel composes associatively:
$\int da_2\,\delta(a_1-a_2)\delta(a_2-a_3)=\delta(a_1-a_3)$ on test
functions. Shape matching, local surface diffeomorphisms, non-Abelian
normal-frame matching, bulk constraints, and their determinants are absent.
This therefore does not prove gluing of complete gravitational diamonds.

## 5. Exact Lorentzian and Euclidean evolution of this charge

For the specified constant $\nu>0$ (the relative rate on a matched pair),

$$
H_C=\nu\hat a,\qquad
U_C(t)=e^{-it\nu\hat a/\hbar},\qquad
T_C(\tau)=e^{-\tau\nu\hat a/\hbar},\quad \tau\geq0 .
$$

$U_C$ is strongly continuous and unitary. $T_C$ is a strongly continuous
positive contraction semigroup. Their kernels in the same area basis are

$$
K_L(a',a;t)=\delta(a'-a)e^{-it\nu a/\hbar},\qquad
K_E(a',a;\tau)=\delta(a'-a)e^{-\tau\nu a/\hbar}.
$$

The same $da$ gluing measure gives
$K_E(\tau)K_E(\sigma)=K_E(\tau+\sigma)$.
An area-window projector commutes with this particular $H_C$, so restriction
to such a window preserves its semigroup property. The Hilbert space in
that window is still infinite. Generic bulk occupation projections need
not commute with the physical Hamiltonian, as shown in the previous report.

For an explicit normalized area packet
$\psi_\lambda(a)=\sqrt{2\lambda}\,e^{-\lambda a}$, $\lambda>0$,

$$
C_E(\tau)=\langle\psi_\lambda,T_C(\tau)\psi_\lambda\rangle
=\frac{2\lambda}{2\lambda+\nu\tau/\hbar},
\qquad
C_L(t)=\frac{2\lambda}{2\lambda+i\nu t/\hbar}.
$$

This gives an exact analytic continuation and a positive-time reflection
kernel:

$$
\sum_{i,j}\bar c_i c_j C_E(\tau_i+\tau_j)
=\int_0^\infty 2\lambda e^{-2\lambda a}
\left|\sum_jc_j e^{-\nu\tau_j a/\hbar}\right|^2da\geq0.
$$

For distinct $\tau_i\geq0$, write $x_i=\lambda+\nu\tau_i/\hbar$.
The kernel matrix is strictly positive; its determinant is the Cauchy
determinant

$$
\det\left[\frac{2\lambda}{x_i+x_j}\right]_{i,j=1}^n
=(2\lambda)^n
\frac{\prod_{i<j}(x_i-x_j)^2}{\prod_{i,j}(x_i+x_j)}>0 .
$$

These are exact operator and integral identities for the derived corner
charge. They do not define stationary vacuum Schwinger functions for the
full gravity theory or a local Euclidean Monte Carlo measure. In particular,
$H_C$ has infimum zero but no normalizable zero-energy eigenvector in this
representation. A vacuum-based reconstruction would need more structure.

## 6. An area bound does not produce a finite state count or Gibbs trace

Let $I=[a_1,a_2]\subset\mathbb R_+$ have $\Delta a=a_2-a_1>0$. The states

$$
\phi_n(a)=\frac{\mathbf1_I(a)}{\sqrt{\Delta a}}
\exp\left(\frac{2\pi in(a-a_1)}{\Delta a}\right),\qquad n\in\mathbb Z,
$$

form an infinite orthonormal family inside the area window. All have
bounded energy under $H_C=\nu a$. Matching two corners as above does not
remove this family.

For every $\tau>0$, each state has the same positive diagonal transfer
matrix element

$$
\langle\phi_n,T_C(\tau)\phi_n\rangle
=\frac{\hbar}{\nu\tau\,\Delta a}
\left(e^{-\nu\tau a_1/\hbar}-e^{-\nu\tau a_2/\hbar}\right)>0.
$$

The sum over $n$ diverges. Hence $T_C(\tau)$ is not trace class, even
after restricting to this area interval. The finite integral
$\int da\,e^{-\beta\nu a/\hbar}$ is not the Hilbert-space Gibbs trace:
it omits the infinite conjugate-coordinate state volume. There is no
derived finite partition function or entropy $\log(2^M)$ or $\log(3^M)$
in this boundary theory.

## Consequence for the shared classical--quantum goal

The uniform corner mode now has an explicit quantum realization and
matching rule, with Euclidean and Lorentzian evolution obtained from the
same charge. That removes an ambiguity in this restricted boundary sector.
It also makes the remaining microscopic question sharper: a finite
gravitational state count requires additional physical reduction or
microscopic dynamics, not just finitely many angular coefficients or an
upper bound on area.

The bulk constraint bracket obstruction remains independent of this corner
calculation. A full construction must supply a dynamically admissible bulk
physical space, its boundary algebra and coupling, its positive quantum
measure, the one-loop scattering observable, and controlled limits.
Neither the kinematical sector above nor its exact transfer semigroup is
substituted for those requirements.
