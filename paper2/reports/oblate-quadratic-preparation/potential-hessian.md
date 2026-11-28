# Spatial-curvature Hessian for the oblate quadratic operator

**Status:** derivation and implementation preparation only. No finite
oblate spectrum, nonlinear closure theorem, or new gravitational experiment is
claimed here.

This note prepares the unit-lapse spatial-curvature block of the stationary
quadratic Hamiltonian on the exact oblate Minkowski background

$$
x=\cos\theta,\qquad d=1-x^2,\qquad
\rho=r^2+\epsilon,\qquad \Sigma=r^2+\epsilon x^2,
$$

$$
\sigma_{ij}dy^idy^j=
{\Sigma\over\rho}dr^2+{\Sigma\over d}dx^2+\rho d\,d\phi^2,
\qquad a=-1,\quad b=V=W=0,\quad p=0,\quad \Lambda=0 .
$$

All radial manipulations below start with smooth compact support on a fixed
shell $r\in I\Subset(0,\infty)$. Boundary terms are recorded rather than
dropped because the physical domain and boundary generator are a separate
input to the full oblate operator assembly.

## 1. Normalization and the block exported here

The Paper I coefficient action is the Galerkin pullback of the
Einstein-Hilbert action. In ADM language near the static unit-lapse,
zero-shift background its normalized action is

$$
\widetilde I={16\pi G\over 4\pi}I,\qquad
\widetilde L_{\mathrm{kin}}
={1\over4}\,
G(T\dot q-\mathcal L_\beta\sigma,
  T\dot q-\mathcal L_\beta\sigma),
$$

with

$$
G(h,k)=\int_I dr\int {d\Omega\over4\pi}\sqrt\sigma\,
\left(h_{ij}k^{ij}-\operatorname{tr}h\,\operatorname{tr}k\right).
$$

The potential whose Hessian is prepared here is

$$
\boxed{V[\sigma]=-\int_I dr\int {d\Omega\over4\pi}
\sqrt\sigma\,R[\sigma].}
$$

The full stationary coordinate Hessian also includes temporal variables
$z=(a,b,W)$, the $T^*GT$ kinetic mass matrix on
$q=(\lambda,X,V,C)$, shift-Schur terms, lapse rows, and the radial
boundary realization. This report supplies only $V''$ on retained spatial
perturbations.

For the shared canonical assembly, define the local spatial mass form
$S=T^*GT$. At the static point and at $V=0$, a shift descriptor
$s=(b,W)$ has vector field $v=(b/\lambda,W)$ and

$$
U_s=\mathcal L_v\sigma,\qquad
\mathcal J s=S^{-1}T^*GU_s,\qquad
A(s,t)=-{1\over2}G(N_\perp U_s,N_\perp U_t).
$$

With the actual normalized canonical momentum $p$ conjugate to $q$, the
stationary quadratic Hamiltonian blocks are organized as

$$
H_2={1\over2}\langle q,V''q\rangle
+\langle p,S^{-1}p\rangle
+\langle p,\mathcal J s\rangle
+\langle\delta a,Lq\rangle
+{1\over2}\langle s,A s\rangle ,
$$

so, in the order $w=(q,p)$,

$$
K_{\mathrm{full}}=\begin{pmatrix}V''&0\\0&2S^{-1}\end{pmatrix},
\qquad
E=\begin{pmatrix}L&0\\0&\mathcal J^*\end{pmatrix},
\qquad
D=\begin{pmatrix}0&0\\0&A\end{pmatrix}.
$$

The lower block $2S^{-1}$ is the kinetic $K$ used in the Paper I oblate
rank derivation; it is not the whole $K_{\mathrm{full}}$.

## 2. Covariant second variation on the flat background

Let $h_{ij}$ and $k_{ij}$ be arbitrary symmetric spatial perturbations, not
necessarily retained. Indices are moved with $\sigma$, $\nabla$ is its
Levi-Civita connection, and

$$
(\delta G[h])_{ij}=
{1\over2}\left(
-\nabla^2h_{ij}-\nabla_i\nabla_jh
+\nabla_i\nabla^a h_{aj}
+\nabla_j\nabla^a h_{ai}
-\sigma_{ij}\nabla^a\nabla^b h_{ab}
+\sigma_{ij}\nabla^2h
\right)
$$

is the flat-background linearized Einstein tensor, with
$h=\sigma^{ij}h_{ij}$. Because $R_{ij}[\sigma]=0$, the bulk first variation
vanishes and the metric-coefficient map has no bulk second-derivative
contribution at this background: if
$g_{ij}(u,v)=\sigma_{ij}+u h_{ij}+v k_{ij}
+uv\,\ell_{ij}+O((u,v)^3)$ then

$$
{d\over ds}V[g(s)]\big|_{s=0}
=\int {\sqrt\sigma\over4\pi}\,G^{ij}[\sigma]\,\dot g_{ij}=0
$$

and the $\ell_{ij}$ term in the mixed second derivative is multiplied by
$G^{ij}[\sigma]=0$. Only explicit boundary terms can feel nonlinear
coefficient-map second derivatives; those terms are assigned to the boundary
domain handoff, not hidden in the bulk operator.

For compact support, the curvature-potential Hessian is therefore

$$
\boxed{
V''_\mathrm{bulk}[h,k]=
\int_I dr\int {d\Omega\over4\pi}\sqrt\sigma\,
h^{ij}(\delta G[k])_{ij}
=
\int_I dr\int {d\Omega\over4\pi}\sqrt\sigma\,
k^{ij}(\delta G[h])_{ij}.}
$$

This sign is fixed by
$\delta V=\int\sqrt\sigma\,G^{ij}\delta\sigma_{ij}$ for
$V=-\int\sqrt\sigma R$ and by the convention
$\delta G^{ij}=\delta R^{ij}-\frac12\sigma^{ij}\delta R$ at $R_{ij}=0$.
For a Cartesian transverse-traceless perturbation with appropriate compact or
periodic integration-by-parts conditions,
$\nabla^i h_{ij}=0$, $h=0$, so

$$
(\delta G[h])_{ij}=-{1\over2}\Delta h_{ij},\qquad
V''_\mathrm{bulk}[h,h]
={1\over2}\int {d^3y\over4\pi}\,\partial_a h_{ij}\partial^a h_{ij}
$$

after compact integration by parts. Combined with
$\widetilde L_{\mathrm{kin}}=\frac14G(\dot h,\dot h)$ and
$G(h,h)=\int h_{ij}h_{ij}$ for TT fields, this gives the expected
positive gradient energy and the standard flat wave speed in the
normalization used here.

## 3. First-derivative weak form and boundary remainder

For implementation the second-derivative operator should be integrated by
parts once. Define

$$
\mathcal B(h,k)=
{1\over2}\nabla_a h_{ij}\nabla^a k^{ij}
-\nabla_i h^{ij}\nabla^a k_{aj}
+{1\over2}\nabla_i h^{ij}\nabla_j k
+{1\over2}\nabla_i k^{ij}\nabla_j h
-{1\over2}\nabla_a h\nabla^a k .
$$

A direct covariant integration by parts gives

$$
\int_M\sqrt\sigma\,h^{ij}(\delta G[k])_{ij}
=\int_M\sqrt\sigma\,\mathcal B(h,k)
+\int_{\partial M}\sqrt\gamma\,n_a\,\Theta^a(h,k),
$$

where

$$
\begin{aligned}
\Theta^a(h,k)=&-{1\over2}h^{ij}\nabla^a k_{ij}
+h^{ai}\nabla^b k_{ib}
-{1\over2}h^{ai}\nabla_i k\\
&-{1\over2}h\,\nabla_i k^{ai}
+{1\over2}h\nabla^a k .
\end{aligned}
$$

The weak bulk form is first order and symmetric in the two real perturbations
pointwise. Complex angular blocks conjugate the first argument. On the initial
compact-support domain the boundary integral vanishes; for the strong bulk
operator the difference of the two Green fluxes controls its adjoint identity. The final
finite-shell realization must choose boundary data for which this remainder
either vanishes or is cancelled by the variation of the chosen boundary
Hamiltonian. Until that choice is made, the report exports both the weak bulk
integrand and $\Theta^a$.

The pure-gradient control is immediate. For $h=\mathcal L_\xi\sigma$ with
compactly supported $\xi$, spatial diffeomorphism invariance of the scalar
curvature functional gives

$$
\delta G[\mathcal L_\xi\sigma]=\mathcal L_\xi G[\sigma]=0,
\qquad
V''_\mathrm{bulk}[\mathcal L_\xi\sigma,k]=0 .
$$

Consequently $V''J_{\ker}=0$ for exactly retained shift-kernel sections:
there $T\mathcal J s=\mathcal L_{v_s}\sigma$ lies in the retained tangent space
without a discarded normal component, so the curvature Hessian annihilates
that column on compact support.

The lapse row has the same sign convention. For a lapse test $f$,

$$
\langle f,Lq\rangle
={1\over2}\int_I dr\int {d\Omega\over4\pi}\sqrt\sigma\,
f\,\delta R[Tq]
={1\over2}G(\nabla^2f,Tq)
$$

on compact tests, equivalently

$$
\langle f,Lq\rangle
=-{1\over2}\int\sqrt\sigma\,
\nabla^i f\left(\nabla^j(Tq)_{ij}-\nabla_i\operatorname{tr}(Tq)\right)
+\text{radial and angular flux}.
$$

The root export should keep that flux beside the $V''$ flux until the same
boundary/domain choice has fixed both adjoints.

## 4. Pullback through the retained spatial map

At fixed $(r,x,\phi)$ the spatial coefficient map used by the oblate
calculation has linear perturbations

$$
h=T_I\,u^I
=\delta\lambda\,dr^2
+\delta X\,q_{S^2}
+2\,\gamma_{AB}\,\delta V^A\,dr\,dy^B
+\delta C_{AB}\,dy^A dy^B ,
\qquad q^{AB}\delta C_{AB}=0 .
$$

The background metric is diagonal; its connection is not. With coordinates
$(r,x,\phi)$ and metric
$(g_{rr},g_{xx},g_{\phi\phi})=(\Sigma/\rho,\Sigma/d,\rho d)$, the nonzero
Christoffel symbols are

$$
\begin{gathered}
\Gamma^r{}_{rr}={r\over\Sigma}-{r\over\rho},\quad
\Gamma^r{}_{rx}={\epsilon x\over\Sigma},\quad
\Gamma^r{}_{xx}=-{\rho r\over d\Sigma},\quad
\Gamma^r{}_{\phi\phi}=-{\rho r d\over\Sigma},\\
\Gamma^x{}_{rr}=-{\epsilon xd\over\rho\Sigma},\quad
\Gamma^x{}_{rx}={r\over\Sigma},\quad
\Gamma^x{}_{xx}={\epsilon x\over\Sigma}+{x\over d},\quad
\Gamma^x{}_{\phi\phi}={\rho xd\over\Sigma},\\
\Gamma^\phi{}_{r\phi}={r\over\rho},\qquad
\Gamma^\phi{}_{x\phi}=-{x\over d}.
\end{gathered}
$$

The $r$-dependence of $T_I$ is part of the derivative:

$$
\nabla_a h_{ij}
=(\partial_aT_{I,ij}
-\Gamma^b{}_{ai}T_{I,bj}
-\Gamma^b{}_{aj}T_{I,ib})u^I
+T_{I,ij}\partial_a u^I .
$$

The shift column uses the same retained spatial map. At $V=0$,
$v=(b/\lambda,W)$ gives

$$
\begin{aligned}
(U_s)_{rr}&=2b'-b\,\partial_r\log\lambda
+W^A\partial_A\lambda,\\
(U_s)_{rA}&=\partial_A b-b\,\partial_A\log\lambda
+\gamma_{AB}(W^B)',\\
(U_s)_{AB}&={b\over\lambda}\partial_r\gamma_{AB}
+(\mathcal L_W\gamma)_{AB}.
\end{aligned}
$$

Thus the $b'$ and $(W^A)'$ pieces are exactly retained tangent injections,
$b'\mapsto\delta\lambda=2b'$ and $(W^A)'\mapsto\delta V^A=(W^A)'$.
The normal response $N_\perp U_s$ has no radial derivatives after these
exact tangent cancellations, but the tangential map still has the form
$\mathcal J s=J_0s+J_1s'$, so $\mathcal J^*$ in the canonical block above must be interpreted
with the chosen radial adjoint.

For one magnetic sector $m\ge0$ and cutoff $L=N-1$, use the archived
associated-Legendre convention

$$
Y_{\ell m}=P_\ell^m(x)e^{im\phi},
\quad
Y^A_{\ell m,\mathrm p}=(d\,\partial_xY_{\ell m},
\,imY_{\ell m}/d),
\quad
Y^A_{\ell m,\mathrm a}=(imY_{\ell m},-\partial_xY_{\ell m}),
$$

with scalar degree $\ell=m,\ldots,L$, vector degree $\ell=\max(1,m),\ldots,L$,
and tensor degree $\ell=\max(2,m),\ldots,L$. A complex $m>0$ block represents
two real blocks:

$$
\operatorname{Re}(A e^{im\phi})=A_c\cos m\phi-A_s\sin m\phi,\qquad
\operatorname{Im}(A e^{im\phi})=A_c\sin m\phi+A_s\cos m\phi .
$$

For each pair of retained coefficient basis tensors $T_I,T_J$, substitute
$h=T_I f(r)$ and $k=T_J g(r)$ into the weak bulk form. Exact angular
integration produces a radial bilinear

$$
V''_{IJ}=\int_I dr\,
\left[
A_{IJ}(r)\overline{f'}g'
+B_{IJ}(r)\overline{f'}g
+\overline{B_{JI}(r)}\overline f g'
+C_{IJ}(r)\overline f g
\right]+B_{\partial I}[f,g].
$$

Use the shared real basis $\sqrt2\operatorname{Re}F$ and
$\sqrt2\operatorname{Im}F$. The Hermitian block becomes
$\left(\begin{smallmatrix}\operatorname{Re}H&\operatorname{Im}H\\
-\operatorname{Im}H&\operatorname{Re}H\end{smallmatrix}\right)$.
It is two identical copies only when the imaginary part vanishes; imaginary
mixed-species entries couple the two real components. This is the conversion
expected by the canonical assembler.

## 5. Analytic controls before a launch

The following are definition-level checks to implement locally before any
operator is passed downstream.

1. **Symmetry.** Check Hermiticity of the weak bulk form and separately verify
   that the strong bulk adjoint difference equals the two Green-flux difference.
2. **Continuum pure-gradient zero.** Insert $h=\mathcal L_\xi\sigma$ and
   compactly supported retained $k$. The strong linearized Einstein tensor is
   zero; the weak radial integral is zero after integration by parts. Its local
   weak density need not vanish coefficient by coefficient before that step.
3. **Exact retained shift kernels.** For the archived oblate kernel sections,
   verify $V''J_{\ker}=0$ in every exactly retained column. This is the same
   continuum identity as item 2, but pulled through the coefficient basis and
   the actual $J$ map.
4. **Cartesian TT normalization.** Independently use the flat Cartesian
   tensor formula and check that a transverse-traceless test gives
   $V''=\frac12\int|\nabla h^{TT}|^2/(4\pi)$.
5. **Finite-deformation convergence.** Independently evaluate the same weak
   form by finite differences of $-\int\sqrt\sigma R$ for a smooth
   compactly supported coefficient perturbation at fixed small $\epsilon$.
6. **Angular and radial convergence.** Double angular quadrature and radial
   Galerkin resolution separately. Angular projection products are not
   truncated before integration, and radial discretization is an additional
   regulator.

The boundary/domain handoff must still name the radial Sobolev space, the
allowed traces at $\partial I$, the treatment of $\Theta^r$, adjoints of
$B_{IJ}\partial_r$ terms, and the compatibility of these choices with shift
elimination, lapse rows, and physical evolution. A launchable stabilized
quadratic experiment also still needs the complete $K,E,D$ assembly and a
canonical momentum normalization.
