# Cubic and quartic interactions of the reduction-first Gaunt candidate

## Status and domain

This derives the complete homogeneous classical functionals $H_3$ and $H_4$
for the [declared reduced canonical candidate](reduced-canonical-gaunt.md).
It retains that construction's local existence assumptions, asymptotically
flat boundary conditions, maximal slicing and dual-Gram momentum map.
It is a change of regulator, not an equivalence theorem for the historical
Paper I or Paper II actions. No repeated binary angular projection is used.

The verification evaluates the functionals on the ten selected
$j=2$ profiles, reconstructs their complete cubic and quartic
configuration potentials by rotational invariants, and checks the local
geometry with an independent finite matrix exponential. The test protocol
and source were recorded before launch. The run passed its specified checks;
the reviewed coefficients and numerical limits are recorded below.

These are interaction vertices, not a one-loop scattering amplitude.
Fifth- and sixth-order vertices, quantization and a physical scattering
limit remain necessary for the requested finite-rank loop calculation.

## 1. Notation and the exact energy

Set $16\pi G=1$ and packet width $a=1$. For a finite flat-TT basis let

$$
Q=\sum_\alpha q^\alpha e_\alpha,\qquad
S=\sum_\alpha(\mathsf G_0^{-1}p)^\alpha e_\alpha,\qquad
s=\operatorname{tr}(S^2),\qquad
\mathcal C=\frac18(-\Delta)^{-1}.
$$

$Q$ and $S$ are first order in weak fields; the scalar $s$ is second order.
The coefficients of $S$ must not be mistaken for the canonical momenta.
All integrals and Green functions below are over $\mathbb R^3$ with
decaying conditions. The free basis cutoffs $N,K_{\rm rad}$ do not cap the
dependent geometry.

Write $\bar\gamma=e^Q$, $\bar R=\sum_{n\geq2}R_n$ and
$\psi=1+\sum_{n\geq2}u_n$. Flat transversality and tracelessness give
$R_1=0$ pointwise. The exact conformal equation is

$$
\bar\Delta\psi-\frac18\bar R\psi+\frac18\sigma^2\psi^{-7}=0.
$$

The geometric decomposition and the longitudinal elliptic operator are
the usual York--Lichnerowicz ingredients; see
[Cook, PDF pp. 8–9, Eqs. (20)–(32)](https://www.alphaxiv.org/pdf/gr-qc/0007085?page=8).
The particular canonical finite basis and the perturbative formulas below
are derived here rather than imported as a pre-existing finite-rank theory.

Since $\det\bar\gamma=1$ and $\bar\gamma-\delta$ is rapidly decaying,
integration of the conformal equation gives the ADM Hamiltonian

$$
H=-8\lim_{r\to\infty}\int_{S_r}\partial_r\psi\,dA
=\int d^3x\,[\sigma^2\psi^{-7}-\bar R\psi].
$$

No conformal-metric surface energy survives for this basis. Orthogonality
and the dual-Gram definition also give the exact identity

$$
\int\sigma^2=p^{\mathsf T}\mathsf G(q)^{-1}p.
$$

## 2. Momentum-constraint response and the Gram matrix

For a fixed direction $S$, the unprojected tangent to the metric is
$E_S=D\exp_Q[S]$. Its left-trivialized expansion is

$$
e^{-Q}E_S
=S-\frac12[Q,S]+\frac16[Q,[Q,S]]+O(Q^3S).
$$

The raw metric norm therefore satisfies

$$
\int\operatorname{tr}[(e^{-Q}E_S)^2]
=\int\operatorname{tr}S^2
-\frac1{12}\int\operatorname{tr}[Q,S]^2+O(Q^3S^2).
$$

The linear term vanishes by cyclicity. The York correction begins at
order $QS$, so it cannot restore a linear term in the Gram matrix.
Thus $\mathsf G_1=0$ in these particular coordinates. This is a coordinate
property, not the absence of gravitational cubic interactions.

To first nonzero order, set

$$
v_i[Q,S]=\frac12\left(
S_{lj}\partial_lQ_{ji}
-S_{jl}\partial_iQ_{jl}
-Q_{jk}\partial_kS_{ij}\right),
$$

$$
\Delta_L W_i=v_i,\qquad
\Delta_L W_i=\Delta W_i+\frac13\partial_i\partial_jW_j,\qquad
(LW)_{ij}=\partial_iW_j+\partial_jW_i-\frac23\delta_{ij}\partial_kW_k.
$$

Here $v_i$ is the first variation of $\bar D^j E_{S,ij}$.
For example, direct expansion gives

$$
v_i=\partial_j\!\left(\frac{QS+SQ}{2}\right)_{ij}
-Q_{jk}\partial_kS_{ij}
-\Gamma_{1\,ji}^{l}S_{lj}-\Gamma_{1\,jj}^{l}S_{il},
$$

which reduces to the displayed expression using the two flat-TT
conditions. Orthogonal subtraction of $LW$ yields

$$
\mathsf G_2[S,S]
=-\frac1{12}\int\operatorname{tr}[Q,S]^2-\|LW[Q,S]\|_2^2.
$$

Expanding the inverse Gram matrix gives the quartic kinetic correction

$$
H_{4,\mathrm{Gram}}=
\frac1{12}\int\operatorname{tr}[Q,S]^2+\|LW[Q,S]\|_2^2.
$$

The local commutator contribution is nonpositive for real symmetric
$Q,S$. The longitudinal norm is nonnegative. Dropping either term changes
the regulated Hamiltonian.

With Fourier convention $e^{-i\mathbf k\cdot\mathbf x}$,

$$
W_i(\mathbf k)=
-\frac1{k^2}\left(\delta_{ij}-\frac{k_i k_j}{4k^2}\right)v_j(\mathbf k),
$$

$$
\|LW\|_2^2
=2\int\frac{d^3k}{(2\pi)^3}
\left(\frac{|v|^2}{k^2}-\frac{|k\cdot v|^2}{4k^4}\right).
$$

For $Q=S=T$ there is a useful independent normalization check:
$v=-\nabla(T:T)/4$, $W=(3/2)\nabla\mathcal C(T:T)$, and

$$
\|LW[T,T]\|_2^2=\frac3{32}\int(T:T)^2.
$$

## 3. Complete third and fourth orders

The scalar constraint at second order gives

$$
u_2=\mathcal C(s-R_2).
$$

Expanding the exact energy and inverse Gram matrix now gives

$$
\boxed{H_2=\int s-\int R_2,\qquad H_3=-\int R_3,}
$$

$$
\boxed{
\begin{aligned}
H_4={}&-\int R_4
+\frac1{12}\int\operatorname{tr}[Q,S]^2
+\|LW[Q,S]\|_2^2\\
&-\left\langle7s+R_2,\mathcal C(s-R_2)\right\rangle .
\end{aligned}}
$$

Self-adjointness of $\mathcal C$ separates the three sectors:

$$
\begin{aligned}
H_4^{(q^4)}&=-\int R_4+\langle R_2,\mathcal C R_2\rangle,\\
H_4^{(q^2p^2)}&=\frac1{12}\int\operatorname{tr}[Q,S]^2
+\|LW\|_2^2+6\langle s,\mathcal C R_2\rangle,\\
H_4^{(p^4)}&=-7\langle s,\mathcal C s\rangle.
\end{aligned}
$$

The last term reduces to the previously checked momentum-direction
coefficient when $q=0$. The functional also fixes every mixed
configuration/momentum coefficient for an arbitrary retained TT basis;
it is not restricted to that direction. It is a classical Taylor
expansion on the domain of the canonical construction.

## 4. Curvature without losing surface derivatives

Let $\Gamma_n$ be the order-$n$ coefficient of the Levi-Civita connection
of $e^Q$. It follows directly from

$$
(e^Q)_n=\frac{Q^n}{n!},\qquad
(e^{-Q})_n=\frac{(-Q)^n}{n!}.
$$

Since $\Gamma^k_{ik}=0$, the integrated Ricci scalar equals the integral
of $\mathcal L=\bar\gamma^{ij}\Gamma^k_{il}\Gamma^l_{jk}$.
The required integrated densities are

$$
\begin{aligned}
\mathcal L_2&=\delta^{ij}\Gamma_1{}^k_{il}\Gamma_1{}^l_{jk},\\
\mathcal L_3&=2\delta^{ij}\Gamma_1{}^k_{il}\Gamma_2{}^l_{jk}
-Q^{ij}\Gamma_1{}^k_{il}\Gamma_1{}^l_{jk},\\
\mathcal L_4&=2\delta^{ij}\Gamma_1{}^k_{il}\Gamma_3{}^l_{jk}
+\delta^{ij}\Gamma_2{}^k_{il}\Gamma_2{}^l_{jk}\\
&\quad-2Q^{ij}\Gamma_1{}^k_{il}\Gamma_2{}^l_{jk}
+\frac12(Q^2)^{ij}\Gamma_1{}^k_{il}\Gamma_1{}^l_{jk}.
\end{aligned}
$$

Thus $\int R_n=\int\mathcal L_n$. However, the scalar Green function in
$H_4$ needs the **pointwise** $R_2$:

$$
R_2=\delta^{ij}\partial_k\Gamma_2{}^k_{ij}
-Q^{ij}\partial_k\Gamma_1{}^k_{ij}
-\delta^{ij}\Gamma_1{}^k_{jl}\Gamma_1{}^l_{ik}.
$$

Replacing $R_2$ here by $\mathcal L_2$ would discard derivatives inside a
nonlocal contraction and generally give the wrong interaction.
The integrated identity
$\int R_2=-\frac14\int|\nabla Q|^2$ checks the quadratic sign.

## 5. Exact evaluation and independent numerical checks

The source
[reduced_gaunt_vertices.py](../../verification/quantum_bridge/reduced_gaunt_vertices.py)
constructs the complete polynomial-Gaussian tensors and the above
functionals with rational arithmetic. Products retain all angular
components. Coordinate-space Gaussian moments evaluate the local terms.
The Fourier transform of each polynomial times $e^{-r^2}$ evaluates the
scalar and longitudinal Green contractions exactly. Both real and
imaginary Fourier components are retained for mixed parities.

For a monomial of degree $2n=2a+2b+2c$, the needed radial/angular moments
are

$$
\begin{aligned}
\int\frac{k_x^{2a}k_y^{2b}k_z^{2c}}{k^2}e^{-k^2/2}\,d^3k
&=4\pi\sqrt{\frac\pi2}
\frac{(2a-1)!!(2b-1)!!(2c-1)!!}{2n+1},\\
\int\frac{k_x^{2a}k_y^{2b}k_z^{2c}}{k^4}e^{-k^2/2}\,d^3k
&=4\pi\sqrt{\frac\pi2}
\frac{(2a-1)!!(2b-1)!!(2c-1)!!}{(2n+1)(2n-1)}.
\end{aligned}
$$

The second expression is used only for $n\geq1$, as ensured by the
$k\cdot v$ numerator. Odd monomials integrate to zero.

Two STF label matrices $E,B$ specify the electric and magnetic
configuration multiplets in the chosen ten-mode example. Rotational
invariance and spatial parity restrict its full cubic potential to
$\operatorname{tr}E^3$ and $\operatorname{tr}(EB^2)$.
The quartic potential has the five-invariant basis

$$
(\operatorname{tr}E^2)^2,\quad
(\operatorname{tr}B^2)^2,\quad
\operatorname{tr}E^2\operatorname{tr}B^2,\quad
(\operatorname{tr}EB)^2,\quad
\operatorname{tr}(E^2B^2).
$$

These exhaust the even-parity invariants: the symmetric square of spin
two has spins $0,2,4$, giving three mixed $E^2B^2$ scalars.
Five dedicated directions determine the quartic coefficients. Two
additional mixed and rationally rotated directions test the result.
The mixed momentum vertices are evaluated by the full functional, not
fitted from a claimed one-parameter equivalence.

The separate source
[reduced_finite_metric.py](../../verification/quantum_bridge/reduced_finite_metric.py)
evaluates $e^{hQ}$ with 70-digit arithmetic, propagating spatial and
directional derivatives through the full matrix-exponential series.
An independent matrix inversion and full Ricci contraction precede
centered differences at $h=0.0002,0.0001,0.00005$.
Richardson extrapolation checks $\mathcal L_2,\mathcal L_3,\mathcal L_4$,
pointwise $R_2$, both first Gram coefficients and the three York-source
components. Four spatial/polarization cases use 32 finite-metric
evaluations in total. Scaled errors divide by $1+|\text{target}|$.
The declared tolerance is $10^{-10}$.

These are finite-metric checks of weak-field ingredients. They are not a
second nonlinear solver for the entire finite-$q,p$ constraint problem.
Exact normalization controls, the aligned York identity, the previously
independent pure-momentum coefficient and the held-out rotations provide
additional checks.

### Reviewed coefficients and independent check

Run 63e0ae36-ddf8-459e-88b1-4bc29a49fcf2 used source
fb41c4bb0b423696ea12eb90cf75afd6214bc64c. The calculation took
71.4825 seconds and passed the seven exact-direction checks, the invariant
reconstruction, the rational rotation check and the four independent
finite-metric checks. The full output is preserved in
[reduced-vertices-results.json](data/reduced-vertices-results.json);
the durable artifact bundle also contains the raw log and frozen source.

The exact cubic potential for the two retained STF multiplets is

$$
H_3=\left(\frac{2\pi}{3}\right)^{3/2}
\left[\frac{446}{3}\operatorname{tr}E^3+
\frac{25004}{9}\operatorname{tr}(EB^2)\right].
$$

Thus the coordinate cancellation of the cubic kinetic term does not remove
the cubic gravitational potential. In the five-invariant quartic basis
listed above,

$$
H_4^{(q^4)}=\frac{\pi^{3/2}}{\sqrt2}
\sum_{a=1}^{5}c_a I_a(E,B),
$$

with

$$
(c_1,c_2,c_3,c_4,c_5)=
\left(
\frac{2978145633}{33554432},
\frac{69776871602481}{18790481920},
\frac{3047052105933}{2348810240},
\frac{1364019815211}{1174405120},
-\frac{40507834761}{146800640}
\right).
$$

The two held-out directions reproduce these coefficients exactly.
The simultaneous rational rotation of configuration and momentum labels
also leaves every separately evaluated Hamiltonian contribution unchanged.
The mixed and pure-momentum directional coefficients, and their separate
local/longitudinal/scalar contributions, are in the result JSON.

The independent finite-metric errors are:

| Configuration | Largest scaled coefficient error |
| --- | ---: |
| Electric axis | $3.14\times10^{-17}$ |
| Magnetic axis | $3.04\times10^{-14}$ |
| Noncommuting labels | $5.64\times10^{-15}$ |
| Held-out mixed labels | $9.83\times10^{-16}$ |

The denominator of each scaled error is $1+|\text{target}|$.
The largest error concerns the magnetic-axis quartic
$\Gamma\Gamma$ density. Metric determinants differ from unity by at most
$4\times10^{-69}$ at the stated precision.
These measurements check the expansion's local geometric ingredients and
selected exact contractions. They are not numerical error bounds for a
finite-$N$ scattering amplitude, which this run did not evaluate.

## 6. Consequences for the remaining quantum calculation

For a one-loop four-point functional, a background Hessian expansion
can involve a sixth-order tadpole, a fifth/cubic pair, a quartic/quartic
pair, quartic/cubic/cubic terms and four cubic vertices.
Finite-basis frequencies provide scales, so continuum dimensional
regularization's disappearance of scaleless tadpoles is not a license
to omit the first two structures. Higher vertices and any ordering or
measure terms must be derived from the same candidate.

The negative $p^4$ term is a local Taylor coefficient of the constraint
solution. Exponentiating the polynomial $H_2+H_3+H_4$ as a global quantum
Hamiltonian would not inherit a lower bound merely from the positive
classical ADM energy. The full domain, quantization and higher terms
matter for interacting Euclidean reconstruction.

The previously derived uniform-area matching and oscillator transfer
kernels remain quadratic results. These interaction vertices can supply
their next corrections, but do not determine the full local corner
symplectic reduction, its microscopic multiplicities or the simulation
measure. No physical large-$N$ loop limit has been evaluated here.
