# A reduction-first candidate connecting the constrained geometry and quantum variables

## Status and change of theory

This is an explicitly proposed **revision** of the regulator, motivated by
the [mixed-constraint obstruction](canonical-mixed-bracket.md). It is not an
equivalence proof for the historical finite-rank actions.

The independent variables are finitely many physical configuration and
momentum coefficients. The dependent metric, conformal factor and
longitudinal corrections are reconstructed without imposing the same angular
cap on them. Their products use complete pointwise/Gaunt contractions.
Consequently the reconstructed geometry can contain harmonics above the
free-data cutoff. Those harmonics are functions of the retained variables,
not additional independent oscillators.

We give a local canonical construction, an explicit family of smooth vacuum
data, its quadratic Hamiltonian, certain fourth- and sixth-order Hamiltonian
coefficients, and the induced **leading uniform-area sector**. The numerical
protocol is committed before its run; no numerical result is asserted here.

The setting is asymptotically flat gravity on $\mathbb R^3$, with maximal
slicing and an asymptotic time translation. This addresses the scattering
side of the project. A finite causal diamond with all its local corner
charges, null boundaries and gluing constraints is additional work.
The complete one-loop finite-$N$ amplitude remains uncomputed.

## 1. Canonical variables on the constraint surface

First set $16\pi G=1$ and choose the extrinsic-curvature sign for which
$\pi^{ij}=\sqrt\gamma(K^{ij}-\gamma^{ij}K)$. Let
$e_{\alpha\,ij}(\mathbf x)$ be a finite set of
linearly independent, rapidly decaying, real symmetric tensors which are
trace-free and transverse in the flat metric. Their total angular momenta
are at most $N-1$; a separate radial basis cutoff is understood.

Define a positive unimodular conformal metric by a matrix exponential,

$$
Q_{ij}(q,\mathbf x)=\sum_\alpha q^\alpha e_{\alpha\,ij}(\mathbf x),
\qquad \bar\gamma(q)=\exp Q.
$$

The exponential here is the ordinary pointwise matrix exponential. No
binary angular projection is inserted into it. Put
$E_{\alpha\,ij}=\partial\bar\gamma_{ij}/\partial q^\alpha$.
Because $\det\bar\gamma=1$, these tangent tensors are trace-free with
respect to $\bar\gamma$.

Let $\Pi_{\bar\gamma}^{\rm TT}$ be the orthogonal York projection of a
trace-free tensor onto the transverse-traceless subspace, with decaying
boundary conditions. Write

$$
T_{\alpha\,ij}=\Pi_{\bar\gamma}^{\rm TT}E_{\alpha\,ij},\qquad
\mathsf G_{\alpha\beta}(q)
=\int d^3x\sqrt{\bar\gamma}\,
T_{\alpha\,ij}T_\beta^{ij}.
$$

Indices in this section are raised with $\bar\gamma$. On a domain where
this Gram matrix is invertible, prescribe the conformal TT curvature

$$
\sigma^{ij}(q,p)
=\sum_{\alpha,\beta}
p_\alpha(\mathsf G^{-1})^{\alpha\beta}T_\beta^{ij}.
$$

Solve the maximal Lichnerowicz equation,

$$
\bar\Delta\psi-\frac18\bar R\psi
+\frac18\sigma_{ij}\sigma^{ij}\psi^{-7}=0,\qquad
\psi\longrightarrow1\quad (r\longrightarrow\infty),
$$

on its positive solution branch. Reconstruct

$$
\gamma_{ij}=\psi^4\bar\gamma_{ij},\qquad
K^{ij}=\psi^{-10}\sigma^{ij},\qquad K=0,\qquad
\pi^{ij}=\psi^{-4}\sqrt{\bar\gamma}\,\sigma^{ij}.
$$

The conformal powers and constraint equation follow the York--Lichnerowicz
construction; see [Cook, PDF pp. 7–9, Eqs. (16), (22), (32)](https://www.alphaxiv.org/pdf/gr-qc/0007085?page=7).
Existence depends on the conformal metric and boundary conditions. In
particular, maximal asymptotically Euclidean data have a Yamabe-condition
restriction; it must not be omitted. See
[Bartnik and Isenberg, PDF p. 17](https://www.alphaxiv.org/pdf/gr-qc/0405092?page=17).

### Local proposition, with explicit hypotheses

Assume that the decaying TT projection and the positive conformal solution
exist smoothly in the chosen coefficient domain, that no relevant decaying
gauge zero mode is present, and that $\mathsf G$ is positive definite.
Then the above map is a symplectic embedding into continuum vacuum
constraint data, with canonical coordinates $(q^\alpha,p_\alpha)$.

**Constraints.** The conformal divergence identity gives

$$
D_jK^{ij}=\psi^{-10}\bar D_j\sigma^{ij}=0.
$$

The scalar constraint is

$$
R[\gamma]-K_{ij}K^{ij}
=-8\psi^{-5}
\left(\bar\Delta\psi-\frac18\bar R\psi
+\frac18\sigma_{ij}\sigma^{ij}\psi^{-7}\right)=0.
$$

**Symplectic form.** The ADM one-form pulls back as

$$
\begin{aligned}
\Theta_{\rm ADM}
&=\int d^3x\,\pi^{ij}\delta\gamma_{ij}\\
&=\int d^3x\sqrt{\bar\gamma}\,\sigma^{ij}
  \delta\bar\gamma_{ij}
  +4\int d^3x\sqrt{\bar\gamma}\,
    \psi^{-1}\delta\psi\,\sigma^{ij}\bar\gamma_{ij}\\
&=\sum_{\alpha,\beta,\gamma}
p_\alpha(\mathsf G^{-1})^{\alpha\beta}
\langle T_\beta,E_\gamma\rangle\,dq^\gamma
=\sum_\alpha p_\alpha\,dq^\alpha.
\end{aligned}
$$

The second term vanishes by tracelessness. Orthogonality of the TT
projection gives $\langle T_\beta,E_\gamma\rangle=\mathsf G_{\beta\gamma}$.
Thus $\Omega=d\Theta=\sum_\alpha dp_\alpha\wedge dq^\alpha$ and
$\{q^\alpha,p_\beta\}=\delta^\alpha_\beta$ in the corresponding convention.

At $q=0$, the TT projection fixes $e_\alpha$ and the Gram matrix is their
positive $L^2$ matrix. Positivity therefore persists on any sufficiently
small domain where the projection is continuous. This observation does not
assert a global coordinate chart or prove uniform-in-$N$ bounds for the
nonlinear conformal problem.

Define the classical regulated Hamiltonian by the reconstructed ADM
energy,

$$
H_{N,K_{\rm rad}}(q,p)
=M_{\rm ADM}[\gamma(q,p),K(q,p)].
$$

Any trajectory generated by this finite canonical Hamiltonian reconstructs
constraint-satisfying data while it stays in the stated domain. In general
the embedded finite-dimensional manifold is **not invariant under the full
continuum Einstein flow**. Its Hamiltonian evolution is a reduced
Galerkin approximation; convergence of that evolution and of its quantum
loops must still be shown.

This prescription also makes a choice: TT momentum components orthogonal
to the span of the $T_\alpha$ are set to zero. They have not been
quantum-mechanically integrated out. The missing Wilsonian effects cannot
be inferred from the classical symplectic calculation.

## 2. An explicit finite free basis with both tensor parities

Measure length in a fixed packet width $a$, initially setting $a=1$.
For a constant symmetric trace-free Cartesian tensor $E_{ij}$ and
$f(\mathbf x)=e^{-r^2/2}$, define

$$
\mathcal P_{ij}=\delta_{ij}\Delta-\partial_i\partial_j,\qquad
T^E_{ij}
=\left(\mathcal P_i{}^k\mathcal P_j{}^l
-\frac12\mathcal P_{ij}\mathcal P^{kl}\right)E_{kl}f.
$$

Commutativity of the Cartesian derivatives proves symmetry, vanishing
trace and vanishing divergence. The operator is rotation-equivariant, so
the five components of $E$ produce a total-$j=2$ multiplet. Its symmetric
curl gives the opposite parity:

$$
T^B_{ij}
=\frac12\left(
\epsilon_{iak}\partial_aT^E_{kj}
+\epsilon_{jak}\partial_aT^E_{ik}\right).
$$

It is also TT and has total $j=2$. We retain five real polarizations in
each parity, giving ten selected modes at angular rank $N=3$.
This is one radial profile per parity, not the entire $j=2$ radial space.

Exact Gaussian moments give, for electric modes,

$$
\int T^{E_1}_{ij}T^{E_2}_{ij}\,d^3x
=\frac{189}{8}\pi^{3/2}\operatorname{tr}(E_1E_2).
$$

The magnetic Gram matrix is $11/2$ times this matrix. The two parity
sectors are orthogonal. For the chosen orthogonal STF polarizations the
gradient matrix is diagonal as well:

$$
\frac{\int|\nabla T^E|^2}{\int|T^E|^2}=\frac{11}{2},
\qquad
\frac{\int|\nabla T^B|^2}{\int|T^B|^2}=\frac{13}{2}.
$$

At quadratic order the TT Einstein Hamiltonian is therefore

$$
H_2=\sum_\alpha\left[
\frac{p_\alpha^2}{g_\alpha}
+\frac{g_\alpha\omega_\alpha^2}{4}(q^\alpha)^2\right],
\quad
\omega_E=\frac1a\sqrt{\frac{11}{2}},\qquad
\omega_B=\frac1a\sqrt{\frac{13}{2}},
$$

with five modes of each frequency. These frequencies characterize the
selected localized basis. They are not graviton masses or a proposed
modified dispersion relation.

The canonical change
$Q_\alpha=\sqrt{g_\alpha/2}\,q^\alpha$,
$P_\alpha=\sqrt{2/g_\alpha}\,p_\alpha$
puts $H_2$ into the usual positive oscillator form.

## 3. Nonlinear reconstruction and specific interaction coefficients

Take $q=0$ and the electric axisymmetric direction
$E=\operatorname{diag}(-1/2,-1/2,1)$, with
$\sigma_{ij}=\varepsilon T_{ij}$.
Then the momentum constraint vanishes exactly for every positive $\psi$,
and the remaining equation is

$$
\Delta\psi=-\frac{\varepsilon^2}{8}\mathcal S\psi^{-7},
\qquad \mathcal S=T_{ij}T_{ij},\qquad \psi\to1.
$$

With $x=\cos\theta$, complete Gaunt contraction yields

$$
\mathcal S(r,x)
=e^{-r^2}\left[S_0(r)P_0(x)+S_2(r)P_2(x)+S_4(r)P_4(x)\right],
$$

where

$$
\begin{aligned}
S_0&=54-72r^2+48r^4-\frac{48}{5}r^6+\frac35r^8,\\
S_2&=72r^2-\frac{444}{7}r^4+\frac{96}{7}r^6-\frac67r^8,\\
S_4&=\frac{108}{7}r^4-\frac{144}{35}r^6+\frac9{35}r^8.
\end{aligned}
$$

The nonzero $j=4$ source already shows why the dependent geometry cannot
be restricted to the free $j=2$ sector.

Define the positive Newton operator

$$
(\mathcal C h)(\mathbf x)
=\frac1{32\pi}\int_{\mathbb R^3}
\frac{h(\mathbf y)}{|\mathbf x-\mathbf y|}\,d^3y.
$$

The equation becomes $\psi=1+\varepsilon^2\mathcal C(\mathcal S\psi^{-7})$.
Let $u_2=\mathcal C\mathcal S$. A radial positive majorant formed from the
absolute polynomial coefficients above proves

$$
0\leq u_2(\mathbf x)\leq K_*=\frac{2301}{56},
\qquad u_2(0)=\frac{87}{40}.
$$

On $1\leq\psi\leq1+\varepsilon^2u_2$, the full integral map is a
contraction if $7\varepsilon^2K_*<1$. It maps that interval into itself,
so this proves a positive solution and convergence in the specified
small-amplitude domain. This bound applies to the **untruncated
integral equation**; the numerical radial/angular approximation is checked
separately.

Writing
$\psi=1+\varepsilon^2u_2+\varepsilon^4u_4+O(\varepsilon^6)$ gives
$u_4=-7\mathcal C(\mathcal S u_2)$.
The ADM energy, in $16\pi G=1$ units, is

$$
H(0,p(\varepsilon))
=\varepsilon^2\int \mathcal S\psi^{-7}\,d^3x.
$$

The relation of the asymptotic conformal factor to ADM energy is given in
[Cook, PDF p. 23, Eq. (75)](https://www.alphaxiv.org/pdf/gr-qc/0007085?page=23).
The physical energy restores the factor $a/(16\pi G)$ when the
dimensionless seed and amplitude above are used.

Self-adjointness of $\mathcal C$ gives the expansion

$$
\boxed{
H(0,p(\varepsilon))
=\varepsilon^2\mathcal I_1
-7\varepsilon^4\mathcal I_2
+77\varepsilon^6\mathcal I_3+O(\varepsilon^8),
}
$$

where

$$
\mathcal I_1=\int\mathcal S=\frac{567}{16}\pi^{3/2},\quad
\mathcal I_2=\int\mathcal S\,\mathcal C\mathcal S>0,\quad
\mathcal I_3=\int\mathcal S(\mathcal C\mathcal S)^2>0.
$$

For example,
$\int\mathcal S u_4=-7\int\mathcal S u_2^2$, which explains
$77=49+28$ without solving for $u_4$ numerically.
Along this canonical direction $p=\mathcal I_1\varepsilon$.
These are specific momentum-direction Taylor coefficients of the proposed
Hamiltonian, not its full multivariate vertices.

The code computes $\mathcal I_2$ exactly by Fourier transformation and
Gaussian moments, and independently through radial Green functions.
It computes $\mathcal I_3$ with two quadratures.
For a source harmonic $s_\ell(r)=e^{-r^2}S_\ell(r)$, the independent
linear solution is

$$
u_{2,\ell}(r)=\frac1{8(2\ell+1)}
\left[r^{-\ell-1}\int_0^r s^{\ell+2}s_\ell(s)\,ds
+r^\ell\int_r^\infty s^{1-\ell}s_\ell(s)\,ds\right].
$$

This is compared with a separate radial finite-difference Poisson solver
before that solver is used for the nonlinear problem.

## 4. Bulk-induced uniform corner spectrum at leading order

This section isolates the uniform area--boost pair of the earlier
[corner derivation](corner-quantization-and-transfer.md). It does not
quantize all local corner fields or the full normal-frame algebra.
Restore $G$ and $\hbar$ here. Restoring $G$ rescales the canonical
normalization of the oscillator coordinates; their frequencies are unchanged.

Consider the geometrically specified sphere $r=R\gg a$ in the chosen
asymptotic coordinate/observer prescription. The retained TT profiles
decay as Gaussians. At quadratic weak-field order the conformal metric
has negligible direct boundary contribution, and the monopole conformal
factor is $1+G H_2/(2R)$. Angular multipoles integrate to zero at this order.
Thus, up to Gaussian-tail corrections at that sphere,

$$
A_R=A_0+8\pi G R H_2+O(\Phi^3),\qquad A_0=4\pi R^2.
$$

Here $O(\Phi^3)$ denotes interactions of third and higher total field
order. This is a far-sphere weak-field relation, not an exact
mass--area relation for arbitrary regions or horizons.

Quantize the canonical oscillators. Fix the vacuum subtraction explicitly:
the normal-ordered quadratic ADM energy of the Minkowski vacuum is zero,
and its reference corner area is $A_0$. Then

$$
\widehat H_2=\hbar\omega_E\widehat n_E+
\hbar\omega_B\widehat n_B,\qquad
\widehat A_R^{(2)}=A_0+8\pi G R\widehat H_2.
$$

The sums $\widehat n_E$ and $\widehat n_B$ each run over five independent
oscillators. Their leading area eigenvalues and multiplicities are

$$
A_{n_E,n_B}^{(2)}
=A_0+\frac{8\pi G\hbar R}{a}
\left(\sqrt{\frac{11}{2}}n_E+\sqrt{\frac{13}{2}}n_B\right),
\quad
d_{n_E,n_B}=\binom{n_E+4}{4}\binom{n_B+4}{4}.
$$

The irrational ratio of the two frequencies makes different pairs
$(n_E,n_B)$ distinct energy levels. The multiplicity follows from the
derived oscillator content, not an assigned two- or three-state corner.

### Matching to the extended area--boost pair

Let $a_c=A/(8\pi G)>0$ and $\eta$ be the uniform boost angle, with
$\{\eta,a_c\}=1$. Extend the bulk phase space by this pair and impose the
geometric matching constraint

$$
C=a_c-\frac{A_R(q,p)}{8\pi G}=0.
$$

Its gauge action shifts $\eta$ and acts on the bulk variables through
their area generator. The section $\eta=0$ pulls the extended
symplectic form back to the original canonical bulk form.

At the stated quadratic order the kinematical Hilbert space is
$\mathcal F_{10}\otimes L^2(\mathbb R_+,da_c)$ and

$$
\widehat C
=a_c-\frac{A_0}{8\pi G}-R\widehat H_2.
$$

For test states with finite Fock support and smooth area wavefunctions,
group averaging produces

$$
\langle\Phi,\Psi\rangle_{\rm phys}
=\sum_{\mathbf n}
\overline{\Phi_{\mathbf n}(a_{\mathbf n})}
\Psi_{\mathbf n}(a_{\mathbf n}),\qquad
a_{\mathbf n}=\frac{A_0}{8\pi G}+R E_{\mathbf n}.
$$

This follows by integrating
$\exp(i\lambda\widehat C/\hbar)$ with measure
$d\lambda/(2\pi\hbar)$, which produces
$\delta(a_c-a_{\mathbf n})$. The resulting physical states are
distributional in the extended area coordinate; after the null quotient
their positive Hilbert space is the bulk Fock space.
The independent half-line area degree of freedom has been removed by
matching. The earlier infinite-area-window result for an *unmatched*
area--boost pair therefore does not imply an infinite count in this
different, reduced sector.

For a chosen area window $A_R^{(2)}\leq A_{\max}$ define

$$
B=\frac{(A_{\max}-A_0)a}{8\pi G\hbar R}.
$$

The leading sector dimension is the finite, calculable number

$$
D(B)=
\sum_{\sqrt{11/2}\,n_E+\sqrt{13/2}\,n_B\leq B}
\binom{n_E+4}{4}\binom{n_B+4}{4}.
$$

The code evaluates this with exact radical inequalities and independently
by recursion over the ten individual modes.

**Limits of this count.** The area window is chosen, not derived as a
fundamental law. The relation to $H_2$ is a weak-field, far-sphere
approximation with a stated vacuum subtraction. The full field theory
has further radial/angular modes and local corner charges. At fixed ten
modes this count grows polynomially with $B$; it is not a derivation of
black-hole entropy, the microscopic degeneracies in Paper I, or a
nontrivial finite-dimensional representation of the undeformed normal
$SL(2,\mathbb R)$ algebra.
There is no uniform claim for arbitrarily high occupation or increasing
radial bases: the fields must remain weak and localized inside the chosen
sphere for the stated area approximation to apply.

## 5. Euclidean and Lorentzian evolution of the same free sector

The positive, self-adjoint number Hamiltonian above gives

$$
T_\tau=e^{-\tau\widehat H_2/\hbar},\qquad
U_t=e^{-it\widehat H_2/\hbar}.
$$

Both act on the same group-averaged states. The uniform corner boost
Hamiltonian, after subtracting its background charge, is
$H_c-\nu a_0=\nu R H_2$ at this order. Thus its evolution differs by the
declared clock factor $\nu R$; this is not an identification of a general
local boost with an exact global time translation.
The free Gibbs trace is finite for every $\tau>0$:
$\operatorname{Tr}T_\tau=
(1-e^{-\omega_E\tau})^{-5}(1-e^{-\omega_B\tau})^{-5}$.

For the Hamiltonian with the usual oscillator zero-point term, each
coordinate-space Euclidean factor is

$$
K_\tau(Q',Q)=
\sqrt{\frac{\omega}{2\pi\hbar\sinh(\omega\tau)}}
\exp\left[-\frac{\omega}
{2\hbar\sinh(\omega\tau)}
\left((Q'^2+Q^2)\cosh(\omega\tau)-2Q'Q\right)\right].
$$

The normal-ordered kernel multiplies the product of these factors by
$\exp(\tau\sum_\alpha\omega_\alpha/2)$.
Gaussian integration verifies sewing; the number basis verifies
self-adjointness, positivity and the semigroup property directly.
Analytic continuation gives the unitary $U_t$ with the usual oscillator
kernel boundary values.

The vacuum covariance is
$\langle Q_\alpha(\tau)Q_\beta(0)\rangle
=\delta_{\alpha\beta}\hbar e^{-\omega_\alpha|\tau|}/(2\omega_\alpha)$.
Its reflected positive-time matrix is a Gram matrix. Higher Gaussian
correlators obey the same Fock-space reflection-positivity construction.

The leading area projector commutes with $H_2$. Restricting to its finite
sector therefore preserves the operator semigroup and unitarity. A
projected coordinate kernel need not be entrywise positive; these
operator statements do not establish a positive sampling weight for the
historical Euclidean simulations.

Interactions change both the Hamiltonian and the area observable.
Their operator ordering, counterterms, domain and reconstruction remain
unproved. This section establishes the shared **quadratic** sector only.

## 6. A controlled free limit and its precise assumptions

There is a direct convergence criterion for the free candidate. Let
$V_{N,K_{\rm rad}}$ be nested TT subspaces whose union is a form core for
$-\Delta$ on the physical TT $L^2$ space. The radial basis cutoff
$K_{\rm rad}$ is distinct from the numerical radial grid used below.
All physical lengths, packet widths and smearings are held fixed.

For $z>0$, let $u=(-\Delta+z)^{-1}f$ and let $u_{N,K_{\rm rad}}$ solve the
Galerkin variational problem on $V_{N,K_{\rm rad}}$. Orthogonality gives

$$
\|u-u_{N,K_{\rm rad}}\|_z
\leq\inf_{v\in V_{N,K_{\rm rad}}}\|u-v\|_z,\qquad
\|w\|_z^2=\|\nabla w\|_2^2+z\|w\|_2^2.
$$

The right-hand side tends to zero under the stated form-core condition.
This proves strong resolvent convergence. With a fixed positive infrared
mass during the comparison, the free covariance is a bounded continuous
function of the spatial operator, so its fixed-smearing Euclidean limit
follows as well. The massless infrared limit is a separate step.

Polynomial-Gaussian tensor seeds followed by the TT differential operator
give one possible route to a form core: smooth TT Fourier data supported
away from $k=0$ have a smooth preimage after dividing by $k^4$, and
Hermite expansions approximate that preimage with its required derivatives.
Density then removes the temporary excluded neighborhood of $k=0$.
The ten modes tested here are only the first selected block, not a
demonstration that all numerical implementations realize this core.

This free variational bound does not control interacting loop sums,
renormalization or a fundamental rank--area cap. Treating $N$ as a
removable regulator here is an additional interpretation to be reconciled
with any retained microscopic area law.

## 7. Tracked verification protocol

The new experiment selects the reduced-data calculation through the fixed
entry point. It uses Python's standard library on local CPU.

- Derive all ten Cartesian TT modes by exact rational differentiation;
  verify their traces, divergences, full Gram and gradient matrices,
  positive frequencies, and dual-Gram canonical pairing.
- Derive the three source Gaunt harmonics from the Cartesian tensor
  square, its norm and the exact Fourier Coulomb energy.
- Compare the quartic energy with two independent Green-function
  quadratures; evaluate the sixth-order coefficient with both.
- Count six leading area windows by the level multiplicities and by an
  independent ten-mode recursion.
- Compare the analytic linear solution with a separate Poisson solver
  on 256, 512 and 1024 radial intervals, with $R=10a$.
- Solve thirteen nonlinear cases: three amplitudes
  $\rho=\varepsilon^2K_*\in\{0.005,0.02,0.08\}$ on those three grids,
  auxiliary angular cutoffs $J=2,4,8,12$ at the largest amplitude, and
  an independent 32-to-48-node angular quadrature comparison.

The full equation is evaluated on separate angular validation points.
A fourth-order radial derivative checks the solutions produced by the
second-order solver. The tests distinguish angular projection error,
radial discretization error, outer-tail bounds and fixed-point error.
The $J=2$ calculation is a deliberate failure control: the dependent
conformal factor is then unable to carry its required $j=4$ response.

The source is
[reduced_gaunt_exact.py](../../verification/quantum_bridge/reduced_gaunt_exact.py)
and
[reduced_gaunt_completion.py](../../verification/quantum_bridge/reduced_gaunt_completion.py).
Results, source commit and raw-log provenance will be recorded after the
tracked run is reviewed.

## 8. What still prevents closing the full goal

The construction provides a specified candidate and several connected
ingredients. It does not supply the complete nonlinear Hamiltonian
through the order needed for a four-point loop, its finite-$N$ physical
amplitude and renormalization, or a controlled interacting limit.
The spectrum and transfer kernels above apply to its leading uniform
boundary sector, not all causal-diamond corner degrees of freedom.
Full local boundary reduction, microscopic state counting, interacting
operator reconstruction, and the observable map to the simulations remain
open.

The [continuum quantum reference](continuum-one-loop-reference.md) is now
available to test any completed amplitude. Passing the present constraint
and free-sector checks would make that later test well motivated; it
would not make the three principal quantum deliverables complete.
