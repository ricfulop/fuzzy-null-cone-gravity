# A constructive positive quantum completion of the local reduced Hamiltonian

## Result and scope

**Result: a constructive finite-mode existence completion is available.**
A smooth local Hamiltonian with a nondegenerate positive quadratic minimum
admits a chosen global extension whose direct Weyl quantization is
self-adjoint, bounded below, has compact resolvent and a ground state.
Subtracting its ground energy gives interacting positive Euclidean transfer
and unitary Lorentzian evolution on the same Hilbert space. This direct
Weyl route preserves the existing quantization prescription exactly.
An alternative sum-of-squares quantization is nonnegative before shifting,
but introduces a specified higher-order quantum correction.

The extension agrees with the exact classical Hamiltonian on a neighborhood
of the origin. Both routes preserve the **formal bare midpoint/Weyl
one-loop vertices** there. The sum-of-squares alternative generally differs
from direct Weyl quantization at order $\hbar^2$. Exact finite-$\hbar$
observables depend on the chosen global extension for either route.

This is expressly an existence completion, not a derivation of the
original global or microscopic theory. Its local input is the Hamiltonian
of [the reduction-first canonical Gaunt candidate](reduced-canonical-gaunt.md),
with the cubic/quartic functionals in [reduced-vertices.md](reduced-vertices.md)
and the higher-order and midpoint/Weyl prescriptions in
[reduced-sixth-order-and-loop-kernel.md](reduced-sixth-order-and-loop-kernel.md).
The [complete quartic and two-point report](complete-quartic-and-two-point.md)
explicitly defines its loop result as a coefficient of $\hbar$, with no
order-$\hbar$ counterterm. That distinction is essential here.

All proofs below are new analytic deductions for this completion. No
experimental run, numerical spectrum, or computed scattering amplitude is
claimed. The separate [local existence and ADM smoothness proof](local-existence-and-adm-smoothness.md)
supplies the classical input for each fixed finite Schwartz TT basis,
using the stated weighted elliptic theorems; finite Taylor coefficients
alone would not establish it. No global
Einstein solution, uniform-in-cutoff neighborhood, or gravitational
boundary representation is assumed to follow from this argument.

## 1. Precise finite-mode theorem

Let $z=(q,p)\in\mathbb R^{2d}$ be fixed canonical coordinates, with
$d<\infty$, and choose fixed units so that the Euclidean norm of $z$
is defined. Let $U$ be an open neighborhood of zero and suppose

$$
H\in C^\infty(U;\mathbb R),\qquad H(0)=0,\qquad dH(0)=0,
\qquad H''(0)>0.                                                \tag{1}
$$

For every fixed $\hbar>0$ the construction below supplies:

1. A smooth $H_{\rm ext}$ on all of $\mathbb R^{2d}$, equal to $H$
   near zero and to its quadratic Taylor polynomial outside a compact set,
   with $c|z|^2\leq H_{\rm ext}(z)\leq C|z|^2$ and a unique global minimum.
2. Real symbols $a_i$, $i=1,\ldots,2d$, linear plus smooth compactly
   supported functions, such that $H_{\rm ext}=\sum_i a_i^2$.
3. Two explicit self-adjoint quantizations on $L^2(\mathbb R^d,dq)$:
   the preferred direct Weyl operator
   $\widehat H_\hbar^W=\operatorname{Op}_\hbar^W(H_{\rm ext})$,
   which is bounded below, and the nonnegative sum-of-squares operator
   from the closed form $\sum_i\|\operatorname{Op}_\hbar^W(a_i)\psi\|^2$.
4. Compact resolvent, discrete spectrum of finite multiplicity tending
   to infinity, a normalizable ground state, and trace-class heat evolution
   for every positive Euclidean time.
5. Exact Weyl symbol $H_{\rm ext}$ for the preferred route.
   The alternative has symbol $H_{\rm ext}+\hbar^2r_2+O(\hbar^4)$,
   with the explicit $r_2$ in Section 5 and no order-$\hbar$ term.
6. A positive Euclidean transfer semigroup, a vacuum reflection-positive
   Gram form, and unitary Lorentzian evolution, all from the same operator.

A unique **quantum** ground state is not asserted. The ground eigenspace
has finite multiplicity, and any normalized vector in it can be used for
the vacuum Gram construction. The theorem proves a finite number of
canonical modes, not a finite-dimensional Hilbert space.

For the normalized ten-mode candidate,

$$
H_2=\frac12\sum_{\alpha=1}^{10}
(P_\alpha^2+\omega_\alpha^2Q_\alpha^2),\qquad
\omega_E^2=\frac{11}{2},\quad\omega_B^2=\frac{13}{2},             \tag{2}
$$

so the positive-Hessian hypothesis follows directly once the exact local
smooth Hamiltonian and its stated quadratic expansion are available.
These are basis frequencies, not masses or a continuum dispersion law.

## 2. Exact classical factorization and a global extension

Choose a ball $B_r\Subset U$ on which $H''(z)$ is uniformly positive.
The integral form of Taylor's theorem gives

$$
A(z)=\int_0^1(1-t)H''(tz)\,dt,\qquad
H(z)=z^TA(z)z,\qquad A_0:=A(0)=\frac12H''(0).                  \tag{3}
$$

For example, if $mI\leq H''\leq MI$ on this ball, then
$mI/2\leq A\leq MI/2$. The positive symmetric square-root map is smooth
on positive definite matrices. Hence

$$
S(z)=A(z)^{1/2},\qquad S_0=A_0^{1/2}
$$

are smooth, symmetric and uniformly positive on a possibly smaller ball.
Choose $\chi\in C_c^\infty(B_r)$ with $0\leq\chi\leq1$ and
$\chi=1$ on a neighborhood of zero. Define globally

$$
B(z)=S_0+\chi(z)(S(z)-S_0),\qquad a(z)=B(z)z.                   \tag{4}
$$

The compact term in (4) is extended by zero outside $B_r$; there is no
need to evaluate $A$ where it has not been defined. Since $B$ is a convex
combination of uniformly positive symmetric matrices,

$$
b_-I\leq B(z)\leq b_+I\quad\hbox{for all }z,
\qquad 0<b_-\leq b_+<\infty.                                  \tag{5}
$$

Consequently

$$
\boxed{
H_{\rm ext}(z)=|a(z)|^2=z^TB(z)^2z,
\qquad b_-^2|z|^2\leq H_{\rm ext}(z)\leq b_+^2|z|^2.
}                                                               \tag{6}
$$

Where $\chi=1$, this is exactly $H$, not a finite Taylor approximation.
Outside the support of $\chi$, it is exactly $H_2=z^TA_0z$.
Its only zero is $z=0$, proving uniqueness of its global minimum.
Moreover

$$
a_i(z)=l_i(z)+k_i(z),\qquad l_i(z)=(S_0z)_i,\qquad
k_i(z)=[\chi(z)(S(z)-S_0)z]_i\in C_c^\infty.                    \tag{7}
$$

For an arbitrary cutoff, (6) alone does not exclude other stationary
points or higher local minima in the transition region. That stronger
property is not needed. It can also be arranged: choose a fixed scaled
cutoff $\chi_r(z)=\chi_*(z/r)$ and shrink $r$. Since
$S(z)-S_0=O(|z|)$, differentiation of (4) shows

$$
\|D^2(H_{{\rm ext},r}-H_2)\|_{L^\infty}=O(r).
$$

Indeed the first two derivatives of
$\chi_r(S-S_0)$ are respectively $O(1)$ and $O(r^{-1})$, while its
size is $O(r)$; the extra factors of $z$ in (6) give the displayed bound.
For sufficiently small $r$, $H_{{\rm ext},r}$ is globally strictly convex.
Neither version proves uniqueness of the quantum ground state.

The map $z\mapsto a(z)$ is **not** asserted to be canonical. No change
of quantum coordinates or phase-space integration measure is performed.
The factors are symbols quantized in the original canonical $z$.
Thus writing a sum of classical squares does not turn the interacting
system into free oscillators.

## 3. Operator domains and closure of the positive form

In Sections 3--5, $\widehat H_\hbar$ denotes the sum-of-squares
alternative. Section 6 then constructs the preferred direct Weyl
operator, with the same domain and spectral properties but a different
higher-order symbol.

Fix the Weyl convention

$$
(\operatorname{Op}_\hbar^W(f)\psi)(x)
=\frac1{(2\pi\hbar)^d}\int
 e^{i(x-y)\cdot p/\hbar}f\!\left(\frac{x+y}{2},p\right)
 \psi(y)\,dy\,dp.                                               \tag{8}
$$

Thus $\widehat q_j=x_j$ and $\widehat p_j=-i\hbar\partial_j$.
Let $L_i=\operatorname{Op}_\hbar^W(l_i)$ and
$K_i=\operatorname{Op}_\hbar^W(k_i)$.

Every $L_i$ is the self-adjoint closure of a real nonzero linear combination
of $\widehat q,\widehat p$ on Schwartz space. This can be seen directly
from its strongly continuous unitary group, a translation followed by a
linear phase, or by a linear canonical transformation to one coordinate.
Because $k_i$ is real, $K_i$ has a Hermitian kernel. Its kernel is a
Schwartz function of $(x,y)$: the partial Fourier transform in $p$ is
Schwartz in $(\frac{x+y}{2},\frac{x-y}{\hbar})$, and the change of
variables is invertible for fixed $\hbar>0$. In particular $K_i$ is
bounded and self-adjoint and maps $L^2$ continuously into Schwartz space.
Therefore

$$
\mathcal A_i:=L_i+K_i
\quad\hbox{is self-adjoint on }D(L_i).                           \tag{9}
$$

This notation distinguishes the quantum factors $\mathcal A_i$ from
the classical matrix $A(z)$. A bounded symmetric perturbation does not
change the domain of a self-adjoint operator.

Their common domain is exactly the oscillator **form** domain

$$
\mathcal Q=
\{\psi\in L^2:\ x_j\psi\in L^2,
\ \partial_j\psi\in L^2\ \hbox{for all }j\}.
                                                                    \tag{10}
$$

Indeed $S_0$ is invertible, so the linear combinations $L_i$ span all
$\widehat q_j,\widehat p_j$. The statement holds distributionally as
well as on Schwartz functions, and each $K_i$ is bounded.
Let

$$
q_0[\psi]=\sum_i\|L_i\psi\|^2,\qquad
q[\psi]=\sum_i\|\mathcal A_i\psi\|^2,
\qquad M_K=\sum_i\|K_i\|^2.
$$

The smallest and largest eigenvalues $\lambda_-,\lambda_+$ of $A_0$
give, on $\mathcal Q$,

$$
\lambda_-\sum_j(\|x_j\psi\|^2+\|\widehat p_j\psi\|^2)
\leq q_0[\psi]\leq
\lambda_+\sum_j(\|x_j\psi\|^2+\|\widehat p_j\psi\|^2).          \tag{11}
$$

No commutation of coordinates and momenta is needed for (11): it is the
positive matrix inequality applied to a vector of Hilbert-space elements.
The elementary inequalities for $u+v$ then imply

$$
\frac12q_0[\psi]-M_K\|\psi\|^2
\leq q[\psi]\leq2q_0[\psi]+2M_K\|\psi\|^2.                    \tag{12}
$$

The norm $[q[\psi]+(M_K+1)\|\psi\|^2]^{1/2}$ is therefore equivalent
to the complete norm of $\mathcal Q$. The form $q$ is densely defined,
closed and nonnegative. Its representation theorem supplies a unique
nonnegative self-adjoint operator $\widehat H_\hbar$ associated with this
form. This is not an unjustified operator sum on an unspecified domain.

### An explicit operator-domain statement

Let $H_{0,\hbar}$ be the positive quadratic oscillator associated with
$q_0$; it is the Weyl quantization of $H_2$. For each $i$, $L_iK_i$ is
bounded, since multiplication or differentiation of a Schwartz kernel
again gives a Schwartz kernel. On $D(L_i)$, $K_iL_i$ has the bounded
extension $(L_iK_i)^*$. Hence

$$
V_\hbar=\sum_i\left[L_iK_i+(L_iK_i)^*+K_i^2\right]              \tag{13}
$$

is bounded and self-adjoint. Its kernel is Schwartz. The form identity
$q=q_0+\langle\psi,V_\hbar\psi\rangle$ implies the stronger result

$$
\boxed{
\widehat H_\hbar=H_{0,\hbar}+V_\hbar,
\qquad D(\widehat H_\hbar)=D(H_{0,\hbar}).
}                                                               \tag{14}
$$

For the already normalized candidate (2), this domain is the usual
oscillator operator domain

$$
D(H_{0,\hbar})=
\{\psi\in L^2:\ \Delta\psi\in L^2,
\ |x|^2\psi\in L^2\},                                         \tag{15}
$$

with distributional derivatives. Positive unequal frequencies give the
same domain. For general positive $H_2$ it is equivalently the domain
of its positive quadratic oscillator; (14) does not require choosing a
particular diagonalization. The form domain (10) is larger than (15).

The operator is still interacting. A bounded smoothing perturbation can
have the prescribed nonzero cubic, quartic and higher phase-space Taylor
coefficients. Also the $\mathcal A_i$ generally have nonconstant mutual
commutators, rather than the free oscillator algebra.

## 4. Compact resolvent, ground state, and heat trace

The inclusion $\mathcal Q\hookrightarrow L^2$ is compact. To prove it,
take a sequence bounded in the norm of (10). It is bounded in $H^1$
on every fixed configuration ball, so the local Rellich compactness
theorem gives a convergent subsequence there. Meanwhile

$$
\int_{|x|>R}|\psi(x)|^2dx
\leq R^{-2}\sum_j\|x_j\psi\|^2                                \tag{16}
$$

uniformly controls its tails. A diagonal subsequence and (16) give
compactness in all of $L^2$. The equivalent form norm (12) proves
compactness of the resolvent of $\widehat H_\hbar$.

It follows that its spectrum consists of finite-multiplicity eigenvalues
$0\leq E_0\leq E_1\leq\cdots\to\infty$, listed with multiplicity, and
that $E_0$ is attained by a normalizable eigenvector. Let
$\Omega\in\ker(\widehat H_\hbar-E_0)$ be any normalized vector.
No positivity-preserving position-space kernel has been proved, so a
Perron--Frobenius argument for a unique vacuum is unavailable.

One also obtains a finite Gibbs trace, which is stronger than compact
resolvent alone. Set
$K_{\rm osc}=\sum_j(\widehat p_j^2+x_j^2)$.
Equations (11)--(12) imply the form inequality

$$
\widehat H_\hbar\geq\frac{\lambda_-}{2}K_{\rm osc}-M_K I.
                                                                    \tag{17}
$$

The min--max principle and the known product oscillator spectrum give,
for $\tau>0$,

$$
\operatorname{Tr}e^{-\tau\widehat H_\hbar/\hbar}
\leq e^{\tau M_K/\hbar}
\sum_{\mathbf n\in\mathbb N^d}
 e^{-\tau\lambda_-(|\mathbf n|+d/2)}<\infty.                    \tag{18}
$$

All constants here may depend on the finite basis, the extension, and
$\hbar$. No uniform large-$d$ or cutoff estimate is asserted.
Trace class and finite-rank spectral projectors below any finite energy
are properties of this chosen confining completion; they are not a
microscopic area law or a finite total number of states.

## 5. Weyl symbol and the first quantum correction

Let $J^{ab}=\{z^a,z^b\}$, so in the $(q,p)$ order
$J=\begin{psmallmatrix}0&I\\-I&0\end{psmallmatrix}$.
The Weyl product convention corresponding to (8) is

$$
f\star_\hbar g
=f\exp\left[\frac{i\hbar}{2}
\overleftarrow\partial_a J^{ab}\overrightarrow\partial_b\right]g
$$

as an asymptotic symbol expansion. From (13)--(14), or from Weyl
composition on Schwartz space followed by closure, the exact symbol is

$$
h_\hbar=\sum_i a_i\star_\hbar a_i.
$$

Antisymmetry of $J$ cancels the odd terms in each self-product; reality
of $a_i$ ensures Hermiticity. The linear part in (7) has no second derivatives, and
$l_i\star k_i+k_i\star l_i=2l_ik_i$ exactly. Thus

$$
\boxed{
\begin{aligned}
h_\hbar(z)&=H_{\rm ext}(z)+\hbar^2r_2(z)+O(\hbar^4),\\
r_2(z)&=-\frac18\sum_i
J^{ab}J^{cd}\,a_{i,ac}(z)a_{i,bd}(z).
\end{aligned}}
                                                                    \tag{19}
$$

This confirms the proposed sign and coefficient for convention (8).
For one canonical pair, each summand contributes
$-\tfrac14(a_{i,qq}a_{i,pp}-a_{i,qp}^2)$ to $r_2$.
For example the self-product of the symbol $qp$ has correction
$+\hbar^2/4$, which independently fixes the sign convention.
This example only checks the Weyl formula; $qp$ is not one of the
required globally linear-plus-compact factors.

The expansion is not an assertion that the Moyal series converges.
Here it has controlled smooth remainders: $k_i\in\mathcal S$,
$k_i\star_\hbar k_i\in\mathcal S$, and the error in (19) is
$\hbar^4$ times a family bounded in Schwartz seminorms for
$0<\hbar\leq\hbar_0$. One way to see the even remainder directly is
to Fourier transform $k_i$. The Weyl-product phase is
$e^{-i\hbar\xi_aJ^{ab}\zeta_b/2}$; symmetrizing $\xi,\zeta$ replaces
it by its cosine. Taylor expansion of the cosine, with the rapidly
decaying Fourier factors, gives (19) and the remainder estimates.
The same estimates and kernel bounds give an $O(\hbar^4)$ bounded
operator remainder in this semiclassical statement.

The exact symbol need not equal $H$ near zero at fixed $\hbar$.
Equation (19) is the precise replacement: its principal symbol equals
$H$ there, its subprincipal symbol vanishes, and its first possible
ordering correction is of order $\hbar^2$. The exact remainder can
have Schwartz tails beyond the support of the classical modification.

## 6. Exactly what is preserved at one loop

Use the same canonical phase-space measure and midpoint/Weyl convention
as Section 7 of [the loop-kernel report](reduced-sixth-order-and-loop-kernel.md).
Do not change variables from $z$ to $a(z)$ in that measure. With
$h_\hbar=H_{\rm ext}+\hbar^2r_2+\cdots$, the formal effective action
has the grading

$$
\Gamma=S_{\rm ext}
+\frac{\hbar}{2}\operatorname{Tr}\log S_{\rm ext}^{(2)}
+O(\hbar^2),                                                   \tag{20}
$$

with the same phase-space kinetic convention and the usual chosen
vacuum subtraction. The $\hbar^2r_2$ term contributes first at order
$\hbar^2$ to the effective action, even as a tree insertion. It cannot
supply an order-$\hbar$ counterterm.

Every classical derivative of $H_{\rm ext}$ at zero equals that of
$H$. Therefore the free propagator and all local classical vertices
$H_3,H_4,\ldots$ are unchanged. In particular:

- The formal one-loop two-point coefficient built from $H_3,H_4$ is
  unchanged by this completion.
- The formal one-loop four-point kernel built from $H_3,\ldots,H_6$
  is unchanged, including the sixth-order tadpole and fifth/cubic
  contribution. This does not evaluate their missing gravitational inputs.
- No field-dependent order-$\hbar$ term is introduced by the
  sum-of-squares quantization. A constant ground-energy subtraction
  does not change the displayed vertex derivatives.

This is a theorem about the **formal loop expansion around the local
minimum**, not an equality of exact correlation functions, and not a
proof of convergence or numerical accuracy of that expansion.
At $\hbar=1$ an order-$\hbar^2$ term is not necessarily numerically small;
it must still be distinguished by its loop grading. Reordering with
normal ordering, adding an independent subprincipal counterterm, or
integrating out momenta with a different measure would require a new
comparison. Such changes have not been made here.

Different cutoffs with $\chi=1$ near zero leave all these local
classical Taylor jets unchanged. Their formal local symbol expansions
also coincide when the same square-root factorization is used there.
At fixed $\hbar$, Weyl quantization and the vacuum are not confined
to that neighborhood, so exact spectra and correlators may differ.
No exponential estimate for this dependence, or theorem equating
complete semiclassical spectral asymptotics, is needed or claimed here.
Other local factorizations of the same classical sum of squares can
also change the order-$\hbar^2$ correction.

### Preferred route: the direct Weyl Hamiltonian

There is a simpler construction that preserves the declared quantization
prescription exactly. Direct Weyl
quantization of $H_{\rm ext}=H_2+R$, $R\in C_c^\infty$, is
$H_{0,\hbar}+\operatorname{Op}^W_\hbar(R)$, a bounded self-adjoint
perturbation with compact resolvent and a finite lower bound. Subtracting
its ground energy also gives a positive transfer theory and preserves
the direct Weyl prescription exactly, apart from that scalar shift.
Explicitly, let

$$
\widehat H_\hbar^W=H_{0,\hbar}+\widehat R_\hbar,\qquad
D(\widehat H_\hbar^W)=D(H_{0,\hbar}),\qquad
E_0^W=\min\operatorname{spec}\widehat H_\hbar^W.
$$

The real compact smooth symbol R has the bounded Hermitian Schwartz
kernel proved above. The eigenvalues therefore obey

$$
\lambda_n(H_{0,\hbar})-\|\widehat R_\hbar\|
\leq\lambda_n(\widehat H_\hbar^W)
\leq\lambda_n(H_{0,\hbar})+\|\widehat R_\hbar\|.
$$

This also proves trace-class heat evolution for every positive time.
There is no field-dependent order-$\hbar^2$ modification for this route:
its Weyl symbol is exactly $H_{\rm ext}$. The subtraction uses its own
$E_0^W$, which need not equal the sum-of-squares ground energy.
Every nonzero-leg formal vertex remains unchanged by a scalar shift.

The sum-of-squares choice additionally guarantees nonnegativity before
this shift, while fixing the specific higher-order correction (19).
Neither choice is selected uniquely by the local gravitational data.

## 7. Interacting Euclidean evolution and vacuum reflection positivity

Choose the direct Weyl route for compatibility with the loop calculation.
The same argument also applies to the sum-of-squares alternative if its
own ground energy is used. Define

$$
K=\widehat H_\hbar^W-E_0^WI\geq0,\qquad
T(\tau)=e^{-\tau K/\hbar},\quad \tau\geq0,\qquad
U(t)=e^{-itK/\hbar}.                                           \tag{21}
$$

For the rest of this section, $\Omega$ denotes a normalized ground
vector of $\widehat H_\hbar^W$, independently chosen from its own ground
eigenspace. It is not assumed equal to a sum-of-squares ground vector.

For the explicitly constructed operator, these are respectively a
strongly continuous self-adjoint positive contraction semigroup and a
strongly continuous unitary group. They use the same spectrum and
spectral measure. Matrix elements of $e^{-zK/\hbar}$ are analytic for
$\operatorname{Re}z>0$ and have the appropriate strong boundary values
on the imaginary axis. Thus Euclidean and Lorentzian evolution are
obtained from one interacting quantum model.

Here is a vacuum Osterwalder--Schrader Gram construction requiring no
position-space path-integral measure. Use bounded observables and write
Euclidean products with increasing times from left to right. For a
positive-time history $F=(O_1,t_1;\ldots;O_n,t_n)$ with
$0\leq t_1<\cdots<t_n$, define linearly

$$
v_F=T(t_1)O_1T(t_2-t_1)O_2\cdots
T(t_n-t_{n-1})O_n\Omega.                                      \tag{22}
$$

Time reflection reverses the order and takes adjoints. For a second
history $G=(B_1,s_1;\ldots;B_m,s_m)$, the corresponding reflected
Euclidean vacuum pairing is

$$
\begin{aligned}
\mathcal S(\theta F,G)
=\langle\Omega,&\,O_n^\dagger T(t_n-t_{n-1})\cdots
O_1^\dagger T(t_1+s_1)\\
&\times B_1T(s_2-s_1)\cdots B_m\Omega\rangle
=\langle v_F,v_G\rangle.
\end{aligned}                                                    \tag{23}
$$

The semigroup law supplies the middle interval. Consequently

$$
\boxed{\sum_{i,j}\overline c_i c_j\,
\mathcal S(\theta F_i,F_j)
=\left\|\sum_jc_jv_{F_j}\right\|^2\geq0.}                       \tag{24}
$$

For example, with a single bounded $O$,

$$
C_E(\tau)=\langle\Omega,O^\dagger T(\tau)O\Omega\rangle,
\qquad
C_L(t)=\langle\Omega,O^\dagger U(t)O\Omega\rangle
$$

are matrix elements of the same spectral measure, and
$C_E(\tau_i+\tau_j)$ is a positive Gram kernel. Vacuum degeneracy
causes no failure of (22)--(24). It can affect clustering and the choice
of state, which are different requirements.

Quotienting the history space by its zero-norm vectors and completing
recovers $\overline{\operatorname{span}\{v_F\}}$, the cyclic space
generated by these histories and their time translations.
Unbounded composite observables require their own domain checks; the
bounded version already gives an exact positive reconstruction.
Equation (18) also permits thermal reflection-positive trace pairings,
but a thermal trace is not needed for the vacuum construction.

“Positive” in (21)--(24) means Hilbert-space/operator and reflection
positivity. A general momentum-dependent Weyl Hamiltonian need not
preserve pointwise nonnegative wavefunctions. These results therefore
do not construct a positive local Boltzmann measure, a stochastic
transition kernel, or an identification with an archived Monte Carlo
measure. They establish quantum-mechanical time reconstruction, not
spatial locality, Lorentz covariance, or all relativistic field-theory
Osterwalder--Schrader axioms.

## 8. What this completion adds, and what it leaves open

This is more than the earlier conditional statement “if a positive
self-adjoint interacting Hamiltonian exists, then its spectral theorem
provides transfer.” Equations (3)--(18) construct such a Hamiltonian and
prove its domains, positivity, confinement and vacuum existence from a
local finite-dimensional classical input. Nonzero local interactions
are preserved, and (19)--(20) connect its formal one-loop prescription
to the one already used for the reduced candidate.

Its nonuniqueness is equally concrete. The cutoff, the extension outside
the local chart, and possible higher-order quantization choices are new
global inputs. The original reconstructed Hamiltonian need not equal
$H_{\rm ext}$ outside the matching neighborhood, or even exist there.
The positive theorem cannot be used to exponentiate the unbounded
quartic Taylor polynomial as if it were this completion.

The [local-corner analysis](local-corner-matching.md) supplies reductions
that, under their stated physical choices, return the same canonical
bulk form. Applying the present theorem after such a reduction gives a
positive interacting **bulk** quantum mechanics with no independently
derived corner multiplicities. It does not prove that the naive local
matching constraints are first class, select the correct gravitational
edge reduction, or construct the full normal-frame charge algebra.

Still unproved are the original global/microscopic Hamiltonian,
finite gravitational state counting or an area law, a physical
large-angular/radial-rank limit, the missing evaluated four-graviton
loop observable, and equivalence to interacting continuum gravity.
The local constraint-existence argument is supplied separately in
[the weighted-space proof](local-existence-and-adm-smoothness.md).
It does not select the global extension used here. No numerical
or uniform-in-$N$ bound on the permitted extension radius follows from
the qualitative finite-mode theorem.
