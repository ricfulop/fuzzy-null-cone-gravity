# Independent audit of sixth-order reduction and the loop prescription

## Scope and verdict

The primary audit is against source commit
`f220c7eb5bbef322c6367553a2f42df51dc4c8ac`. It compares
[the sixth-order report](reduced-sixth-order-and-loop-kernel.md) with
[the canonical construction](reduced-canonical-gaunt.md) and
[the third/fourth-order derivation](reduced-vertices.md). Line numbers below
refer to the sixth-order report at that commit. Paper II's README was read.
This is an independent algebraic audit, not a rerun of the historical
calculations or a numerical validation of general gravitational H5/H6.

No algebraic sign, normalization, or combinatorial error was found in the
displayed H5/H6, York, pointwise scalar, phase propagator, or partition
formulas, **conditional on the stated decaying inverses, smooth positive
conformal branch, and invertible Gram matrix**. The derivations below specify
what was checked rather than treating the whole construction as approved.

One consequential scope clarification is needed: four derivatives of the
static trace-square-root expression test **four configuration legs at
zero external frequency and zero background momentum**. They do not test
arbitrary zero-frequency phase-space legs. In particular, the static check
cannot test vertices with four or more momentum factors. The original
wording at lines 427–428 should name this restriction explicitly.

Other implementation distinctions that must survive are: two York solves
mean two perturbative orders per supplied direction, not two right-hand
sides for the entire vertex tensor; eliminating W3 does not eliminate b3;
and the angular numbers are sufficient upper bounds, not a proof of radial
accuracy or a finite-amplitude angular cap.

At the parent agent's subsequent request, Section 10 also audits the
27-invariant representation and analytic two-point formulas in the newly
added source files. Those files are not part of the pinned primary source.
Their numerical execution and coefficient values are outside this audit.

## 1. Scalar recursion and H5/H6: confirmed

Let s denote the full scalar sigma squared, rather than only its leading
term. The defining constraint is

$$
8\bar\Delta(1+u)-\bar R(1+u)+s(1+u)^{-7}=0.
$$

Using `(1+u)^(-7)=1-7u+28u^2-84u^3+...` and
$\mathcal C=(-8\Delta)^{-1}$ gives

$$
u=\mathcal C(s-\bar R)+8\mathcal C(\bar\Delta-\Delta)u
-\mathcal C[(\bar R+7s)u]+28\mathcal C(su^2)+O(su^3).
$$

Because R and s start at total field order two, collecting orders two,
three, and four gives exactly the three displayed responses at lines
43–46. In particular, the sign of the D1 term is positive, whereas the
order-four B2 u2 term has a minus sign. Unimodularity gives
$\bar\Delta=\partial_i[(e^{-Q})^{ij}\partial_j]$, so the stated Dn
also has the correct sign and factorial.

The exact energy expands as

$$
H=\int(s-\bar R)-\int(\bar R+7s)u+28\int su^2+O(\lambda^8).
$$

Consequently

$$
H_5=I_5-\int R_5-\langle B_3,u_2\rangle-\langle B_2,u_3\rangle,
$$

$$
H_6=I_6-\int R_6-\langle B_4,u_2\rangle
-\langle B_3,u_3\rangle-\langle B_2,u_4\rangle
+28\int s_2u_2^2.
$$

These are ordinary homogeneous coefficients, without an implicit n!
normalization. The first possible u5 contribution has order seven. There
is no missing u5 or u6 term at the orders asserted.

At orders three and four this same expansion gives
$H_3=-\int R_3$ and
$H_4=I_4-\int R_4-\langle R_2+7s_2,\mathcal C(s_2-R_2)\rangle$,
matching the earlier vertex report. Pointwise sources cannot be replaced
by densities equal only after spatial integration.

## 2. Inverse Gram expansion: confirmed

For $A=\mathsf G_0^{-1}$ and
$E=\mathsf G_2+\mathsf G_3+\mathsf G_4+\cdots$,

$$
(\mathsf G_0+E)^{-1}=A-AEA+AEAEA+\cdots.
$$

Multiplication by p on both sides, with $a=Ap$, produces

$$
\begin{aligned}
I_2&=p^TAp, & I_3&=0,\\
I_4&=-a^T\mathsf G_2a, & I_5&=-a^T\mathsf G_3a,\\
I_6&=-a^T\mathsf G_4a+a^T\mathsf G_2A\mathsf G_2a.
\end{aligned}
$$

The order of the noncommuting matrices in I6 is correct. This uses the
exact canonical identity $\int\sigma^2=p^T\mathsf G^{-1}p$, not an
assumption that the coefficients of S are canonical momenta. It legitimately
avoids pointwise s5 and s6.

## 3. York subtraction, N4, and the raw fourth-order norm: confirmed

The tangent expansion follows from the Frechet derivative of the ordinary
matrix exponential:

$$
D\exp_Q[S]=\sum_{n\geq0}\frac1{(n+1)!}
\sum_{a=0}^nQ^aSQ^{n-a}.
$$

It is trace-free with respect to $\bar\gamma$ because
$D\log\det(e^Q)[S]=\operatorname{tr}S=0$. With contravariant W, the
longitudinal tensor is its metric Lie derivative minus
$(2/3)\bar\gamma\,\bar D_kW^k$. Since det of the conformal metric is
one, $\bar D_kW^k=\partial_kW^k$; this proves lines 111–115.

For a trace-free covariant tensor T, integration by parts gives

$$
\int V^i\bar D^jT_{ij}
=-\frac12\int(L_{\bar\gamma}V)_{ij}T^{ij}.
$$

The pairing here is the natural vector/covector contraction with coordinate
volume. It is not a separately metric-weighted vector inner product. Thus
each coefficient of $\mathcal A=\mathcal D L$ is symmetric in that
fixed pairing, and
$\mathcal N=-2\langle b,\mathcal A^{-1}b\rangle$ has the correct sign.
The flat operator is negative on the decaying domain, so N2 is nonnegative.

To check N4 without assuming the claimed cancellation, temporarily retain
the third response:

$$
\begin{aligned}
\mathcal A_0W_1&=b_1,\\
\mathcal A_0W_2&=b_2-\mathcal A_1W_1=r_2,\\
\mathcal A_0W_3&=b_3-\mathcal A_1W_2-\mathcal A_2W_1.
\end{aligned}
$$

The q-fourth coefficient of $\langle b,W\rangle$ is

$$
\begin{aligned}
\langle b_1,W_3\rangle+\langle b_2,W_2\rangle
+\langle b_3,W_1\rangle
&=2\langle b_3,W_1\rangle
+\langle b_2-\mathcal A_1W_1,W_2\rangle
-\langle W_1,\mathcal A_2W_1\rangle\\
&=2\langle b_3,W_1\rangle
+\langle r_2,\mathcal A_0^{-1}r_2\rangle
-\langle W_1,\mathcal A_2W_1\rangle.
\end{aligned}
$$

Multiplication by -2 gives precisely the asserted N4. The same calculation
one order lower gives N3. There is no missing A1-square term: it is inside
the r2 quadratic form. b3 remains necessary even though W3 is not solved.
The index on Nn counts powers of Q at fixed S; N4 is total field order six.

For an independent raw-norm derivation, diagonalize Q pointwise, with
eigenvalues $\lambda_i$. Then

$$
(E_S)_{ij}=S_{ij}\frac{e^{\lambda_i}-e^{\lambda_j}}
{\lambda_i-\lambda_j},
$$

using the continuous limit on the diagonal. Hence

$$
\operatorname{tr}[(e^{-Q}E_S)^2]
=\sum_{ij}S_{ij}^2
\left[\frac{\sinh((\lambda_i-\lambda_j)/2)}
{(\lambda_i-\lambda_j)/2}\right]^2.
$$

The multiplier is $1+x^2/12+x^4/360+O(x^6)$. Since
$\operatorname{tr}[Q,S]^2=-\sum_{ij}x_{ij}^2S_{ij}^2$ and
$\operatorname{tr}[Q,[Q,S]]^2=\sum_{ij}x_{ij}^4S_{ij}^2$,
the reported M2 and M4 coefficients, including their signs, follow.
All odd raw coefficients vanish. There is no analogous reason for N3,
and therefore G3, to vanish.

The Gram identities $\mathsf G_n=\mathcal M_n-\mathcal N_n$ follow
from orthogonal York projection. Polarization is legitimate, since these
expressions remain quadratic in the supplied S. Computing the full matrix
still requires sufficient independent right-hand sides or polarizations.

## 4. Pointwise s3 and s4: confirmed, with a useful sign check

The actual coefficient vector in sigma is
$\mathsf G^{-1}p=a-\mathsf G_0^{-1}\mathsf G_2a+O(q^3p)$.
The resulting flat-tensor correction first enters Z2, exactly with the
minus sign at lines 206–208. Corrections to its curved tangent or York
projection begin one order higher.

As mixed-index matrices,

$$
e^{-Q}\sigma=S+U_1+U_2+\cdots,
\qquad \sigma_{ij}\sigma^{ij}=\operatorname{tr}[(e^{-Q}\sigma)^2].
$$

Multiplying the matrix series proves both U formulas and

$$
s_3=2\operatorname{tr}(SU_1),\qquad
s_4=\operatorname{tr}(U_1^2)+2\operatorname{tr}(SU_2).
$$

In particular, write $A=[Q,S]$ and $F=L_0W_1$. Then
$U_1=-A/2-F$. A is antisymmetric and F is symmetric, giving

$$
s_3=-2S:F,\qquad
\operatorname{tr}(U_1^2)=\frac14\operatorname{tr}A^2+\operatorname{tr}F^2.
$$

This is a matrix-square trace, **not** the Euclidean Frobenius square
$\operatorname{tr}(U_1^TU_1)$. Replacing it by the latter would flip
the commutator term and be an implementation error.

Flat TT orthogonality gives $\int S:L_0W_1=0$, proving I3=0 without
setting its integrand to zero. Before the dual-Gram correction, the
integrated q-squared coefficient of the projected norm is
$a^T\mathsf G_2a$. The Z2 correction contributes
$-2a^T\mathsf G_2a$ through $2S:Z_2$, leaving
$\int s_4=-a^T\mathsf G_2a$. This independently checks the sign and
factor of the canonical correction.

## 5. Pure-momentum and configuration reductions: confirmed

At q=0, $u_2=\mathcal C s_2$, $u_3=0$, and
$u_4=-7\mathcal C(s_2u_2)$. H6 becomes

$$
49\langle s_2,\mathcal C(s_2u_2)\rangle+28\int s_2u_2^2
=(49+28)\int s_2(\mathcal C s_2)^2.
$$

Self-adjointness of C accounts for the 77. The odd fifth order vanishes.
This confirms the relation to the earlier pure-momentum formula; it is
not an independent numerical evaluation of that historical coefficient.

For p=0, $\mathcal L_Q\psi_0=0$ and
$w=\psi_0-1=-\mathcal L_Q^{-1}\bar R$. Therefore
$V=-\int\bar R+\langle\bar R,\mathcal L_Q^{-1}\bar R\rangle$.
Set $A_1=-8D_1$ and $A_2=-8D_2+R_2$, where R2 here acts by
multiplication. The inverse through q-squared is

$$
\mathcal L_Q^{-1}
=\mathcal C-\mathcal C A_1\mathcal C
+\mathcal C A_1\mathcal C A_1\mathcal C
-\mathcal C A_2\mathcal C+O(q^3).
$$

With $v=\mathcal C R_2$, the q-sixth terms before completing the square
are

$$
2\langle R_4,v\rangle+\langle R_3,\mathcal C R_3\rangle
+16\langle R_3,\mathcal C D_1v\rangle
+64\langle D_1v,\mathcal C D_1v\rangle
+8\langle v,D_2v\rangle-\int R_2v^2.
$$

They give exactly the reported V6 after adding $-\int R_6$. The
q-fifth terms give the reported V5. Integration by parts gives the
stated positive sign for $\langle v,D_1v\rangle$ and negative sign for
$\langle v,D_2v\rangle$.

## 6. Angular support: the upper bounds are valid

A rotationally equivariant flat derivative maps the appropriate tensor
representation to another tensor representation without changing total j.
The decaying scalar and vector Green operators likewise preserve total j.
Products of n free tensors have total j at most nL. Recursion therefore
gives u2/u3/u4 upper bounds 2L/3L/4L. W1 is bilinear in Q and S;
W2 is trilinear. Their bounds are 2L and 3L, respectively.

The finite-basis dual-Gram correction does not spoil the spatial bound:
it is a linear combination of retained tensors, each of total j at most
L, with spatially constant coefficients. The complete selected multiplets
also preserve rotational equivariance in coefficient space.

Thus for L=2, scalar solves through j=8 and vector solves through j=6
are sufficient for the stated procedure. This is not a claim that every
allowed j is nonzero or that all those channels must contribute to every
final contraction. b3 can carry j up to 8 even though no W3 solve is
needed. The bound concerns total j, not the orbital harmonic index of an
individual Cartesian component. A componentwise scalar implementation
must respect that distinction.

The statement proves finite angular bookkeeping at fixed perturbative
order. It supplies neither radial convergence nor a finite-dimensional
radial closure after inverse solves, and it does not cap the exact
finite-amplitude geometry.

## 7. Euclidean propagator and ordered partitions: confirmed

With Fourier convention $e^{i\nu\tau}$ and kinetic term
$-iP\dot Q$, the symmetrized quadratic Hessian is

$$
\mathcal K_i(\nu)=
\begin{pmatrix}\omega_i^2&-\nu\\\nu&1\end{pmatrix}.
$$

Direct multiplication shows

$$
\mathcal K_i(\nu)^{-1}
=\frac1{\nu^2+\omega_i^2}
\begin{pmatrix}1&\nu\\-\nu&\omega_i^2\end{pmatrix}.
$$

In particular $\mathcal G(\nu)^T=\mathcal G(-\nu)$. Changing just
one off-diagonal sign would violate the inverse identity. Completing the
momentum square also gives
$G_{PP}=1-\nu^2/(\nu^2+\omega^2)=\omega^2/(\nu^2+\omega^2)$;
the contact term and induced velocity term must be combined. This is
consistent with the declared phase-space prescription, not a positive
real joint Q/P sampling covariance.

The one-loop functional is formally

$$
\frac\hbar2\operatorname{Tr}\log(\mathcal K_0+\mathcal B)
=\text{constant}+\frac\hbar2\sum_{k\geq1}
\frac{(-1)^{k+1}}k\operatorname{Tr}(\mathcal G\mathcal B)^k.
$$

Each nonempty external block differentiates one Hessian factor; at the
zero background it selects H of degree block-size plus two. Four
distinct derivatives distributed among k nonempty factors are precisely
ordered set partitions. There is no additional block factorial, since
the B definitions already differentiate labeled external legs.

| k | Block sizes | Ordered count | Count after cyclic identification |
| --- | --- | ---: | ---: |
| 1 | 4 | 1 | 1 |
| 2 | 1+3 | 8 | 4 |
| 2 | 2+2 | 6 | 3 |
| 3 | 1+1+2 | 36 | 12 |
| 4 | 1+1+1+1 | 24 | 6 |

The report's totals 1,14,36,24 and logarithm factors are correct. In the
commuting single-coordinate limit the resulting integrand is

$$
\frac\hbar2\left[
g\,V^{(6)}-g^2\{4V^{(3)}V^{(5)}+3(V^{(4)})^2\}
+12g^3(V^{(3)})^2V^{(4)}-6g^4(V^{(3)})^4\right],
$$

which is independently the fourth derivative of
$\hbar\log(\nu^2+V'')/2$. This checks all five vertex families and
their multiplicities without gravitational coefficient data.

For a background block with frequency $e^{-i\Omega_S\tau}$,
$\langle\nu|B_S|\mu\rangle$ has support
$\mu=\nu+\Omega_S$. Therefore the report's ordered-product routing
$\nu_a=\nu+\sum_{b<a}\Omega_{S_b}$ has the correct sign. Closure of
the trace uses the external frequency-conservation condition.

At fixed finite mode number and positive free frequencies, these particular
frequency integrals have no ultraviolet divergence: for k at least two,
the product is at worst order $|\nu|^{-k}$; for k=1, symmetry of B
cancels the antisymmetric order-$1/\nu$ part of G, leaving order
$1/\nu^2$. This elementary fixed-mode observation does not control
growth of mode sums as either physical cutoff is removed.

The formula is the determinant contribution with the declared Weyl
prescription. Any separately chosen order-hbar Hamiltonian term must be
added at that order; it is not produced automatically by the determinant.
No four-point gravitational number, on-shell scattering limit, or
interacting reconstruction follows just from this correct bookkeeping.

## 8. Static kinetic matrix and determinant: confirmed, but restricted

At fixed q, sigma is linear in p. For $p\mapsto tp$, the constraint
through t-squared is

$$
\mathcal L_Q\chi=\sigma_1^2\psi_0^{-7},\qquad
\mathcal L_Q\psi_0=0.
$$

Green's identity gives

$$
\int\sigma_1^2\psi_0^{-6}
=-8\lim_{r\to\infty}\oint
(\psi_0\partial_r\chi-\chi\partial_r\psi_0)
=-8\lim_{r\to\infty}\oint\partial_r\chi.
$$

Here $\psi_0=1+O(1/r)$ and $\chi=O(1/r)$ make the extra surface
terms vanish; the rapidly decaying conformal metric introduces no new
ADM surface term. Consequently

$$
H_{p^2}(q,p)=p^T\mathsf G^{-1}\mathsf J\mathsf G^{-1}p,
\quad
\mathsf J_{\alpha\beta}=\int\psi_0^{-6}T_{\alpha ij}T_\beta^{ij},
$$

so $\mathsf M=2\mathsf G^{-1}\mathsf J\mathsf G^{-1}$, including
the factor two and exponent -6, is correct. Indices on T are raised by
the conformal metric. Positive psi0 and positive G also imply positive J
and M in the stated domain; positivity of the potential Hessian is an
additional condition.

There is a lower-order consistency check. Since $\psi_0=1-v+O(q^3)$,
the q-squared correction to the momentum energy is

$$
-a^T\mathsf G_2a+6\int v\,S:S,
$$

which is exactly the earlier mixed quartic functional. At q=0,
$\mathsf J=\mathsf G_0$ and $\mathsf M=2\mathsf G_0^{-1}$,
recovering the original quadratic kinetic normalization.

At a static background (q,0), evenness in p sets the mixed Hessian to
zero. The phase Hessian is

$$
\begin{pmatrix}\mathsf K&-\nu I\\\nu I&\mathsf M\end{pmatrix},
\quad \mathsf K=V''(q),
$$

with determinant

$$
\det\mathsf M\,\det(\mathsf K+\nu^2\mathsf M^{-1})
=\det(\nu^2I+\mathsf M\mathsf K).
$$

For positive M and K, MK has the same positive eigenvalues as
$\mathsf M^{1/2}\mathsf K\mathsf M^{1/2}$. Frequency integration of
the logarithm gives the stated
$\hbar\operatorname{tr}\sqrt{\mathsf M^{1/2}\mathsf K\mathsf M^{1/2}}/2$
up to a background-independent subtraction and specified counterterms.
Adding an extra field-dependent determinant after this Schur complement
would double count the momentum integration.

**Necessary scope repair:** this function knows only V(q) and the p-squared
coefficient M(q). Four q derivatives test the qqqq vertex at zero external
frequency. They cannot see $q p^4$ in H5, or $q^2p^4$ and $p^6$ in H6:
two internal momentum derivatives still leave positive powers of p at
p=0. General phase-space four-point checks require external momentum legs
or nonzero momentum backgrounds as well. Nonzero-frequency checks are
separately needed even for four configuration legs.

## 9. Concrete acceptance requirements and audit checks actually performed

The owning report should retain all formulas checked above. Before treating
its prescription as implemented, the following distinctions should be
explicit in its acceptance record:

1. Name the static check as qqqq, at p=0 and all external frequencies zero;
   add separate momentum and mixed-leg checks. Do not infer full H6
   validation from the 77 control plus the static potential alone.
2. Keep b3 in G4 even though W3 is not solved. Record how the full retained
   G2 contraction and matrix polarizations are obtained.
3. Compute s4 as a trace of a matrix square, with the dual-Gram correction;
   compare its integral with the independent inverse-Gram expression.
4. Distinguish total-j bounds from component orbital degrees and from
   radial/tail error. Verify all pointwise scalar sources, rather than
   substituting integrated curvature densities.
5. Treat nonlinear canonical-coordinate controls carefully: fixed-coordinate
   off-shell effective-action derivatives are not automatically invariant
   numbers under nonlinear redefinitions. A comparison must specify the
   transformed sources, observable, and quantum ordering, or use linear
   symplectic changes as the simpler control.

In addition to the analytic derivations, one inline standard-library Python
sanity check used exact rational arithmetic and left no calculation output
file. With

$$
\mathcal A_0=\operatorname{diag}(-2,-3),\quad
\mathcal A_1=\begin{pmatrix}1&2\\2&-1\end{pmatrix},\quad
\mathcal A_2=\begin{pmatrix}2&-1\\-1&1\end{pmatrix},
$$

and $b_1=(1,2)$, $b_2=(-2,1)$, $b_3=(3,-1)$, the explicit W3 expression
and reduced residual expression both give N4=551/108. These matrices do
not commute, so this check does not accidentally assume scalar inverses.
The same lightweight check recovered 1,1/12,1/360 from the sinh multiplier
and enumerated the ordered-partition counts and block sizes above.
These are algebraic sanity checks, not gravitational experiments.

No new gravitational solve, tracked run, numerical H5/H6 evaluation,
build, branch change, or commit was performed for this audit. The source
report's statements about earlier successful runs were not revalidated.

## 10. Supplement: quartic completeness and the proposed two-point formulas

This section responds to a later request to inspect
[complete_quartic.py](../../verification/quantum_bridge/complete_quartic.py)
and [quartic_loop_check.py](../../verification/quantum_bridge/quartic_loop_check.py)
as newly added working files. The representation and analytic factors can
be checked without executing their gravitational coefficient calculation.

### The count 27 and the five cross invariants are correct

Let $A=q_E$, $B=q_B$, $C=S_E$, $D=S_B$ be STF spin-two label matrices.
Time reversal allows only even momentum degree. Spatial parity allows only
even total magnetic-species degree. These restrictions split H4 into:

| Sector | Independent invariant count | Reason |
| --- | ---: | --- |
| q4 | 5 | One E4, one B4, three E2B2 |
| S4 | 5 | Same representation decomposition |
| qE2 SE2, qE2 SB2, qB2 SE2, qB2 SB2 | 12 | Three invariants for each of four species pairs |
| qE qB SE SB | 5 | Four distinct spin-two factors |
| Total | 27 | All parity- and time-reversal-allowed sectors |

For a repeated pair, $\operatorname{Sym}^2(V_2)=V_0\oplus V_2\oplus V_4$,
so coupling two such pairs gives three singlets. A trace basis is
$\operatorname{tr}A^2\operatorname{tr}C^2$,
$(\operatorname{tr}AC)^2$, and $\operatorname{tr}(A^2C^2)$.
For four distinct factors, $V_2\otimes V_2$ contains each spin
0,1,2,3,4 once, hence there are five singlets in the fourfold product.

The six natural trace expressions are the three double traces and three
single traces modulo cyclicity and reversal. Polarizing the 3x3 traceless
Cayley–Hamilton identity $\operatorname{tr}X^4=(\operatorname{tr}X^2)^2/2$
gives their relation:

$$
\operatorname{tr}(ABCD)+\operatorname{tr}(ABDC)+\operatorname{tr}(ACBD)
=\frac12\left[
\operatorname{tr}(AB)\operatorname{tr}(CD)
+\operatorname{tr}(AC)\operatorname{tr}(BD)
+\operatorname{tr}(AD)\operatorname{tr}(BC)\right].
$$

Thus the three double traces plus TrABCD and TrACBD used in the new file
are an appropriate five-dimensional basis. The omitted TrABDC is
determined by this identity. Exact full-rank direction selection in the
implementation is a useful additional independence check; numerical
coefficient reconstruction has not been executed in this audit.

The code's canonical factors also agree with the original normalization:
$q_\alpha=\sqrt{2/g_\alpha}\,Q_\alpha$ and the coefficient of S is
$p_\alpha/g_\alpha=P_\alpha/\sqrt{2g_\alpha}$. S is not the
unnormalized canonical p. Its local commutator factor 1/24 matches the
prior implementation's units: the local Gaussian integral unit is half
the common quartic unit, so it implements the functional's 1/12.

This establishes the completeness of the allowed representation and the
normalization map, not the correctness of as-yet-unrun fitted coefficients.

### Tadpole and cubic bubble factors are correct

Write $B_{12}=\partial_1\partial_2 H''_4(0)$ and
$C_a=\partial_aH''_3(0)$. The two-leg version of the trace-log expansion is

$$
\frac{\Gamma^{(1)}_{12}}\hbar
=\frac12\int\frac{d\nu}{2\pi}\operatorname{tr}(G B_{12})
-\frac12\int\frac{d\nu}{2\pi}
\operatorname{tr}[G(\nu)C_1G(\nu+\Omega)C_2].
$$

For this candidate C has only configuration entries. Symmetry of B cancels
the antisymmetric mixed propagator in the tadpole. Since
$\int d\nu/[2\pi(\nu^2+\omega_i^2)]=1/(2\omega_i)$, its contribution is

$$
\mathrm{Tad}_{12}=\sum_i\left[
\frac{(B_{12})_{Q_iQ_i}}{4\omega_i}
+\frac{\omega_i(B_{12})_{P_iP_i}}4\right].
$$

The scalar two-propagator integral is

$$
\int\frac{d\nu}{2\pi}
\frac1{(\nu^2+\omega_i^2)((\nu+\Omega)^2+\omega_j^2)}
=\frac{\omega_i+\omega_j}
{2\omega_i\omega_j[(\omega_i+\omega_j)^2+\Omega^2]}.
$$

One independent derivation multiplies the time-domain kernels
$e^{-\omega_i|t|}/(2\omega_i)$ and
$e^{-\omega_j|t|}/(2\omega_j)$ and Fourier transforms their product.
Therefore

$$
\mathrm{Bubble}_{12}=-\frac12\sum_{ij}
(C_1)_{Q_iQ_j}(C_2)_{Q_jQ_i}
\frac{\omega_i+\omega_j}
{2\omega_i\omega_j[(\omega_i+\omega_j)^2+\Omega^2]}.
$$

This agrees with the new analytic routine. The sums include all ordered
internal mode pairs; adding another factor for i unequal to j would double
count. These expressions are coefficients of hbar in normalized variables.

### The configuration inverse correction needs the mixed-block condition

If the corrected phase Hessian is

$$
\begin{pmatrix}
\omega^2+\hbar\Pi_{QQ}&-\Omega+\hbar\Pi_{QP}\\
\Omega+\hbar\Pi_{PQ}&1+\hbar\Pi_{PP}
\end{pmatrix},
$$

its configuration Schur complement through one loop is

$$
\omega^2+\Omega^2+
\hbar[\Pi_{QQ}-\Omega^2\Pi_{PP}
+\Omega\Pi_{PQ}-\Omega\Pi_{QP}].
$$

For the present H3/H4, the mixed one-loop blocks vanish: a momentum
external leg kills the cubic bubble, while a QP external differentiation
of the even-momentum H4 leaves only mixed internal Hessian entries in
the tadpole, whose contraction with G vanishes by symmetry. Thus the
reported simplification $\Pi_{QQ}-\Omega^2\Pi_{PP}$ is valid at this
order for this candidate. Direct mixed-leg checks remain worthwhile before
using that simplification in numerical results. It is not a universal
identity for arbitrary momentum-dependent cubic Hamiltonians.

These two-point formulas use only H3/H4. Successfully evaluating them
would constitute a finite-mode two-point loop calculation for the declared
candidate, not the requested four-graviton amplitude. Nothing in this
audit closes the general H5/H6 evaluation, local corner-state construction,
microscopic state counting, or interacting Euclidean–Lorentzian
reconstruction.
