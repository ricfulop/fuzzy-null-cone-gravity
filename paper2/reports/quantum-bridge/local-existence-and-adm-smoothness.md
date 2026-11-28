# Local existence and smooth ADM energy for the finite Gaunt chart

## Result and scope

For any fixed finite family of linearly independent real Schwartz TT
profiles on flat three-space, the reduction-first construction has a
smooth positive conformal solution and a positive invertible York Gram
matrix in a sufficiently small neighborhood of the Minkowski coefficient
origin. The reduced ADM Hamiltonian is smooth there, has zero value and
gradient at the origin, and has positive definite Hessian.

This fills the local analytic hypothesis used in
[the reduced canonical candidate](reduced-canonical-gaunt.md) for this
specified class of profiles. It uses the standard weighted elliptic
isomorphisms below and a finite-parameter implicit-function argument.
It is an analytic derivation, not an additional numerical result.
It supplies no global coefficient domain, no uniform neighborhood as the
angular or radial rank grows, and no invariance of the finite image under
the full Einstein flow.

In particular, this result can supply the local Hamiltonian input to
[the positive quantum completion](positive-local-quantum-completion.md).
That completion still makes a separate choice of global dynamics and
quantum ordering beyond the local classical chart.

## 1. Spaces and the standard elliptic input

Use Cartesian coordinates on all of R3, no punctures or inner boundary,
and write $\langle r\rangle=(1+r^2)^{1/2}$. Use the following weighted
Sobolev convention, which is the convention of the source cited below:

$$
\|u\|_{W^4_{k,\delta}}^4
=\sum_{|\alpha|\leq k}\int_{\mathbb R^3}
 |\partial^\alpha u|^4
 \langle r\rangle^{4(\delta+|\alpha|)}\,d^3x .
$$

Fix an integer $k\geq5$ and set

$$
X=W^4_{k,0},\qquad Y=W^4_{k-2,2}.
$$

Tensor and vector spaces have the corresponding component norms.
The scalar Laplacian and the flat York vector operator are isomorphisms
from X to Y:

$$
\Delta:X\longrightarrow Y,\qquad
\mathcal A_0 W=\Delta W+\frac13\nabla\operatorname{div}W:
X^3\longrightarrow Y^3.
$$

The weight lies strictly in the allowed interval
$-3/4<\delta<1/4$. Decaying scalar harmonic functions on all of R3
vanish; the decaying York kernel is also zero. The latter follows from
integration by parts and the absence of a decaying conformal Killing
field. No constant or growing gauge vector is admitted in X.

The weighted convention and multiplication/embedding estimates are
given in [Choquet-Bruhat, Isenberg and York, PDF p. 5](https://www.alphaxiv.org/pdf/gr-qc/9906095v2?page=5).
The vector solvability statement is on
[PDF p. 6](https://www.alphaxiv.org/pdf/gr-qc/9906095v2?page=6);
the elliptic isomorphism and scalar injectivity statements are in
[Appendices A--B, PDF pp. 24--25](https://www.alphaxiv.org/pdf/gr-qc/9906095v2?page=24).
Only these standard analytic inputs are imported. Their application to
the particular finite dual-Gram chart is derived below.

For example, the embedding allows the decay exponent $\beta=2/3$:

$$
|\partial^j u(x)|\leq C\|u\|_X\langle r\rangle^{-\beta-j},
\qquad 0\leq j\leq2 .
$$

The strict inequality $\beta>1/2$ will make the ADM surface remainder
vanish. Its strict inequality $\beta<3/4$ respects the source's embedding
range rather than asserting an endpoint estimate.

## 2. Smooth York projection and invertible Gram matrix

Let $Q(q)=\sum_\alpha q^\alpha e_\alpha$, with the finite Schwartz
profiles trace-free and divergence-free in the Euclidean metric, and put

$$
\bar\gamma=e^Q,\qquad
E_\alpha=\partial_{q^\alpha}\bar\gamma .
$$

The map $q\mapsto\bar\gamma-I$ is smooth into every fixed weighted space
needed here. The conformal metric is uniformly elliptic near zero, has
determinant one, and differs from the Euclidean metric by Schwartz
functions. Its scalar curvature is also Schwartz and depends smoothly
on q in all required norms.

Write $\mathcal D_q$ for the curved divergence of a trace-free covariant
tensor and $L_q$ for the curved longitudinal operator on a contravariant
vector. Tangency to $\det\bar\gamma=1$ gives
$\operatorname{tr}_{\bar\gamma}E_\alpha=0$. Set

$$
\mathcal A_q=\mathcal D_q L_q,\qquad
\mathcal A_q W_\alpha=\mathcal D_q E_\alpha,\qquad
T_\alpha=E_\alpha-L_qW_\alpha .
$$

Here $\mathcal A_q-\mathcal A_0$ tends to zero in bounded-operator
norm from X to Y as $q\to0$. The inverse theorem, or the convergent
Neumann inverse near $\mathcal A_0$, therefore gives a smooth
$q\mapsto W_\alpha\in X^3$. The source is Schwartz, but the resulting
vector need not be Schwartz. Crucially,

$$
T_\alpha\in W^4_{k-1,1}\subset L^2 .
$$

The inclusion has the explicit bound

$$
\|T\|_2^2
\leq \|\langle r\rangle T\|_4^2
\left(\int_{\mathbb R^3}\langle r\rangle^{-4}d^3x\right)^{1/2}.
$$

Raising tensor indices with the smooth bounded metric preserves these
properties. Thus the finite Gram matrix

$$
\mathsf G_{\alpha\beta}(q)=
\int T_{\alpha ij}T_\beta^{ij}\,d^3x
$$

is smooth. At zero, $T_\alpha=e_\alpha$, and linear independence makes
$\mathsf G(0)$ strictly positive. Continuity gives a neighborhood with
positive invertible $\mathsf G(q)$ and smooth inverse. This proves,
rather than assumes, the absence of a Gram degeneracy in a sufficiently
small fixed-rank neighborhood.

The dual-Gram conformal momentum

$$
\sigma^{ij}(q,p)
=\sum_{\alpha,\beta}p_\alpha
(\mathsf G^{-1})^{\alpha\beta}T_\beta^{ij}
$$

then depends smoothly on $(q,p)$ in $W^4_{k-1,1}$ and in L2.
It is linear in p at fixed q. Its squared curved norm
$s=\sigma_{ij}\sigma^{ij}$ belongs smoothly to Y and L1.
The weighted multiplication estimate applies with weights 1 and 1
and target weight 2; the derivative order $k-2$ leaves ample regularity.

## 3. Positive local conformal solution

Set $\psi=1+u$ and define, on the open set $\|u\|_\infty<1/2$,

$$
\mathcal F(q,p,u)=
-8\bar\Delta_q(1+u)
+\bar R(q)(1+u)-s(q,p)(1+u)^{-7}.
$$

This is a smooth map from a neighborhood in finite coefficient space
times X into Y. Multiplication and composition are continuous on the
stated spaces, and the negative power is smooth because the conformal
factor stays uniformly positive. At the origin,

$$
\mathcal F(0,0,0)=0,\qquad
D_u\mathcal F(0,0,0)=-8\Delta:X\longrightarrow Y .
$$

The derivative is an isomorphism. The Banach implicit-function theorem
therefore supplies a unique nearby branch $u(q,p)\in X$ with
$\mathcal F(q,p,u(q,p))=0$. It is smooth in the finite coefficients.
Shrink the coefficient neighborhood so that $1/2<\psi<3/2$ everywhere.
The reconstructed physical metric $\psi^4\bar\gamma$ is then uniformly
equivalent to the Euclidean metric and complete.

This is local uniqueness in the selected Banach neighborhood. It is
not a claim about every solution for arbitrary free coefficients.
The maximal asymptotically Euclidean Yamabe restriction has not been
discarded: the implicit-function solution only constructs a small
admissible neighborhood of flat data. In particular, the p=0 solution
already provides its positive scalar-flat conformal representative.

## 4. Smooth volume energy and equality to the ADM flux

In the project's normalization $16\pi G=1$, define

$$
H(q,p)=\int_{\mathbb R^3}
\left[s(q,p)\psi(q,p)^{-7}-\bar R(q)\psi(q,p)\right]d^3x .
$$

The squared momentum is L1, the curvature is Schwartz, and the conformal
factor and its inverse are uniformly bounded on a smaller parameter
neighborhood. The maps used to form this integrand are smooth into L1.
Consequently H is smooth, with finite-parameter derivatives obtained by
differentiating the integrand in L1. This statement does not rely on
differentiating a nonuniform asymptotic expansion.

Unimodularity gives
$\bar\Delta\psi=\partial_i(\bar\gamma^{ij}\partial_j\psi)$.
Integrating the exact constraint over a ball therefore gives

$$
-8\oint_{S_R}\bar\gamma^{ij}\partial_j\psi\,n_i\,dS
=\int_{B_R}(s\psi^{-7}-\bar R\psi)\,d^3x .
$$

The right side converges by absolute integrability. The physical ADM
surface integrand satisfies

$$
\oint_{S_R}(\partial_j\gamma_{ij}-\partial_i\gamma_{jj})n_i\,dS
=-8\oint_{S_R}\psi^3\partial_r\psi\,dS+o(1).
$$

All omitted terms contain the rapidly decaying conformal-metric
perturbation or its derivative. The difference between this expression
and the constraint flux is bounded by a constant times

$$
R^2\,O(R^{-\beta})\,O(R^{-\beta-1})
=O(R^{1-2\beta})\longrightarrow0 .
$$

Thus the ADM flux exists and equals H. A separate pointwise $1/r$
expansion of u is not needed for this conclusion.

## 5. The local energy minimum

The flat TT conditions give zero first-order scalar curvature, and
the momentum energy begins quadratically. Hence $H(0)=0$ and
$dH(0)=0$. The already derived quadratic functional is

$$
H_2(q,p)=p^T\mathsf G(0)^{-1}p
+\frac14\int|\nabla Q(q)|^2d^3x .
$$

Both terms are strictly positive for nonzero corresponding coefficients:
the Gram matrix is positive; a Schwartz tensor with zero spatial
gradient must vanish; and the finite profile family is linearly
independent. Therefore the full phase-space Hessian $H''(0)$ is
positive definite. Its normalization in the ten-mode example is
recorded in [the verified interaction report](reduced-vertices.md).

Continuity of the Hessian makes H strictly convex on a sufficiently
small ball. It follows by integrating the Hessian along line segments
that H is a positive quadratic form with smooth, position-dependent
coefficients on that ball. This is the input to the sum-of-squares
quantum construction; it is not a statement that a truncated quartic
Hamiltonian is bounded below globally.

## 6. Limits that remain open

Every neighborhood and inverse bound above can depend on the chosen
angular and radial basis. No uniform-in-rank estimate has been obtained.
The proof concerns exact dependent solves on R3, not finite auxiliary
grids, angular projections of nonlinear products, or an occupation-budget
compression. Numerical dependent-field approximations still need their
own convergence estimates.

The resulting constrained geometries are legitimate local vacuum data,
and the local Hamiltonian and symplectic chart are smooth. The finite
image need not be invariant under the full Einstein evolution. Neither
continuum scattering convergence, microscopic corner counting, nor a
unique global quantum completion follows from local existence.
