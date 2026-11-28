# York operator: finite-metric and perturbative implementation

This is analytic implementation guidance, not a numerically verified evaluator.
It specializes the [sixth-order reduction](reduced-sixth-order-and-loop-kernel.md),
Section 3, to inputs already available from metric and Green-function jets.
No calculations were run and no current implementation was changed for this note.

## 1. Index types and required inputs

Let $g_{ij}=(e^Q)_{ij}$ in three Cartesian coordinates, with $\det g=1$.
The unknown is the **contravariant vector** $W^k$. The source and operator output
$b_i=(\mathcal A_gW)_i$ are **covectors**. The pairing is
$\int V^i b_i\,d^3x$, with coordinate volume and no extra metric factor.
When varying $g=e^{\lambda Q}$ to extract $\mathcal A_1,\mathcal A_2$, hold the
contravariant components of the supplied test field fixed. Holding $W_i=g_{ik}W^k$
fixed instead produces a different expansion.

A concrete pointwise interface is
`york_apply(g, dg, ddg, W, dW, ddW) -> {L, dL, A}`.
Use `dg[a][i][j]=partial_a g_ij`, `ddg[a][b][i][j]=partial_a partial_b g_ij`,
`dW[a][k]=partial_a W^k`, and `ddW[a][b][k]=partial_a partial_b W^k`.
Compute $g^{-1}$ and $\Gamma$ from these inputs; no derivative of $\Gamma$ or
third spatial derivative is needed. This formula assumes an exactly unimodular
metric; a general determinant requires additional divergence terms.

## 2. Explicit finite-metric formula

Put $d=\partial_kW^k$ and $d_{,a}=\partial_a\partial_kW^k$. First form
$$
T_{ij}=(L_gW)_{ij}
=W^k g_{ij,k}+g_{ik}W^k{}_{,j}+g_{jk}W^k{}_{,i}
-\frac23g_{ij}d.
$$
Its derivative is a direct product rule:
$$
\begin{aligned}
T_{ij,a}={}&W^k{}_{,a}g_{ij,k}+W^k g_{ij,ka}
+g_{ik,a}W^k{}_{,j}+g_{ik}W^k{}_{,ja}\\
&+g_{jk,a}W^k{}_{,i}+g_{jk}W^k{}_{,ia}
-\frac23(g_{ij,a}d+g_{ij}d_{,a}).
\end{aligned}
$$
Then, with $\Gamma^m{}_{ai}=\tfrac12g^{m\ell}
(g_{\ell i,a}+g_{\ell a,i}-g_{ai,\ell})$,
$$
\boxed{(\mathcal A_gW)_i
=g^{ja}\left(T_{ij,a}-\Gamma^m{}_{ai}T_{mj}
-\Gamma^m{}_{aj}T_{im}\right).}
$$
Both connection terms are required: $T$ has two covariant indices. In flat space
this reduces to $\Delta W_i+\tfrac13\partial_i\partial_kW^k$, using the flat
metric only to identify the vector and covector components in this check.

## 3. Extracting the two perturbative operators

Write $g=\sum_n\lambda^ng_n$, $h=g^{-1}=\sum_n\lambda^nh_n$,
$\Gamma=\sum_{n\ge1}\lambda^n\Gamma_n$, with
$g_0=h_0=I$, $g_1=Q$, $g_2=Q^2/2$, $h_1=-Q$, and $h_2=Q^2/2$.
Define $T_n(W)$ by replacing $g$ by $g_n$ in the expression for $T$ above;
the same replacement gives $\partial_aT_n$.
For any covariant tensor field $T$, define
$$
(\mathcal D_0T)_i=\delta^{ja}\partial_aT_{ij},\qquad
\mathcal B(\Gamma,T)_{aij}=\Gamma^m{}_{ai}T_{mj}+\Gamma^m{}_{aj}T_{im},
$$
$$
\begin{aligned}
(\mathcal D_1T)_i&=h_1^{ja}\partial_aT_{ij}
-\delta^{ja}\mathcal B(\Gamma_1,T)_{aij},\\
(\mathcal D_2T)_i&=h_2^{ja}\partial_aT_{ij}
-h_1^{ja}\mathcal B(\Gamma_1,T)_{aij}
-\delta^{ja}\mathcal B(\Gamma_2,T)_{aij}.
\end{aligned}
$$
These $\mathcal D_n$ are tensor-divergence coefficients; they are
not the scalar-Laplacian perturbations called $D_n$ in the sixth-order report.
The required operator coefficients are
$$
\boxed{\mathcal A_1W=\mathcal D_0T_1(W)+\mathcal D_1T_0(W),\qquad
\mathcal A_2W=\mathcal D_0T_2(W)+\mathcal D_1T_1(W)+\mathcal D_2T_0(W).}
$$
Thus $\Gamma_1,\Gamma_2$ and metric coefficients through two suffice. The output
is an ordinary series coefficient: the second derivative in $\lambda$ is
$2\mathcal A_2W$, not $\mathcal A_2W$.

For Gaussian-polynomial $g_n$, evaluate its physical first/second spatial jets
using Gaussian rate $n$. The Green field $W_1$ is generally non-Gaussian:
use its independently evaluated value/gradient/Hessian in the product rules.
Do not apply a Gaussian polynomial derivative to $W_1$ as though it shared
the metric envelope. In amplitude checks, freeze the reference $W_1$ field;
varying its own amplitude simultaneously mixes different operator orders.

## 4. Use and independent checks for the tracked harness

Apply the strong formula to obtain $r_2=b_2-\mathcal A_1W_1$, then solve
$W_2=\mathcal A_0^{-1}r_2$. The pairings entering $\mathcal N_3,\mathcal N_4$
also have an independent weak-form evaluation:
$$
\int V^i(\mathcal A_nW)_i
=-\frac12\sum_{a+b+c+d=n}\int
h_a^{ik}h_b^{j\ell}[T_c(V)]_{ij}[T_d(W)]_{k\ell}\,d^3x.
$$
This follows by expanding $\langle V,\mathcal A_gW\rangle
=-\tfrac12\int g^{ik}g^{j\ell}(L_gV)_{ij}(L_gW)_{k\ell}$.
For production $\langle W_1,\mathcal A_2W_1\rangle$, it can avoid second
derivatives of the Green field; pointwise $\mathcal A_1W_1$ is still needed.
Boundary terms vanish under the stated decaying-domain hypotheses; numerical
integration must include exterior-tail errors rather than silently dropping them.

Check flat recovery, $g^{ij}T_{ij}=0$, $T_{ij}=T_{ji}$, spatial Hessian symmetry,
and strong/weak pairing agreement for distinct test vectors. Compare $\mathcal A_1W$
and $\mathcal A_2W$ with centered finite-metric evaluations using independent
matrix-exponential spatial jets, multiple amplitudes, and increased precision.
When $W_1$ is counted as total field order two, $\mathcal A_1W_1$ has order three
and $\mathcal A_2W_1$ order four. These checks validate this operator ingredient;
they do not by themselves evaluate H5/H6 or the gravitational loop coefficient.
