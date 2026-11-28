# Exact consolidation of Gaussian Fourier kernels

This is a new exact-algebra implementation for the
[Gaussian Fourier backend](gaussian-fourier-bilinear-backend.md), implemented
only in [gaussian_kernel_consolidation.py](../../verification/quantum_bridge/gaussian_kernel_consolidation.py).
It returns the existing immutable `GaussianBilinearKernel` type. The source
compiler and numerical integrator are unchanged. No numerical evaluation,
consolidation call, source construction, or quadrature was executed while
writing this implementation; its checks await a tracked calculation.

## Contract and normalization

```python
consolidate_kernel(kernel, *,
    max_expansion_products=50_000_000,
    max_terms=200_000,
    max_degree=128,
    max_coefficient_bits=4096,
    progress=None) -> (new_kernel, statistics)
```

The input is retained untouched for provenance and independent comparison.
All operations use exact rational coefficients. The result still represents
the same complete bilinear divided by $\pi^{3/2}$. It is passed directly to
the existing integrator with the same axis/precision grids, tolerances, and
evaluation-work cap. The new rewrite-work cap is separate from both the
source compiler's caps and the numerical evaluator's cap.

## Derivation and sparse construction

For each endpoint block write $A=\beta_L+a_Lt^2$, $B=\beta_R+a_Ru^2$ and
$S=A+B$. A present endpoint has $a>0$, so

$$
t^{2m}=a_L^{-m}(A-\beta_L)^m.
$$

This is an algebraic identity in the original quadrature variable, not a
change of variables. There is no Jacobian. A pure endpoint has fixed width
$\beta$ and no heat weight; its Laurent power is absorbed into the exact
coefficient first.

Expand only the heat binomials and merge their coefficients by $(i,j,k)$.
After this exact cancellation, choose $I=\min i$, $J=\min j$, $K=\max k$.
The remaining kernel is

$$
\frac{A^IB^J}{S^{K+1/2}}P(A,B),\qquad
P=\sum_k P_k(A,B)S^{K-k},
$$

where every exponent in $P_k$ is nonnegative. Starting at the smallest
occupied $k$, build $P$ by the sparse recurrence $P\leftarrow SP+P_k$ as
$k$ increases, including missing layers. Multiplication by $S$ has two
monomial scatters: $A+B$ for two heat endpoints, or $\beta+B$ / $A+\beta$
for one fixed endpoint. This avoids separately expanding a large power of
$S$ for every original Laurent term.

A numerator monomial $p_{rs}A^rB^s$ becomes the existing kernel term with
powers $(I+r,J+s,K)$. Present endpoint heat weights become zero; absent
endpoint weights retain sentinel $-1$. With two pure endpoints, the whole
block reduces directly to one rational coefficient divided by $\sqrt S$.
Exact zero coefficients and zero blocks disappear by algebraic equality,
without a tolerance or pruning rule.

## Equal endpoints and identity checks

The existing kernel class canonically orients equal-endpoint terms using
their heat weights and Laurent powers. Removing heat weights can change
which monomials are swapped. Therefore the returned integrand need not
match the original pointwise at an unsymmetrized pair $(t,u)$.

For equal descriptors, test the symmetrized stripped rational integrand

$$
\tfrac12\sqrt{A+B}\,[F(t,u)+F(u,t)].
$$

For unequal descriptors, compare the stripped rational functions directly
in the fixed canonical endpoint order. Rational $t,u$, rates and widths
make these comparisons exact. Square-grid integrals, or paired rectangular
grids exchanged under $t\leftrightarrow u$, respect the same symmetry.
Individual rectangular-grid values and axis-specific diagnostics can differ
between representations while their target integral remains identical.

## Diagnostic and numerical limits

The existing integrator retains every independent axis refinement and its
final precision refinement. These now act on the **compact Laurent basis**.
Its absolute-term sum, cancellation ratio, and termwise refinement/precision
indicators are not asserted to equal or dominate those of the original
basis. Both representations and their diagnostics must remain identifiable
in the result archive.

Expanding $(A-\beta)^m$ introduces alternating coefficients. Exact merging
can remove large algebraic cancellations, but remaining binomial cancellation
can also worsen numerical conditioning. No numerical pass, smaller support,
or runtime improvement follows from exact equivalence alone. The same
unchanged precision and convergence gates may reject the compact result.
The implementation does not enlarge budgets, relax tolerances, select a
smaller observable, or change the gravitational fields.

## Bounded work and recorded evidence

`statistics['work']['total']` is the deterministic counter limited by
`max_expansion_products`. Subcounts record binomial coefficients generated,
expanded heat products, Laurent scatters, numerator scatters, Horner
scatters, output scatters, and fixed-width products. These are declared
rewrite work units, not processor instructions or a runtime estimate.
All caches are local to one call, so a recorded total $W$ supports
deterministic cap-$W$ success and cap-$(W-1)$ rejection controls.

The effective live-map cap is the minimum of `max_terms` and the input
kernel's existing term cap. It includes concurrent intermediate maps and
already completed output blocks. A separate cache holds at most 4096
binomial coefficients. The degree cap bounds heat powers and shifted
numerator degrees. Numerator/denominator bit lengths are checked, including
a conservative preflight before potentially large rational powers.
Exceeded caps raise; they never delete a nonzero term.

Statistics retain per-block input/output term counts, denominator power and
Laurent shifts, orientation swaps, work totals, dimensional term counts,
and peak live support, degree and coefficient size. Progress callbacks emit
block-start, block-complete and final-complete records. Numerical work must
be preflighted separately by the unchanged integrator after consolidation.

The caller's validation should cover fixed widths, unequal heat endpoints,
equal-endpoint orientation changes, negative Laurent powers and mixed
denominator powers, exact $A+B$ and heat-weight cancellations, input
immutability, idempotence, and deterministic cap rejection. Small genuine
two-parameter integrals should compare raw and compact kernels against
independent scalar/York spatial identities. Production export should retain
both exact coefficient representations and their hashes before integration.
