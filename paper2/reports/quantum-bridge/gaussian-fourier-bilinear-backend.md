# Gaussian Fourier compilation of nested scalar and York bilinears

This report gives a new analytic derivation and the implementation contract for
[gaussian_fourier_bilinears.py](../../verification/quantum_bridge/gaussian_fourier_bilinears.py).
At creation the implementation had only AST parsing and import preflight.
The subsequent [tracked ingredient suite](sixth-order-backend-validation.md)
passed nine normalization, adjoint, derivative and spectral-bound controls,
including two-dimensional heat-parameter integrals. Production mixed
gravitational bilinears and the full sixth-order Hamiltonian remain
unevaluated; the derivation below is distinct from that bounded validation.

The conventions and the need for these pairings come from the
[sixth-order reduction](sixth-order-bilinear-and-tail-reduction.md),
[original recursion](reduced-sixth-order-and-loop-kernel.md), and
[contracted four-point prescription](contracted-four-point-implementation.md).
The Fourier compilation and numerical design below are new derivations, not
claims sourced from an external calculation. In particular this backend does
not by itself evaluate the sixth-order Hamiltonian or a gravitational loop.

## 1. Source domain and the spatial integrations that disappear

Write $\mathcal C=(-\Delta)^{-1}/8$ and $\mathcal Y=\Delta_L^{-1}$. The remaining
nested terms in the stated sixth-order energy are

$$
-\langle K,\mathcal CF\rangle
-2\langle r+2t,\mathcal Yr\rangle.
$$

The endpoints are finite sums of pure polynomial Gaussians and terms

$$
f(x)=P(x)e^{-\beta |x|^2}I_m(a|x|^2),\qquad
I_m(z)=\int_0^1t^{2m}e^{-zt^2}\,dt,
\quad \beta,a>0,\quad m\in\mathbb N_0.
$$

The [source representation](../../verification/quantum_bridge/gaussian_heat_sources.py)
stores exact real rational coefficients and positive rational rates. It groups
polynomials by $(\beta,a,m)$; the pure Gaussian sentinel is precisely
$(a,m)=(\mathrm{None},\mathrm{None})$. Thus a heat group is

$$
\int_0^1t^{2m}P(x)e^{-A(t)|x|^2}\,dt,
\qquad A(t)=\beta+at^2\in[\beta,\beta+a].                 \tag{1}
$$

For the current sixth-order sources the pure rate is $3/2$, and a heat group
has $A(t)=1/2+t^2$. The implementation permits other strictly positive rates.
Neither endpoint contains two heat factors. Spatial polynomial degree is
bounded by a configurable compilation cap; no angular projection is made.

Every group is localized:

$$
|f(x)|\le\sum_p\frac{|P_p|\,|x|^{|p|}e^{-\beta |x|^2}}{2m+1}.
                                                               \tag{2}
$$

Widths in (1) stay in a compact positive interval. Their Fourier polynomials
therefore have uniform Gaussian decay, with an integrable $|k|^{-2}$
singularity for the scalar inverse. The York longitudinal numerator supplies
two powers of $k$, so its apparent $|k|^{-4}$ denominator has the same local
integrability. These bounds justify exchanging the finite heat-parameter and
spatial/Fourier integrations. The Fourier integration can be done analytically
before either heat integral. There is no spatial cutoff, omitted spatial tail,
inverse heat parameter, or reconstructed second-response field in this method.
The separately required three-HEAT local term with $\beta=0$ is outside its
domain and remains governed by the exterior analysis in the reduction report.

## 2. Fourier polynomial, phases, and the York numerator

Use $\widehat f(k)=\int e^{-ik\cdot x}f(x)\,d^3x$ and inverse measure
$d^3k/(2\pi)^3$. For a polynomial Gaussian,

$$
\widehat{P e^{-A|x|^2}}(k)
=\left(\frac\pi A\right)^{3/2}e^{-|k|^2/(4A)}Q_P(A;k),
$$
$$
Q_P(A;k)=\sum_\alpha P_\alpha i^{|\alpha|}
 e^{|k|^2/(4A)}\partial_k^\alpha e^{-|k|^2/(4A)}.            \tag{3}
$$

Compile the derivative polynomial into rational coefficients of
$k^\gamma A^{-u}$. A Cartesian derivative acts on a monomial by

$$
k^\gamma A^{-u}\longmapsto
\gamma_i k^{\gamma-e_i}A^{-u}
-\tfrac12 k^{\gamma+e_i}A^{-(u+1)}.                        \tag{4}
$$

The factor $i^{|\alpha|}$ places a term in a real or imaginary bin, with sign
$(-1)^{\lfloor|\alpha|/2\rfloor}$. The scalar pairing integrates the product
of the two real bins plus the product of the two imaginary bins. The omitted
cross-bin terms are odd under $k\mapsto-k$ and integrate to zero.

The flat York symbol is

$$
\widehat{\mathcal Y}_{ij}(k)
=-\frac{\delta_{ij}}{|k|^2}+\frac{k_i k_j}{4|k|^4}.          \tag{5}
$$

Its first term uses componentwise scalar contractions with coefficient $-1$;
its second contracts $k\cdot Q_b$ with $k\cdot Q_c$ with coefficient $1/4$.
Multiplication by $k_i$ preserves the **original** real/imaginary phase bin in
(3). Reclassifying the phase by the raised polynomial degree would incorrectly
insert the Fourier factor of a spatial derivative and reverse the longitudinal
sign. Every surviving longitudinal monomial has degree at least two. A
constant numerator multiplying $|k|^{-4}$ is rejected as nonintegrable.

This API accepts real exact rational source coefficients only. For real
sources $\widehat f(-k)=\overline{\widehat f(k)}$, which explains the two-bin
rule. A future complex continuation must implement the bilinear expression
$\widehat f(-k)\widehat g(k)$ and must not conjugate the complex source
coefficients; the physical functional here is symmetric bilinear.

## 3. Exact moments and Laurent kernels

Let $p$ have even Cartesian entries, $|p|=2N$, and
$D_p=\prod_i(p_i-1)!!$, with $(-1)!!=1$. The angular integral and the elementary
radial Gaussian integral give

$$
\int_{\mathbb R^3}\frac{k^p e^{-\lambda |k|^2}}{|k|^{2s}}\,d^3k
=\frac{2\pi D_p}{(2N+1)!!}
\Gamma(N+3/2-s)\lambda^{-N-3/2+s}.                        \tag{6}
$$

Odd Cartesian powers integrate to zero. Here $s=1$ permits $N\ge0$;
$s=2$ requires $N\ge1$. The common Fourier prefactor is
$1/[8(AB)^{3/2}]$ and $\lambda=(A+B)/(4AB)$.

For a product term $c\,k^p A^{-u}B^{-v}$, canceling the half-integer Gamma
factor in (6) gives the following **normalized** spatial moments, with
$\pi^{3/2}/\sqrt{A+B}$ factored out:

$$
M_1=c\,\frac{2^{N-1}D_p}{2N+1}
\frac{A^{N-1-u}B^{N-1-v}}{(A+B)^N},\qquad N\ge0,           \tag{7}
$$
$$
M_2=c\,\frac{2^{N-2}D_p}{(2N+1)(2N-1)}
\frac{A^{N-2-u}B^{N-2-v}}{(A+B)^{N-1}},\qquad N\ge1.       \tag{8}
$$

The compiler applies $M_1/8$ for $\mathcal C$, and $-M_1$ plus $M_2/4$ to
the corresponding York numerators. Consequently a compiled block is a finite
rational Laurent sum with no transcendental function besides $\sqrt{A+B}$.

Direct algebraic normalization controls follow without quadrature. For
$g_A=e^{-A|x|^2}$ and constant vectors $c,d$,

$$
\langle g_A,\mathcal Cg_B\rangle
=\frac{\pi^{3/2}}{16AB\sqrt{A+B}},\qquad
\langle c g_A,\mathcal Y(d g_B)\rangle
=-\frac{11\pi^{3/2}(c\cdot d)}{24AB\sqrt{A+B}}.            \tag{9}
$$

Indeed the $N=0$ scalar inverse moment is $1/(2AB)$, while the $N=1$
longitudinal moment is $1/(6AB)$. Thus the vector constant is
$-1/2+1/24=-11/24$. Odd-parity source tests are also necessary: constants
alone cannot detect a lost imaginary Fourier bin.

## 4. Immutable representation and exact merging

A `GaussianBilinearKernel` stores the read-only mapping

```text
blocks[((beta_left, a_left_or_None), (beta_right, a_right_or_None))]
      [(m, n, i, j, k)] = Fraction
```

The term means, without $\pi^{3/2}$,

$$
c\int\!dt\int\!du\;
\frac{t^{2m}u^{2n}A(t)^iB(u)^j}{(A(t)+B(u))^{k+1/2}}.      \tag{10}
$$

A pure endpoint has moment sentinel $-1$: its weight is one and its integral
is absent. This is never an actual negative power of a heat variable.

Endpoint descriptors are ordered canonically using reciprocity. Swapping
endpoints also swaps both heat weights $(m,n)$ and powers $(i,j)$, preserving
their association with the rates. When the two descriptors coincide, an
individual term can likewise be exchanged under $t\leftrightarrow u$.
This preserves its integral, although the canonically stored integrand need
not equal the original integrand pointwise. Exact fractions with the same
resulting key merge before quadrature; zeros are deleted. This is a sparse
representation, not a complete normal form for every rational identity among
$A,B,A+B$ or every heat-moment recurrence.

Transforms of repeated polynomials are cached within a compilation. Hermite
coefficients and angular moment factors use bounded immutable caches. Cartesian
parity grouping avoids products whose angular moments vanish. York divergence
groups combine all three vector components before longitudinal contraction.
There is no dependence on numerical first-Green evaluation or on `heat_integrals`.

The public calls are:

```python
compile_scalar_bilinear(left, right, *,
    max_pair_products=20_000_000, max_fourier_terms=200_000,
    max_terms=200_000, max_degree=64, progress=None)
compile_york_bilinear(vector_left, vector_right, *,
    max_pair_products=20_000_000, max_fourier_terms=200_000,
    max_terms=200_000, max_degree=64, progress=None)

kernel.scale(exact_factor)
kernel.add(other, factor=Fraction(1))
kernel.merge(other, factor=Fraction(1))
kernel.statistics()
kernel.integrate(**integration_options)
```

Both compilers return a new immutable kernel. Both addition methods return a
new kernel; neither mutates an operand. Therefore signed Hamiltonian families
can be combined exactly before a single integration. Family kernels may be
retained separately for diagnostics without recompiling their sources.

Compilation caps bound emitted Fourier monomials, parity-eligible pair
products, live output terms, and source degree. The pair-product cap also
independently caps source-group pair traversal, including zero contractions.
These are explicit operation/term counts, not bounds on rational bit complexity
or measured runtime. An exceeded cap raises rather than silently truncating.
No actual $F,K,r,r+2t$ source-size benchmark has yet been performed for this
backend; the defaults have not been increased in anticipation of one.
An optional `progress(record)` callback receives JSON-compatible cumulative
counts on entering an endpoint-block visit and at successful completion.
The block record includes its canonical endpoints and inverse denominator
power; a later revisit can emit the same descriptor again. The completion
record and `kernel.statistics()` report actual compilation counts, including
Fourier terms emitted, pair products, source-group pairs, and peak output terms.
The callback does not alter the computation or suppress exceptions.

## 5. Decimal integration and the tolerance convention

```python
integrate_heat_bilinear(kernel, *,
    orders=(12, 24, 36), precision=64, check_precision=96,
    absolute_tolerance=1e-10, relative_tolerance=1e-9,
    max_term_evaluations=20_000_000, max_power=4096, max_order=256)
```

**All absolute tolerances and authoritative Decimal outputs are in normalized
units without $\pi^{3/2}$.** Relative tolerance is unchanged. A physical
absolute budget $\epsilon$ for a single unscaled functional must therefore
be supplied as $\epsilon/\pi^{3/2}$ in these units; any further family factors
also belong in the caller's budget or in exact kernel scaling. Returned
`settings["ABSOLUTE_TOLERANCE_UNITS"]` explicitly records this convention.
`normalized_value_decimal` is the high-precision result string. The full
`value` includes $\pi^{3/2}$ using binary64 arithmetic and is a display or
binary64 assembly value, not a high-precision representation of $\pi$.
If conversion overflows or rounds a nonzero Decimal to binary64 zero, the
display is `None` and `value_float_out_of_range` is true. The Decimal result
and its acceptance arithmetic remain authoritative and unchanged.

The evaluator uses Decimal Gauss--Legendre quadrature on $[0,1]$. The shared
`run.gauss_nodes` returns **half** of the standard $[-1,1]$ weights. Under
$t=(x+1)/2$ those returned weights are used unchanged. The precision must be
at least 56 digits to satisfy that routine's normalization guard.

For each neighboring pair of requested orders $L<H$, evaluate the four grids
$(L,L),(H,L),(L,H),(H,H)$, reusing grids already present. A pure endpoint
collapses its axis to one evaluation. The final diagonal grid is reevaluated
at `check_precision`, which must exceed `precision`. Widths, integer powers,
square roots, summation, differences, and acceptance arithmetic all use Decimal.
No square-variable substitution or endpoint-singular weight is introduced.

Let $v_j^{HH}$ denote the integral of one retained signed Laurent term, and
define

$$
E_t=\sum_j|v_j^{HH}-v_j^{LH}|,\qquad
E_u=\sum_j|v_j^{HH}-v_j^{HL}|,\qquad
E_d=\sum_j|v_j^{HH}-v_j^{LL}|.                              \tag{11}
$$

Absent axes have zero change. The final parameter indicator is
$E_{\rm par}=\sum_{\rm blocks}\max(E_t+E_u,E_d)$. The final precision
indicator is likewise the sum of absolute changes of individual terms,
$E_{\rm prec}$. The diagnostic `converged` requires

$$
E_{\rm par}+E_{\rm prec}
\le\epsilon_{\rm abs}+\epsilon_{\rm rel}|V|,                \tag{12}
$$

and a reported working-precision rounding proxy no larger than the same budget.
The proxy for a block is
$10^{1-p}\sum_j|v_j|(1+N_{\rm nodes}+N_{\rm terms})$; it is a heuristic,
not an interval-arithmetic bound. Per-grid signed values, termwise changes,
signed changes, precision changes, and cancellation ratios
$\sum_j|v_j|/|\sum_jv_j|$ are all returned as JSON-compatible strings or
scalars. A zero sum with nonzero terms is flagged separately. Positive rates
make each individual Laurent basis function positive, so its coefficient
determines its sign throughout the integration domain.

Termwise indicators prevent a small change in the final cancelling sum from
being the only acceptance test. They still do not bound unresolved quadrature
error, missed narrow structure, or all rounding effects. Consequently every
result sets `accuracy_certified=False`. Intermediate refinement rows are
retained for inspection. A nonconverged result is not an accepted energy.

Before requesting any quadrature nodes, the integrator counts all planned
term-node evaluations, including independent axis grids and the final precision
check. It rejects excessive work, invalid rate domains, excessive powers, or
invalid orders. The exact zero kernel requires no quadrature. General positive
rates can have large condition numbers; this interface does not promise that
the default orders resolve every permitted rate ratio.

## 6. Controls for the future tracked validation

The following are mathematical targets and proposed checks, not measured
results of this task:

1. Verify (9) at unequal positive rates, in both endpoint orders, with scalar
   constants and several nonparallel vector pairs. Check zero cross terms and
   pure/pure zero-dimensional evaluation.
2. For polynomial Gaussians of mixed parity, verify exact reciprocal compiled
   representations, bilinearity, coordinate permutations, and exact merging
   against separately retained families. Include pure/heat and heat/heat
   endpoints with different rates and moment orders so misplaced heat weights
   cannot pass by symmetry alone.
3. Use the Fourier multiplier identities
   $$
   \sum_i\langle\partial_i f,\mathcal C\partial_i g\rangle
   =\tfrac18\int fg,\qquad
   \langle\nabla f,\mathcal Y\nabla g\rangle
   =-\tfrac34\int fg.                                    \tag{13}
   $$
   For pure polynomial Gaussians, the right sides have independent exact
   Gaussian moment formulas. These tests exercise odd Fourier bins and the
   York longitudinal sign. A divergence-free vector source independently
   reduces the York pairing to minus the componentwise scalar inverse.
4. Check scalar self-pairing nonnegativity and the vector self bounds
   $$
   -E\le\langle b,\mathcal Yb\rangle\le-\tfrac34E,
   \qquad E=\sum_i\langle b_i,(-\Delta)^{-1}b_i\rangle,      \tag{14}
   $$
   against the declared numerical error budget. Cross bilinears need not have
   a definite sign. These spectral bounds follow directly from the transverse
   and longitudinal eigenvalues in (5).
5. Compare selected heat-source scalar pairings with an independently
   implemented radial/harmonic Green integral. With real spherical harmonics
   normalized by $\int Y_{\ell m}Y_{\ell'm'}d\Omega
   =\delta_{\ell\ell'}\delta_{mm'}$, its exact scalar radial formula is
   $$
   (\mathcal Cg)_{\ell m}(r)=\frac1{8(2\ell+1)}
   \left[r^{-\ell-1}\int_0^r s^{\ell+2}g_{\ell m}(s)ds
   +r^\ell\int_r^\infty s^{1-\ell}g_{\ell m}(s)ds\right].  \tag{15}
   $$
   Pair with $f_{\ell m}(r)r^2dr$ and sum the harmonics. This check needs its
   own radial resolution and tail accounting; agreement with another caller
   of the same Fourier compiler is not independent validation.
6. Exercise work-cap failures, invalid source domains, zero kernels, strongly
   cancelling sums, independent parameter refinements, and precision
   refinement. Record source/kernel statistics before selecting larger budgets.

The current module supplies a bounded finite-source evaluation mechanism for
the two nested inverse families. It neither changes the reduced canonical
candidate nor settles finite microscopic state counts, local-corner charges,
large-rank limits, physical scattering, or the missing validation of the full
sixth-order energy and loop assembly.
