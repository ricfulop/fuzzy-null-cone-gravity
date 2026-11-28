# Noncompact inverse and actual first-response protocol

**Acceptance pending.** This protocol follows the passing
[independent inverse-primitive result](independent-inverse-results.md).
The [initial attempt](noncompact-inverse-memory-repair.md) stopped at its
declared memory ceiling during the final angular slice. Its completed
subchecks are preserved; the repeat changes only the memory ceiling.
Neither this protocol nor a passing first-response check
establishes the independent mixed sixth-order functional or four-point loop.

## Question and operator conventions

The finite radial solver's exterior remainder bound does not change its
central Robin datum. A polynomial-Gaussian source therefore needs an explicit
outside-source response. The new adapter adds that response, its regular
origin shift, and its asymptotic homogeneous coefficient without replacing
the coefficient matched from the numerical interior solve.

For $f_{\ell m}=r^\ell F_{\ell m}(r^2)$, the scalar inverse is
$\mathcal P=(-\Delta)^{-1}$ or $\mathcal C=\mathcal P/8$. York uses
$V=\mathcal P f$, $h=\operatorname{div}V$, $\chi=\mathcal P h$,
$W=-V-\nabla\chi/4$. The full localized $h_{00}$ is retained. All magnetic
components within the declared band survive; guard modes are checked before
restriction. Exact rational anchors cannot silently become different
binary64 grid anchors. Source moments proportional to $\sqrt{\pi/\beta}$
remain distinguished from matched numerical coefficients.

The Gaussian provider has Decimal arithmetic enclosures for its own moments
and source corrections. Those enclosures do not certify radial discretization,
interpolation, angular projection, the homogeneous coefficient, or an energy.
The full field is the algebraic exterior plus its Gaussian remainder.
Both Kelvin trace levels, including the genuine $r^{-3}$ term, are retained.

## Declared work and acceptance

1. Eighteen pointwise algebra groups check the original unreduced recursion:
   manufactured TT metric derivatives; raw noncommuting tangents; nonzero
   $T_1,L_0W_2,E_c$ and $b_2-\mathcal A_1W_1$; rational expansion of the
   original Lichnerowicz equation; nine independent spatial energy impulses;
   a single separate global Gram-square term; and oriented flux densities.
   The scaled tolerance is $2\times10^{-12}$. These are not integrated flux
   identities or constraint solves.
2. Public Gaussian moment enclosures are compared with independent positive
   Decimal quadrature and its explicit truncation tail. The six $(\beta,R)$
   cases are $(1,0),(1,1/8),(3/2,1),(1,4),(1,8),(1,12)$; moment powers are
   $0,1,2,3,8,20,40,68$. Reference precision is 96 decimal digits, quadrature
   orders 24/48, and relative moment/refinement/tail gates are $10^{-22}$.
   Boundary corrections, derivative channels, exact-anchor rules and work
   rejection controls are checked separately.
3. Corrected finite-grid fields are tested for $\mathcal C e^{-3r^2/2}$,
   $\mathcal Y[(1,2,-1)e^{-r^2}]$, and $\mathcal Y\nabla e^{-r^2}$.
   Independent positive quadrature supplies analytic Cartesian targets.
   The gradient case checks $h=-e^{-r^2}$, $\chi=-\mathcal P e^{-r^2}$
   and $W=-3\nabla\mathcal P e^{-r^2}/4$. The grid is $R=4$,
   $N=128,256,512$, fixed origin patch 0.5, and sphere 8x16, with all
   $h$ modes through 2 and a full degree 3 guard band. Structural tests
   apply at every level; physical tests apply at the finest level and the
   latest radial refinement. Value tolerance is $3\times10^{-4}$,
   derivative tolerance $3\times10^{-3}$, and structural tolerance $10^{-10}$.
4. The actual source direction is raw
   $Q=E(\operatorname{diag}(-1/2,-1/2,1))$, $S=B(XZ)$.
   The physical first-source rate is $\beta=1$; geometry's rate2 denotes
   the same envelope. Local polynomials are matched exactly to the archived
   lower-order geometry. Independent mass targets, after dividing by
   $\pi^{3/2}$, are
   $$\pi^{-3/2}\int R_2=-6237/128,\qquad \pi^{-3/2}\int s_2=2079/8,
     \qquad \pi^{-3/2}\int b_1=(0,2079/32,0).$$
   The independent route computes $u=\mathcal C(s_2-R_2)$,
   $v=\mathcal C(R_2+7s_2)$ and $W_1=\mathcal Yb_1$.
   Existing analytic first-HEAT fields supply comparison targets only.
   They supply no source, origin, boundary or grid-solution coefficient.
   Values, gradients and Hessians are compared at origin, both patch sides,
   mixed spatial points, both boundary sides, and exterior points.

The actual-source schedule is $R=4$, $N=256,512,1024$, sphere 8x16, with
an isolated 10x20 sphere refinement at the final $N$. Scalar source support
is $\ell\le4$, vector source support $\ell\le5$, $h$ support $\ell\le6$,
and the full guard band is $\ell\le7$. Exact Fischer reconstruction rejects
insufficient source support. The same value/derivative/structural thresholds
as the synthetic field suite are declared before execution. Passing is not
assumed; underresolution is a failed gate. Coarse discrepancies remain visible.

For this actual-source schedule alone, the adapter plans 435,344 divergence
shell points, 55,724,032 projection products, 435,476 charged point evaluations
and 1,306,296 Cartesian calls. Its FD count is
$$2820\left(\sum_{f\in\{u,v,b_{1x},b_{1y},b_{1z}\}}
 \#\text{source modes}(f)+49\right),$$
computed and checked after exact source projection, before solving. The
independent synthetic suite shares the same cumulative budget. The planned
caps are 500,000 FD unknowns, 1,000,000 point evaluations, 80,000,000 projection
products, 3,000,000 Cartesian calls, 20,000,000 reserved exact products and
5,000,000 stored values. Provider-internal moment and algebra caps are
separate, finite, and recorded; these counters are not hardware FLOPs.
The local invocation has an unchanged 1800-second wall limit. The initial
2-GiB sampled memory ceiling is raised to 8 GiB for the repeat, following a
recorded 2.027-GiB guard interruption. No grid, tolerance, work counter or
scientific expression changes. Full-run timing remains unmeasured; the
interrupted attempt is not a runtime promise for the completed workload.

## Evidence and remaining scope

The driver checkpoints complete rows, adapter/source diagnostics, prepared
comparison targets, completed or stopped levels, work preflights, exceptions,
and exact source-polynomial records. The final conjunction includes every
required comparison and every suite verdict. Completed coarse failures remain
diagnostic; none weakens a finest-grid or structural gate. Frozen parent
snapshots and result bytes remain unchanged.

No $W_2,u_3,u_4$, canonical spatial integrals, integrated finite-ball identities,
or mixed $H_6$ are evaluated by this protocol. A separate exterior-radius
slice for the actual gravitational source also remains pending. A successful
run supports connecting these fields to the
[unreduced mixed-functional design](unreduced-sixth-order-reference.md), with
its forced $u_4$ tails and physical error budget. It does not establish a
one-loop amplitude, a physical large-rank limit, microscopic corner states,
or the original theory's Euclidean–Lorentzian reconstruction.

Source: [driver](../../verification/quantum_bridge/noncompact_inverse_check.py),
[configuration](../../verification/quantum_bridge/noncompact_inverse_config.json),
[actual first responses](../../verification/quantum_bridge/noncompact_first_response.py).
