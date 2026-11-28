# Which radial endpoint terms survive the lapse adjoint?

**Analytical preparation, unexecuted.** For the original Paper I bulk action, the smallest endpoint test is the linearized scalar-curvature constraint and its formal adjoint on one spacelike shell. A Dirichlet lapse test leaves a term proportional to the perturbation of the endpoint area element. Even when that term vanishes, the flat Hessian identity used in lapse stabilization has a second endpoint term proportional to the mean curvature of the sphere. The two tests below distinguish these requirements before choosing a domain for the full generator algebra. The source revision is `773ccc4e3e8a1a89e06086163042fbf45e171c76`. [Definitions and source hashes](provenance.json)

## What each domain supplies

1. **Compact bulk tests.** Reconstruct the bandlimited metric, form the ordinary Einstein–Hilbert density, and vary all retained coefficient fields. Compact spacetime support eliminates boundary terms. On a spatial shell $I=[r_-,r_+]$, $0<r_-<r_+$, smooth radial tests supported in its interior justify the integrations by parts defining the canonical bulk expressions. The first-order ADM representative uses this convention.
2. **Reduced lapse form.** At the static flat oblate solution, eliminating the nondegenerate shift pairs gives the compact-test bracket $-G(Z_f,Z_g)/2$. The normal projection removes the second radial derivative, so this bilinear form extends continuously to $H^1_0(I)$ coefficient profiles. On the complementary angular channels the appropriately signed form is coercive for sufficiently small nonzero oblateness on a fixed shell. The resulting weak map acts from $H^1_0$ to $H^{-1}$. A general member of $H^1_0$ has a value trace but need not have a normal-derivative trace or an $L^2$ Hessian. The continuous extension transports the value of the reduced form; specifying domains for its constituent differential operators requires additional work. [Reduced lapse evidence](https://github.com/ricfulop/fuzzy-null-cone-gravity/blob/773ccc4e3e8a1a89e06086163042fbf45e171c76/paper1/reports/oblate-lapse-consistency-2026-11-27/lapse-summary.json)
3. **Physical boundary.** A fixed-$r$ world tube of this static background is timelike. An endpoint sphere on its spacelike slice is a test surface for a bulk Green identity. A causal-diamond variational principle would require specified null boundary data, a normalization of its null generators, a null boundary and joint action, and a check of the resulting variational flux and charges. Paper I leaves its boundary quantum algebra and gluing measure open. The manuscript and canonical definitions are indexed in the [source record](provenance.json).

## Green identity for the unreduced lapse secondary

Use the stationary oblate spatial metric $\sigma$, $x=\cos\theta$, $\rho=r^2+\epsilon$, $\Sigma=r^2+\epsilon x^2$, with $\rho,\Sigma>0$. Let $U=T h$ be a retained spatial coefficient variation pulled back to a covariant metric perturbation. All traces and raised indices below use $\sigma$, except that $\operatorname{tr}_{\gamma} U_{AB}$ uses the induced endpoint metric. Write $dV$ and $dA$ for the physical volume and area measures and $n$ for the outward unit normal on each component of $\partial I\times S^2$. Use real descriptors. For the normalized action,
$$
G(A,U)=\frac1{4\pi}\int(A_{ij}U^{ij}-\operatorname{tr}A\,\operatorname{tr}U)\,dV.
$$
The test $f$ varies the actual coefficient $a=h_{tt}$. At $a=-1$, vanishing shift and unit lapse, $\partial\alpha/\partial a=-1/2$. With Paper I's Hamiltonian curvature sign, the original flat-background lapse pairing is $+(8\pi)^{-1}\int f\,\delta R[U]dV$. Set
$$
P^{ij}(U)=U^{ij}-\sigma^{ij}\operatorname{tr}U,\qquad
\delta R[U]=\nabla_i\nabla_jP^{ij}(U).
$$
Two integrations by parts give the precise defect between the secondary and the Hessian test:
$$
\frac1{8\pi}\int f\,\delta R[U]dV-\frac12G(\nabla^2 f,U)
 =\frac1{8\pi}\int_{\partial(I\times S^2)}
 n_i\left(f\nabla_jP^{ij}(U)-P^{ij}(U)\nabla_j f\right)dA. \tag{1}
$$
At $f|_{\partial I}=0$, all tangential derivatives of its smooth trace vanish and $P^{nn}=-\operatorname{tr}_\gamma U_{AB}$. Thus (1) reduces to
$$
\frac1{8\pi}\int_{\partial(I\times S^2)}(\partial_n f)\operatorname{tr}_\gamma U_{AB}\,dA
=\frac1{8\pi}\left[\rho\int_{S^2}(\partial_r f)\,\gamma^{AB}U_{AB}\,dx\,d\phi\right]_{r_-}^{r_+}. \tag{2}
$$
Here $n^r=\pm\sqrt{\rho/\Sigma}$ and $dA=\sqrt{\rho\Sigma}\,dx\,d\phi$, with $x$ integrated from $-1$ to $1$. The metric trace in (2) is the variation of twice the logarithm of the endpoint area density. It involves both $\delta X$ and $\delta C_{AB}$ on an oblate sphere, even though $\delta C$ is trace-free for the reference round metric.

For arbitrary retained Dirichlet lapse tests with independent endpoint slopes, (2) vanishes for a given smooth $h$ exactly when the **retained scalar moments** of $\gamma^{AB}U_{AB}$ vanish at each endpoint. This follows by choosing the slope of each harmonic independently at just one end. Pointwise fixed area density is a sufficient stronger condition. Another sufficient choice for this pair is to impose both $f=0$ and $\partial_r f=0$ while allowing arbitrary smooth $h$. These are candidate trace restrictions for (1); they specify neither a maximal adjoint domain nor invariant evolution.

## A second boundary form in the compact-test reduction

On this flat spatial metric, $P^{ij}(\nabla^2g)=\nabla^i\nabla^j g-\sigma^{ij}\Delta g$ is divergence-free. Integrating just once gives
$$
G(\nabla^2f,\nabla^2g)
=\frac1{4\pi}\int_{\partial(I\times S^2)}n_iP^{ij}(\nabla^2g)\nabla_j f\,dA. \tag{3}
$$
For smooth $f,g$ with zero endpoint values, the tangential trace of $\nabla^2g$ is $H\partial_ng$, with $H=\operatorname{div}_\sigma n$. Consequently (3) equals
$$
-\frac1{4\pi}\int_{\partial(I\times S^2)} H(\partial_n f)(\partial_n g)\,dA,\qquad
 H=n^r\partial_r\log\sqrt{\rho\Sigma}
  =n^r r(1/\rho+1/\Sigma). \tag{4}
$$
For the outer sphere $H>0$. Taking $f=g$ with independent nonzero outer slope and zero inner slope makes (4) strictly negative. This explains why extending the reduced bracket to $H^1_0$ does not establish $G(\nabla^2f,\nabla^2g)=0$ for smooth Dirichlet profiles evaluated by the original Hessian expression. Imposing both endpoint values and slopes zero suffices for (3). Equations (1)–(4) concern smooth traces; any use of these conditions on Sobolev completions must state the required regularity.

## Smallest proposed check and decision criterion

Use $N=3$, the real $m=0$, $\ell=2$ scalar $Y=P_2(x)$ and the actual retained $X_2$ metric variation. These exist on the exact static oblate background; a small-deformation assumption is needed only when invoking the previously proved shift reduction and reduced coercivity. Work at both ends of a shell with independent jets. As an outer-end test take $F(r)=(r-r_-)^2(r-r_+)$, $f=FY$, and $U_{AB}=Yq_{AB}$ with the other components zero ($q_{AB}$ is the unit-sphere metric). Then (2) has positive outer value since $F'(r_+)>0$ and $\gamma^{AB}q_{AB}>0$. The inner value is zero. With $g=f$, (4) has negative outer value. These are analytic predictions for two sign-sensitive controls.

1. In coefficient variables independently form the strong curvature pairing minus its formal-adjoint volume pairing, with the angular measure and single density factor above. Compare its endpoint coefficients with (1), preserving both radial endpoints and all derivatives. Compare the two-Hessian pairing with (3) using the same shell and tests. Keep the distinct factors of $1/2$.
2. Repeat with inner-end slopes only (the orientation reverses), then with both value and slope zero at both ends (both defects vanish). With Dirichlet slopes free and the retained area moments fixed, only the first defect is required to vanish. Keep a nonzero endpoint control; compact test functions alone cannot distinguish these trace choices. Earlier shared checks of the original action supplied sampled local Green flux comparisons. The proposed test instead pairs the actual endpoint traces and tests a candidate radial domain. [Shared weak-form evidence](https://github.com/ricfulop/fuzzy-null-cone-gravity/blob/773ccc4e3e8a1a89e06086163042fbf45e171c76/paper1/reports/shared-oblate-weak-validation-2026-11-28/evidence-reference.json)
3. A pass would identify trace conditions cancelling these two boundary forms and expose any sign or projection error in the unreduced lapse block. Before declaring a common domain, one must specify spaces and dual pairings for metric, momentum and descriptors; check $J_R,J_R^*$ (the retained shift map contains radial derivatives), $V''$, and compositions including $K_{\rm eff}L^*$; and show the chosen domains are preserved by Hamiltonian evolution and stabilization. A boundary action or boundary charge would need its own variation if it is used to cancel either flux. The endpoint calculations alone neither select physical null data nor yield a Dirac count.

No checker, domain realization, numerical experiment or scientific import was performed for this preparation.

## Preparation review and source access

Independent read-only review found no remaining mathematical or scope blocker (msg_cedbd066). Source links pin the reviewed repository revision; the local [validation record](validation.json) records source hashes, review scope and packaging. No proposed scientific calculation has been executed, and neither report establishes a common physical domain or a final constraint classification.
