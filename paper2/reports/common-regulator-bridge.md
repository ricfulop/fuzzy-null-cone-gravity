# Aligning the two papers' finite-rank actions

**Status:** proposed common convention, with exact algebraic checks. It is a repair of the product/evaluation mismatch, not a completed finite-rank gravitational quantization. The manuscripts have not been rewritten to adopt this proposal.

**Starting repository:** `fa9c3a7` on `main`, including Paper I's nonlinear-consistency checkpoint and Paper II v0.2. The relevant definitions are [Paper I's bulk action](../../paper1/sections/nonlinear_action.tex) and [Paper II's Route C](../paper2_einstein_gravity.tex). Historical calculations retain their original meanings.

## 1. Recommendation

Use **one angular Galerkin bulk action, evaluated with complete Gaunt interaction tensors**, for both papers. Keep a finite set of independent angular coefficients, but do not truncate the angular content of products, inverses, determinants, or field redefinitions during action evaluation. Paper I already prescribes this bulk operation order.

Paper II's binary projected product can remain a useful three-field tensor and a diagnostic. It must not be iterated as the definition of every nonlinear vertex. Replace that use with the derivatives of the common action itself. This changes the old Route C prescription; it does not prove that the old prescription was already equivalent.

There is no unique, fully specified finite-N quantum action in Paper II to compare at all orders yet. Its projected-product proposal does not fix all parenthesizations, inverse operations, constraints, boundary terms, and measures. The precise result below is that repeated binary projection is insufficient to reproduce Paper I's nonlinear operations, and that a complete tensor prescription can remove this particular mismatch.

## 2. Where agreement ends: an exact identity

Let L=N-1. Let P retain spherical-harmonic degrees at most L, Q=1-P, and use the ordinary real L2 pairing on the sphere. All independent input fields below are in ran(P).

The binary operation is

$$
f\star_N g=P(fg).
$$

For a retained third field h, orthogonality gives

$$
\langle h,P(fg)\rangle=\langle h,fg\rangle.
$$

Thus scalar cubic contractions agree exactly. This alone says nothing about the next order. For four real fields,

$$
\boxed{\langle fg,hk\rangle
=\langle P(fg),P(hk)\rangle
+\langle Q(fg),Q(hk)\rangle.}
$$

The proof is the orthogonal P+Q decomposition; both cross terms vanish. The last term is missing when the two pair products are truncated prematurely. It is fixed by the input fields, not a new fitted coupling. Its absolute value is bounded by ||Q(fg)|| ||Q(hk)||. It vanishes if either pair has no discarded component, but is not zero for arbitrary retained inputs.

For real orthonormal scalar harmonics Y_i, define C_ij,A=<Y_A,Y_i Y_j>. The full quartic tensor is

$$
T_{ijkl}=\int_{S^2}Y_iY_jY_kY_l\,d\Omega
=\sum_{\ell=0}^{2L}\sum_m C_{ij,\ell m}C_{kl,\ell m}.
$$

The repeated binary-product prescription keeps only the terms with ell<=L. Restoring the remaining allowed channels recovers the full tensor and its pairing/permutation identities. For complex or spin-weighted harmonics, the corresponding conjugations, spin selection rules, and tensor contractions must be included; the displayed scalar formula is not itself the gravitational vertex.

The auxiliary summation indices up to 2L are **not additional independent field amplitudes or propagating states**. They represent the product of existing fields. No extra physical mode is introduced by evaluating that product exactly.

## 3. Exact counterexamples and a gravitational restriction

Use x=cos(theta), u=P_2(x), and normalized sphere average <f>=(1/2) integral[-1,1] f(x) dx. P_j here denotes a Legendre polynomial, whereas P above is the band projector.

### Quartic mismatch at L=2, N=3

$$
u^2=\frac15P_0+\frac27P_2+\frac{18}{35}P_4.
$$

Consequently,

$$
\langle u^4\rangle=\frac3{35},\qquad
\langle[P(u^2)]^2\rangle=\frac{69}{1225},\qquad
\langle[Q(u^2)]^2\rangle=\frac{36}{1225}.
$$

The missing fraction is 12/35 of this particular scalar quartic contraction. This is an exact counterexample to general equality, **not a numerical prediction for a graviton amplitude**.

### Quartic agreement does not guarantee sixth-order agreement

Keep the same u but use L=4, N=5. The quartic contraction now agrees, since u^2 has no degree above 4. However,

$$
Q(u^3)=\frac{18}{77}P_6,
$$

and the sixth-order contraction is

$$
\langle u^6\rangle=\frac{53}{1001},\qquad
\langle 1,(P M_u P)^6 1\rangle=\frac{289}{5929},\qquad
\Delta=\frac{324}{77077}>0,
$$

where M_u multiplies by u. A small tree-level check therefore cannot validate all vertices needed for a one-loop calculation.

### The issue also appears in a known restriction of Paper I's Einstein action

Paper I derives the angular bulk density of g=(1+epsilon u) g_Schwarzschild, with a common dimensional prefactor suppressed:

$$
\mathcal L_{\rm ang}
=\frac32\left\langle\frac{|\nabla(1+\epsilon u)|^2}{1+\epsilon u}\right\rangle.
$$

The family is representable for N>=3, and |epsilon|<1 is sufficient for positivity. Here u=P_2 rather than Paper I's sqrt(5) P_2; this only changes the amplitude normalization.

Define a specific projected comparator: set r=|grad u|^2, project r, and define v by P[(1+epsilon u)v]=1 within ran(P). Its Taylor coefficients are v_k=(-P M_u P)^k 1. This is one explicit binary-product inverse rule; it is **not claimed to be a uniquely defined implementation already present in Paper II**.

At L=2 the angular-density coefficients are:

| Order in epsilon | Full pointwise density | Specified projected comparator |
|---|---:|---:|
| 2 | 9/5 | 9/5 |
| 3 | -9/35 | -9/35 |
| 4 | 9/35 | 531/1225 |

The quartic correction is -216/1225, exactly (3/2)<Q(r),Q(u^2)>. The check starts from Paper I's already-derived conformal density; it is not a new independent four-dimensional Ricci derivation or a complete spin-two vertex comparison.

## 4. The common prescription

Adopt Paper I's actual coefficient-to-metric map, with all ten real metric components retained before any justified gauge reduction:

$$
ds_N^2=h_{ab}dx^a dx^b
+(Xq_{AB}+C_{AB})(dy^A+A_a^A dx^a)(dy^B+A_b^B dx^b).
$$

The independent coefficients of h, X, C, and A have their specified spin-harmonic bands. The metric and its nonlinear functions need not remain within those bands. Use the common bare bulk action

$$
S_N^{\rm bulk}[a]
=\frac1{16\pi G_0}\int d^4x\sqrt{-g_N[a]}
\left(R[g_N[a]]-2\Lambda_0\right).
$$

For an admissible background a_bar, define every perturbative vertex by functional differentiation of this same action:

$$
V^{(n)}_{I_1\dots I_n}
=\left.\frac{\delta^n S_N^{\rm bulk}}
{\delta a_{I_1}\cdots\delta a_{I_n}}\right|_{\bar a}.
$$

**Representation equivalence:** evaluating these derivatives by reconstructing the angular functions and performing the full angular integral, or by expanding all their angular contractions into complete Gaunt tensors, gives the same vertices. For a spherical background, finite-order polynomial coefficients have finite harmonic expansions, and the equality follows from linearity of integration and the complete product expansion. For a general angularly varying background, inverse background fields can require infinite harmonic sums or separately converged angular integration; a finite padding rule alone need not be exact. Both methods must use the same metric map, background, normalization, derivative operators, and boundary convention. This is an identity between two evaluations of one specified action, not a proof equating two different actions.

In a spherical-background expansion, finite-order products require only finitely many auxiliary angular channels. Scalar quartic pair contractions require up to 2L. More general n-field contractions must retain the complete band of each intermediate subset; a safe upper bound for n factors of degree at most L is nL. Spin and derivative contractions require their own appropriate selection rules.

For finite-amplitude inverses and determinants, a finite polynomial-degree bound does not make a fixed quadrature exact. Use converged angular quadrature or a locally convergent expansion on the nondegenerate metric domain. The formal Taylor tower specifies interactions order by order; it is not a proof of convergence at arbitrary amplitudes.

The natural higher products are P(f_1 ... f_n), evaluated before the final projection, or directly their integrated tensors. They are not obtained by repeatedly applying only the binary product P(fg). Their use removes a multiplication ambiguity without claiming that P(fg) has become associative.

## 5. Field coordinates must also match

Truncating a metric component and truncating a nonlinear parametrization of that component are different choices. A scalar illustration is sufficient to expose the issue:

$$
A=1+\epsilon u,\qquad b=\log A
=\epsilon u-\frac12\epsilon^2u^2+O(\epsilon^3).
$$

At L=2, u=P_2, b contains the discarded component -9 epsilon^2 P_4/35. Re-exponentiating its truncated version gives

$$
\exp(Pb)-A=\frac{9}{35}\epsilon^2P_4+O(\epsilon^3)\ne0.
$$

Similarly, projection and inversion do not commute: projecting the pointwise inverse of 1+epsilon u differs from its within-band multiplication-operator inverse by -36 epsilon^3 P_2/245 at this cutoff.

Paper II's exponential shear coordinates may be used as continuum computational coordinates, but they cannot be given an independent cutoff and automatically identified with Paper I's linear X,C coefficient ansatz. Their relation must be implemented with its induced harmonics and, for quantization, its measure/Jacobian. Some continuum gauge transformations also leave the finite ansatz; affine Bondi gauge must not be imposed as an exact finite-N gauge freedom without the constraint analysis.

## 6. What this repair does not supply

1. **A Wilsonian effective action.** Restoring angular channels inside a local product is algebraic evaluation, not integrating out dynamical modes with angular momentum above L. Physical high-mode exchange and loop matching are separate calculations, with propagators, determinants, and generally additional induced operators.
2. **A finite-dimensional Hilbert space.** A finite number of angular coefficients leaves continuous coefficient fields and their occupation numbers. The occupation budget, radial regulation, and corner sectors remain separate requirements. Their compatibility with dynamics and physical constraints must be proved or bounded.
3. **The full finite-N constraint algebra.** Paper I's existing Legendre-transform counterexample still applies: angular restriction of a Lagrangian and restriction of the continuum Hamiltonian are not interchangeable. The canonical formulation must be derived from the action actually adopted.
4. **Boundary quantization or Euclidean reconstruction.** The displayed action is a bulk action. A common null-boundary/joint variational principle and physical measure have not yet been derived. The archived HMC model does not become its Euclidean continuation by changing a multiplication rule.
5. **A microscopic derivation of Einstein's action.** This route takes Einstein--Hilbert as the common parent action. Recovering its classical low-energy behavior is a regulator-consistency test; it is not a derivation of Einstein gravity solely from finite-state counting.

These are limits on the present result, not a theorem that the broader program is impossible.

## 7. Connecting classical Einstein gravity and the quantum theory

The target is one consistent construction:

**Common action and boundary data -> canonical constraints/reduction -> physical quantum measure or Hamiltonian -> classical and quantum observables.**

After a justified quantization, its effective action would have the expansion Gamma_N=S_N+hbar Gamma_N^(1)+..., with the same S_N supplying the classical stationary equations and the loop vertices. This expression is the target, not a path integral already constructed in this repository.

The shortest useful implementation sequence is:

1. **One vertex source.** Implement complete spin-weighted interaction tensors through fourth order, and through sixth where required for one-loop four-point diagrams, from the common metric map. Compare independent pointwise and tensor evaluations at small N. Include deliberately omitted-channel controls. The scalar checks supplied here are the first verification layer, not completion of this step.
2. **One constraint analysis.** Derive primary and secondary constraints from that action, classify them, identify zero modes and boundary data, and count physical degrees of freedom. If finite-N symmetry is broken, determine whether a consistent reduction and controlled low-energy restoration are possible. Do not import the continuum Dirac algebra or ghost determinant by assertion. If the action fails this test, revise the common parent in both papers together.
3. **One quantum prescription.** Derive gauge-fixing determinants or second-class/reduced measures as appropriate to those actual constraints. Define a physical inner product and contour, or a self-adjoint physical Hamiltonian. Apply and test the occupation restriction only on that specified structure; the existing compression identities quantify why it is not innocuous.
4. **One-loop matching.** Compute at least a cut-bearing channel and a rational channel with compatible infrared subtraction. Verify the needed quantum symmetry identities and budget-dependent intermediate-state sums. A tree comparison alone does not establish this.
5. **Controlled physical limits.** Separate numerical angular-quadrature convergence from changes in physical rank, diamond size, and occupation budget. Respect the proposed area/rank cap; an unconstrained N->infinity limit at fixed area would be a different proposal. Test low-energy predictions and classical behavior along an explicitly allowed scaling family.
6. **Corner states and Euclidean observables from the same theory.** Derive the boundary representation, area spectrum/degeneracies, and sewing rule. A simulation intended as evidence about this gravity theory must use a derived measure and observable map; retain the older HMC results as results of their stated auxiliary model.

Matching classical and one-loop Einstein effective theory would establish a substantial low-energy bridge. Full finite-state ultraviolet predictivity, corner entropy, and strong-gravity behavior require further results. They should not be counted as completed by this prescription change.

## 8. Context from original research

Donoghue makes the low-energy target precise: "To any given order in the energy expansion there are only a finite number of parameters, which can then be used in making predictions." A shared Einstein effective action can connect classical behavior to low-energy quantum predictions, while leaving ultraviolet completion open. [General relativity as an effective field theory, p. 3](https://www.alphaxiv.org/pdf/gr-qc/9405057?page=3). His background-field calculation includes gauge-fixing and ghost contributions. [p. 6](https://www.alphaxiv.org/pdf/gr-qc/9405057?page=6).

Bahr and Dittrich emphasize that a discrete action's symmetries and its canonical constraints must match. They explicitly state: "These results do not exclude that discrete models with diffeomorphism symmetry can be constructed." Their results concern Regge models, not a no-go theorem for this angular regulator; they motivate testing the actual common action rather than assuming continuum constraints. [(Broken) Gauge Symmetries and Constraints in Regge Calculus, p. 2](https://www.alphaxiv.org/pdf/0905.1670?page=2).

Lehner et al. specify the null-surface and joint contributions and explain the remaining dependence on null-generator normalization. Their classical additivity prescriptions are relevant input for a common boundary action; they do not themselves derive the finite quantum corner states sought here. [Gravitational action with null boundaries, p. 36](https://www.alphaxiv.org/pdf/1609.00207?page=36).

The decisive evidence for the product repair is the algebra above, rather than a figure from these papers. The literature supplies the quantum and constraint context, not a proof of this proposed finite-N completion.

## 9. Reproduction and present deliverables

From `paper2/`, run:

```bash
python3 verification/common_regulator_identities.py
```

The checker uses only the Python standard library and exact rational arithmetic. It checks Legendre orthogonality, every axisymmetric cubic and quartic basis tuple at L=2 and L=3, all three quartic pairings with complete channels, the explicit fourth- and sixth-order counterexamples, inversion and field-coordinate mismatches, and the stated conformal-action comparison. Results are archived in [common-regulator-identities.json](common-regulator-identities.json); code is in [common_regulator_identities.py](../verification/common_regulator_identities.py).

**Established now:** an exact obstruction to the naive identification; the missing quartic term; a common bulk vertex convention that eliminates the product mismatch when adopted; and a reproducible check of the algebra and one conformal sector.

**Still required:** adoption in the manuscripts and full spin-two vertex implementation using identical field coordinates, finite-N canonical analysis, quantum measure and loop amplitudes, dynamical corner states, and Euclidean--Lorentzian reconstruction.
