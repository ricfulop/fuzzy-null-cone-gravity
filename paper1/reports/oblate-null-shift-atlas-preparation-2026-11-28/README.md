# Preparation of the all-magnetic null-shift coefficient atlas

This package prepares one exact coefficient-map calculation for the original
Paper I metric-coefficient action. It supplies no executed atlas, new
admissible solution, Euler-residual pairing, or nonlinear classification.
Scientific runtime and outcomes remain unmeasured at this checkpoint.

At $N=3$, the calculation considers the three real angular descriptors
$\xi_i=w_i(r)n_iK_z$, for $i=x,y,z$, at the exact static oblate background.
It will compute $Q_cDF_{\xi_i}(c_0)c_1$ on all 78 retained coefficient
columns, including lapse, radial shift, and angular temporal shift. The
discarded catalog has 160 labels at degrees three and four. These are
dimensions of coefficient spaces, not physical degrees of freedom.

## Equations and independent checks

The coefficient map uses scalar transport for $a,b,\lambda$, vector Lie
brackets for $V,W$, and the tensor Lie derivative for
$\gamma=Xq_{\rm round}+C$. Taking the fixed round half-trace of the latter
determines the transformed $X$; it must include the trace induced by $C$.
For each signed Fourier branch, exact round projection through degree four
is checked by reconstruction of every coefficient tensor.

A separate four-coordinate metric Lie derivative checks

$$
DR(c_0)DF_\xi(c_0)c_1+D^2R(c_0)[c_1,F_\xi(c_0)]
=\mathcal L_\xi(DR(c_0)c_1).
$$

The affine background generator contains $\partial_r\xi$ in the radial
angular shift. Its profile $w'$ stays independent of $w$ in this identity.
Both the shift-times-shift and angular-metric-times-shift terms of $D^2R$
are retained. The coefficient and metric routes share harmonic conventions
and exact Fourier arithmetic; they are not separate computer algebra systems.

Independent coefficient controls will check:

- Six explicit low-harmonic fixtures, including both signs of magnetic
  indices one and two.
- The previously accepted axisymmetric $z$ map on all 22 axisymmetric
  input columns, with nine nonzero coefficients in eight output rows.
  The two temporal-shift rows and the $X_4\leftarrow C_{{\rm a},2}$
  trace coefficient are included.
- Fourier selection, complex conjugacy, and the Cartesian transverse phase
  relation on all columns; all three rotations must have zero discarded map.
- Rejection of three deliberate defects: a changed axial-tensor-to-trace
  coefficient, a deleted temporal-shift coefficient, and a flipped transverse
  phase.

The checker separately fixes two transverse highest-mode scalar coefficients
and the existing $(4/105)(6,-1)$ odd-coordinate-leg normalization. These
are tests of the coefficient map. They do not evaluate the continuum Euler
residual or its pairing on a new first-order solution.

## Run and evidence contract

The prepared invocation uses pinned Python 3.12 and SymPy 1.14.0, one local
scientific child, one thread, a 900-second wall budget, and a 2 GiB sampled
process-group RSS threshold. The threshold is not an operating-system hard
memory cap. The supervisor checks that the invocation and cleanup reserve
fit before the recorded approval cutoff.

The supervisor authenticates the committed sources before dispatching the
checker from those copied bytes. The older axisymmetric result and its
frozen defining source are separately pinned. The scientific checker imports
only its local atlas core, independent controls, SymPy, and standard-library
modules; it does not import the old reducer or execute an archived calculation.

Each completed column is written atomically after its three Cartesian maps,
three rotation checks, coefficient reconstruction and metric chain checks.
Partial columns are write-only diagnostics. They are neither resumed nor
accepted as a complete atlas. Only a final index can authenticate the catalog
and all 78 completed columns, together with the independent controls.

The supervisor validates the final index, every referenced payload hash,
ordered input/output labels, sparse rational complex entries, source and
runtime identity, controls, and scope. Passing telemetry follows the terminal
marker binding the execution and metrics records. This infrastructure
validation checks the evidence structure; the scientific checker performs
the algebraic identities.

## Remaining scientific obligations

Even a successful invocation would supply only the discarded coefficient
tests needed for further Noether compatibility work. A dynamical conclusion
still needs the full linearized Euler pairing on admissible histories and
the corresponding constraint analysis. This package establishes no new
nonintegrable tangent, nonlinear gauge classification, common radial boundary
domain, final Dirac surface, physical degree count, or boundary quantum theory.

The independently reviewed stationary mixed-tangent obstruction and selected
lapse endpoint certificate keep their existing scopes. The present atlas
extends neither result merely by increasing the angular catalog.

## Validation and frozen sources

The scientific source and supervisor/schema received independent static reviews
with no blockers. The author and root each passed 31 standard-library mocked
runner tests; the recorded suite times are 3.879 and 3.507 seconds. Neither
suite imports the scientific checker. Shell syntax, Python AST, JSON and
view-fragment parsing passed. The seven scientific control groups, including
three deliberate defects, remain unexecuted.

[Preparation provenance](preparation-provenance.json),
[author runner validation](runner-validation.json), and
[root runner validation](root-runner-validation.json) preserve the reviewed
hashes and exact captured test commands/output. The [frozen checker README](sources/verification/oblate_null_shift_atlas_20261128/README.md)
records the full basis and wire conventions; its original relative links are
preserved as source bytes. The [manifest](archive-manifest.json) covers every
payload in this package except itself.
