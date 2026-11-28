# Exact targets for the independent spherical Green reference

These are analytic, independently checked fixtures for the proposed
[unreduced sixth-order reference](unreduced-sixth-order-reference.md).
They have not been evaluated by a new numerical implementation. They define
whole-space inverse and boundary targets without calling the production
HEAT, Fourier-bilinear or tail backends. No mixed H6 validation follows from
this design alone.

Write P=(-Delta)^(-1), C=P/8 and Y=Delta_L^(-1), with
Delta_L W=Delta W+(1/3)grad(div W). The following scalar radial coefficient
multiplies one specified angular harmonic. Its normalization cancels from the
radial equation: use the same convention on the source and answer. Retain
all magnetic components in the implementation and its reconstruction tests.
Useful nonaxial choices are Y_1(n)=n_y and Y_2(n)=n_x n_z; neither is asserted
to have unit angular norm.

## Scalar sources with exact exterior data

For a fixed harmonic degree l, p>=2 and p+l>2, set

    f_l(r) = r^l       for 0 <= r <= 1,
             r^(-p)   for r >= 1.

The regular whole-space C response inside the unit ball is

    u_in(r) = r^l [1/2 + 1/(p+l-2)] / [8(2l+1)]
              - r^(l+2) / [16(2l+3)].

Outside it is

    u_out(r) = {r^(-l-1) [1/(2l+3) + J(r)]
                 + r^(2-p)/(p+l-2)} / [8(2l+1)],

where J(r)=(r^(l+3-p)-1)/(l+3-p) unless p=l+3, in which case
J(r)=log(r). Both u and u' match at r=1. The breakpoint is explicit;
interpolation or differentiation must not assume global smoothness there.
Independent PDE comparisons exclude the join itself when the required
derivative is discontinuous.

The radial Green formula proves these expressions by integrating the source
inside and outside the evaluation radius separately. No fitted boundary
coefficient is used. In particular, for R>=1 the boundary datum is

    R u'(R)+(l+1)u(R) = R^(2-p) / [8(p+l-2)].

The condition p+l>2 alone only ensures convergence of the outer harmonic
integral. Requiring p>=2 also excludes growing responses from this backend's
bounded/decaying inverse domain. At p=2, l must be at least one.

| Degree and exterior source | Interior response | Exterior response |
|---|---|---|
| l=0, p=4 | 1/8-r^2/48 | 1/(6r)-1/(16r^2) |
| l=1, p=4 | 5r/144-r^3/80 | [log(r)/24+1/45]/r^2 |
| l=2, p=4 | 3r^2/160-r^4/112 | 1/(32r^2)-3/(140r^3) |
| l=1, p=2 | r/16-r^3/80 | 1/16-1/(80r^2) |

The p=4 cases test the forced r^-2 particular solution and its resonant
logarithm. The l=1,p=2 case tests an allowed bounded nondecaying dipole:
its boundary datum is exactly 1/8, so homogeneous decaying Robin data are
incorrect. By contrast, l=0,p=2 has a logarithmically divergent exterior
Green integral and must be rejected in this inverse domain.

The general formula also supplies higher-degree controls, including l=8 and
the proposed guard band. Cartesian origin values and derivatives must come
from the regular solid-harmonic expansion, not division by a noisy r^l.

## York source with a traced second moment

Let c be a fixed real Cartesian vector and f=c inside the unit ball, zero
outside. Its mass and second moment are

    M = (4*pi/3)c,
    integral y_a y_b f_j d^3y = (4*pi/15) delta_ab c_j.

The exact vector response is

    W_in = -11c/24 + 19r^2 c/120 - x(c dot x)/60,
    W_out = -[7c+n(c dot n)]/(24r)
            +[3n(c dot n)-c]/(120r^3),

where n=x/r. The r^-3 term is genuine traced-moment information from the
biharmonic Kelvin kernel. Scalar harmonic multipoles alone would miss it.
The auxiliary potential in the York factorization is

    chi_in  = (c dot x)(-1/6+r^2/30),
    chi_out = (c dot n)(-1/6+1/(30r^2)).

Its bounded dipole is required even though W decays. Values and first radial
derivatives match at r=1; Delta_L W=c in the interior and zero in the exterior.
The source jumps at the join, so a PDE residual there is not an ordinary
pointwise smoothness test.

This fixture also supplies exact global integral and boundary checks:

    integral f dot W = -22*pi*|c|^2/45,
    integral (L_0 W):(L_0 W) = 44*pi*|c|^2/45,

    integral_{r>R} (L_0 W):(L_0 W)
      = pi*|c|^2 [5/(6R)-1/(45R^3)+1/(150R^5)]
      = -2 integral_{S_R} W dot ((L_0 W)n),  R>=1.

The source pairing, full longitudinal norm and exterior flux must be
computed through their distinct numerical paths. Do not construct one side
from the expected identity and count that as a comparison.

## A localized auxiliary monopole that must survive

Let g=(1-r^2)^2 inside the unit ball and zero outside, and take f=grad g.
Its scalar P response is

    q_in  = 1/6-r^2/6+r^4/10-r^6/42,
    q_out = 8/(105r).

The York factorization then has h=div(Pf)=-g, chi=-q, and

    W = -(3/4)grad q,
    W_in  = x [1/4-3r^2/10+3r^4/28],
    W_out = 2x/(35r^3).

Consequently integral h=-32*pi/105, chi(0)=-1/6, and
partial_j W_i(0)=delta_ij/4. Deleting the full l=0 channel of h would fail
these targets. Only the forbidden algebraic r^-2 monopole must vanish;
a localized monopole and its nonzero integral are allowed.

## Implementation and acceptance boundary

The numerical API, radial representation and work limits are being aligned
before implementation. A later committed configuration must specify its
refinement levels, comparison points, absolute/scaled tolerances, origin and
join handling, and cumulative work/wall limits before any tracked run.

Require separate radial, angular and exterior-radius refinement, with all
other resolution choices fixed in each slice. Include nonaxial real modes,
regular-origin jets, source/PDE checks away from joins, the full exterior
boundary data, logarithmic resonance, the Kelvin traced term, and separately
computed flux/norm targets. Preserve failures and complete rows before the
next case. Agreement with these fixtures validates inverse machinery only;
the gravitational W2/u3/u4 assembly, independent mixed H6 comparison and
four-point observable remain separate requirements.
