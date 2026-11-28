# Finer scalar grid with unchanged inverse-control gates

This child follows the measured
[initial scalar-resolution failure](independent-inverse-initial-result.md).
The original source and experiment remain frozen. It tests whether the
observed approximately second-order decrease in regular third derivatives
continues to meet the existing 0.002 gate.

Only scalar radial resolution changes: R=4 now uses 1024,2048,4096 intervals;
the R=2,4,6 slice consequently holds h=1/1024. The physical origin patch is
still r=1/2. The scalar solver, interpolation, analytic targets, all seven
fixtures, every acceptance tolerance, all structural gates, and the
200,000-unknown work cap are unchanged. Planned scalar unknowns total
140,331, including the deliberately incorrect exterior-source ablation.

Angular controls, the complete two-stage York schedule, cumulative vector
work cap, 1800-second wall limit and sampled 2-GiB RSS limit are unchanged.
The driver merely pins the new explicit scalar schedule. The shared
1024-interval R=4 cases will be compared with their saved parent records.

Extrapolating the parent trend suggests that 4096 intervals could meet the
regular-derivative gate. This is a testable prediction, not an accepted
result or a certified error estimate. If the run reaches the York stage,
that will be its first evaluation; no vector outcome is presumed.

Independent read-only review located the amplification explicitly. For the
interior exact regular solution G=a+b s, write e1=G_s-b and let delta f'
be the interpolated source-slope error at a grid node. The ODE-based bridge
obeys

    G_ss = -(2l+3)e1/(2r^2),
    G_sss = (2l+3)(2l+5)e1/(4r^4) - alpha*delta f'/(8r^(l+3)).

Here alpha=1/8 for C. At r=1/2 the first G_sss coefficient is 1596 for
l=8 and 2300 for l=10. The recorded first-regular-derivative errors therefore
explain the failed third-derivative gate without invoking a new sign error.
This also explains why small raw field errors alone are insufficient.
The same reviewers checked the byte-identical archived result, log hashes,
and exact planned work count. This review launches no new computation.

The task remains validation of independent inverse primitives for the
future unreduced mixed-H6 reference. Four-point scattering, physical
large-rank control, microscopic corner dynamics, and reconstruction of the
original microscopic theory remain outside the conclusion of this child.
