# Original third-response sources: preparation and error contract

**Implementation preparation; no numerical acceptance yet.** The
[passed first-response calculation](noncompact-first-response-results.md)
supplies the tested independent route for $u,v,W_1$. The next calculation
must construct actual $F_3,r$ and solve $u_3=\mathcal CF_3$, $W_2=\mathcal Yr$.
This note specifies that connection to the
[unreduced sixth-order reference](unreduced-sixth-order-reference.md).
It does not supply a measured second response or full mixed energy.

## Central sources and independent point checks

The supplied exact local geometry retains its physical Gaussian rates. With
$T=L_0W_1$, the original recursion gives
$$s_3=2\operatorname{tr}\{S[(QS+SQ)/2-T-QS]\},$$
$$F_3=s_3-R_3+8\,\partial_i[(g^{-1})_1^{ij}\partial_j u],
\qquad r=b_2-\mathcal A_1W_1.$$
The point evaluator consumes the computed physical field jets and supplied
local geometry; it constructs no inverse and imports no production reduced
source or bilinear backend. It also exposes the fourth-order sources and the
original nine spatial densities, with the global Gram square excluded.

Inline wiring controls use direct Decimal polynomial/product-rule jets,
metric-derived connections, and a second York construction that first lowers
the vector and then forms covariant derivatives. Artificial nonsymmetric
first vector derivatives and symmetric derivative-pair Hessians expose index
transpositions. Separate source/energy wiring checks, direct $E_c,u_4$
impulses, momentum reversal and joint degree scaling supplement the previously
measured original-algebra controls. They do not constitute inverse or
integrated-functional validation.

## Origin data without fitting

Store finite radial arrays multiplying $H_{\ell m}(x)$, rather than converting
every radial power to a high-degree Cartesian polynomial. Only angular
products need exact Fischer decomposition. For a first PG source
$$f_{\ell m}(r)=r^\ell\sum_k c_k r^{2k}e^{-\beta r^2},$$
the independent exact anchors are
$$B_0=\sum_k\frac{c_k k!}{2\beta^{k+1}},\qquad
C_0=\sum_k\frac{c_k(k+1)!}{2\beta^{k+2}}.$$
For $V=\mathcal Pf$, its regular homogeneous origin coefficient is
$a_0=B_0/(2\ell+1)$; $\mathcal C$ adds $1/8$.

These anchors are compared with the **full computed** coefficient, including
the exterior-provider shift. They supply an error estimate, never a replacement
central coefficient. The second York potential requires the full origin
series, not just an independently corrected constant. Put $d=2\ell+1$ and
write $f=H_{\ell m}\sum_k d_k s^k$, $s=r^2$, for its full PG Taylor series.
The independent exact comparison series uses a formal biharmonic primitive
$\Phi=H_{\ell m}J(s)$ with
$$J(s)=\frac{C_0}{2d(2\ell-1)}
       -\frac{B_0s}{2d(2\ell+3)}
       +\sum_{k\ge0}\frac{d_ks^{k+2}}
       {4(k+1)(k+2)(2\ell+2k+3)(2\ell+2k+5)}.$$
For vector source component $j$, $\chi=\partial_j\Phi$. Differentiating this
full series retains the effect of the computed $V$ homogeneous coefficient
on the particular part of $\chi$. The irrelevant constant of $\Phi$ at
$\ell=0$ drops out under differentiation.

The source builder retains computed origin coefficients and explicit Taylor
remainder/conversion bounds. Inside the declared origin patch it may use the
regular source expression derived from those same computed coefficients,
after checkpointing its comparison with raw pointwise angular projections.
This avoids dividing projection noise by $r^\ell$; it is not a fitted source.
Right-hand patch traces remain explicit. Intermediate-radius FD and
interpolation errors remain separate refinement diagnostics.

## Localized exterior source moments

Outside the matching radius, computed first fields contain their complete
matched Laurent coefficients plus explicit Gaussian remainders. The third
sources are linear in first responses and remain Gaussian-localized; their
positive envelope rates include $1/2$ and $3/2$. There is no bare algebraic
third-source tail. No Green-times-Green product belongs to this linear source
representation; such products are rejected rather than truncated.

For each source harmonic, four signed outer moments suffice:
$$B(R)=\int_R^\infty r^{1-\ell}f_{\ell m}(r)\,dr,\quad
C(R)=\int_R^\infty r^{3-\ell}f_{\ell m}(r)\,dr,$$
$$T(R)=\int_R^\infty r^{\ell+2}f_{\ell m}(r)\,dr,\quad
T_1(R)=\int_R^\infty r^{\ell+4}f_{\ell m}(r)\,dr.$$
The existing scalar Green and derivative-channel identities then supply the
outside-source response and the auxiliary York correction without nested
inverse quadrature. The extra homogeneous contribution from the lower
integration endpoint must be included exactly once.

The new companion uses composite midpoint values with a declared bound,
not quadrature differences presented as an enclosure. For
$g(r)=r^kF(r^2)$, $k\in\{1,3,2\ell+2,2\ell+4\}$,
$$g''=k(k-1)r^{k-2}F+(4k+2)r^kF_s+4r^{k+2}F_{ss}.$$
Positive Gaussian--Laurent envelopes therefore bound each panel error by
$h^3\sup|g''|/24$. Add weighted point-arithmetic errors, directed arithmetic
and conversion errors, and an explicit Gaussian tail beyond the cutoff.
Finite work caps and a missed tolerance cause rejection, not coefficient
pruning, source omission or tolerance relaxation.

These enclosures concern the computed first-field surrogate. Inherited
first-field coefficient errors are recorded separately. For the exact PG
first sources, their exterior remainders are known analytically, so comparison
of the complete computed Laurent coefficients with independent scalar/Kelvin
moments bounds their exterior coefficient error, including both Kelvin trace
levels. Those exact moments do not replace the computed central coefficients.

## Transfer to an eventual physical energy budget

Holding first fields and their canonical correction fixed, second-response
errors obey
$$\delta s_4=-2S:L_0\delta W_2,\qquad
\delta F_4=-2S:L_0\delta W_2+8D_1\delta u_3,\qquad
\delta u_4=\mathcal C\delta F_4.$$
Let $E_W,E_L,E_3,E_H$ bound the absolute component errors of
$W_2,L_0W_2,u_3,\nabla^2u_3$, and let $b_+\ge|B_2|$ be a positive radial
envelope. With $v_+=\mathcal Cb_+$, the original expression gives
$$|\delta H_6|\le\int\left[
2|r|\cdot E_W+14|u|\,|S|:E_L+|B_3|E_3
+v_+\bigl(2|S|:E_L+8|Q|:E_H\bigr)\right]d^3x.$$
Here
$$v_+(r)=\frac18\left[r^{-1}\int_0^r s^2b_+(s)\,ds
                              +\int_r^\infty s b_+(s)\,ds\right].$$
The bounds must cover the integration domain and its exterior. This is an
error-transfer inequality only: the central $u_4$ and original energy remain
explicit. First-field/core FD errors, quadrature uncertainty, canonical
integration and the forced algebraic fourth-source tail need their own
contributions before any common physical-unit acceptance budget can pass.

The preparation APIs are capped and preserve all magnetic modes, origin and
exterior anchors, raw samples, before-restriction support checks and completed
work records. A future tracked child must first measure actual support/work
and then apply declared source, PDE, moment, flux and refinement gates. No
timing, convergence, second-response acceptance or full-$H_6$ result is implied
by this design.
