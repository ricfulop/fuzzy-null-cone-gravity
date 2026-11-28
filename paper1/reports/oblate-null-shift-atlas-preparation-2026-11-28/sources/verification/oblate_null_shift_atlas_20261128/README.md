# Discarded angular coefficient maps for the three oblate null shifts

**Preparation only.** The proposed checker generates three complete linear maps at $N=3$ for the actual Paper I coefficient action. Their inputs are arbitrary retained coefficient variations, including temporal ones. The outputs are the discarded coefficients of $DF_{\xi_i}(c_0)c_1$, per unit radial profile $w_i(r)$, for $\xi_i=w_i n_i K_z$. No scientific import, execution or admissible-seed search was performed in preparing this source.

## Coordinates and labels

Use angular coordinates $(x,\phi)$ with $x=\cos\theta$, $d=1-x^2$, and round metric $q=\mathrm{diag}(1/d,d)$. The Cartesian functions are $(n_x,n_y,n_z)=(\sqrt d\cos\phi,\sqrt d\sin\phi,x)$. The rotations have contravariant components $K_i=(-\partial_\phi n_i,\partial_x n_i)$, giving $K_z=\partial_\phi$. The exact oblate chart uses $r>0$, $\rho=r^2+\epsilon>0$, $\Sigma=r^2+\epsilon x^2>0$, $\lambda_0=\Sigma/\rho$, $\gamma_0=\mathrm{diag}(\Sigma/d,\rho d)$, $a_0=-1$ and $b_0=V_0=W_0=0$. Coefficient and metric reconstruction checks leave $r,\epsilon$ symbolic.

A label is the triple **(field, ell, m)**. Unit scalar amplitudes multiply
$$Y_{\ell m}=P_\ell^{|m|}(x)e^{im\phi},$$
with the Condon–Shortley sign and no spherical-harmonic normalization factor. At negative $m$ we use the conjugate of the positive mode. In particular there is no negative-order factorial conversion. A real field has conjugate coefficients for opposite $m$. The signed complex basis spans the complexification of the stated real space.

For scalar amplitude $y=P_\ell^{|m|}$, suppress the common exponential. Unit contravariant polar and axial vectors are $(d y',im y/d)$ and $(im y,-y')$. The covariant round tracefree polar tensor $E$ has
$$ E_{xx}=y''-xy'/d+\ell(\ell+1)y/(2d),\quad E_{x\phi}=im(y'+xy/d),\quad E_{\phi\phi}=-d^2E_{xx}. $$
Its axial partner has independent entries $(E_{x\phi}/d,-dE_{xx},-dE_{x\phi})$ in the order $xx,x\phi,\phi\phi$. With round area average $dx\,d\phi/(4\pi)$ the scalar squared norm is $s_{\ell m}=(\ell+|m|)!/[(2\ell+1)(\ell-|m|)!]$. Vector and tensor squared norms are respectively $\ell(\ell+1)s_{\ell m}$ and $(\ell-1)\ell(\ell+1)(\ell+2)s_{\ell m}/2$. These conventions reproduce the positive-$m$ basis in [all_modes.py](../oblate_constraint_rank_20261127/all_modes.py).

| Field strings, in catalog order within each $(\ell,m)$ | Input degrees | Input columns |
|---|---|---:|
| `a`, `b`, `lambda`, `X` | $0,1,2$ | 36 |
| `V_polar`, `V_axial`, `W_polar`, `W_axial` | $1,2$ | 32 |
| `C_polar`, `C_axial` | $2$ | 10 |

The catalog orders degrees increasingly, then $m=-\ell,\ldots,+\ell$, then fields in the table order. All ten fields occur at discarded degrees 3 and 4, giving **160 possible rows per generator and 78 input columns**. Rows with $|m|=4$ are included as structural zeros. The axial descriptor preserves input $m$; the transverse descriptors change it by $+1$ or $-1$. Each $n_iK_z$ lies in the vector band through degree two. The equivariant Lie product with a retained field thus has degree at most four. The checker separately tests that reconstructing through degree four reproduces the full coefficient field. This degree bound belongs to the test map only.

## Coefficient map and independent identities

For time-independent descriptors use
$$
DF_\xi(c_0)c_1=\big(\xi(a_1),\xi(b_1),\xi(\lambda_1),
\mathcal L_\xi\gamma_1,[\xi,V_1],[\xi,W_1]\big),\qquad \gamma_1=X_1q+C_1 .
$$
In particular $DF_\xi X=\tfrac12q^{AB}(\mathcal L_\xi\gamma_1)_{AB}$ contains $\xi(X_1)+X_1\operatorname{div}_q\xi$ and the half-trace of $\mathcal L_\xi C_1$. Taking a fixed round trace and then subtracting it supplies the $X$ and tracefree $C$ outputs. Applying scalar transport alone to $X$ would miss a measured axisymmetric coupling. Use the positive round coefficient pairings for projection. The deformed DeWitt metric occurs in a separate canonical calculation.

[atlas_core.py](atlas_core.py) represents mode $m$ by its rational coefficient after removing $d^{p/2}e^{im\phi}$, $p=|m|\bmod 2$. Thus an $x$ derivative acts on the stored coefficient by $\partial_x-px/d$, and a product of two odd modes contributes one factor $d$. The round Hermitian pairing restores $d^p$ before integration. The integrand must cancel to a polynomial in $x$; the checker integrates its even monomials by $\tfrac12\int_{-1}^1 x^{2k}dx=1/(2k+1)$. It reconstructs both retained and discarded components in the angular chart and requires the difference from the directly computed Lie action to vanish identically. Only the discarded coefficients are exported.

At the background $F_\xi V=\partial_r\xi$ and $F_\xi W=0$. The checker verifies the complete affine generator is retained, including its radial-shift component. Varying the fields at fixed descriptor eliminates this affine term from $DF$; it still enters the differentiated metric identity. Write $A=(W,V)$, $h=\left(\begin{smallmatrix}a&b\\b&\lambda\end{smallmatrix}\right)$ and
$$ R(c)=\begin{pmatrix}h+A^t\gamma A&A^t\gamma\\\gamma A&\gamma\end{pmatrix}. $$
For every input and each Fourier branch of all six descriptors, test
$$ DR(c_0)DF_\xi c_1+D^2R(c_0)[c_1,F_\xi(c_0)]=\mathcal L_\xi(DR(c_0)c_1). $$
The right side is a four-coordinate tensor Lie derivative with a separate radial Jacobian column proportional to $w_i'$. On the left, the second derivative has a base block $A_1^t\gamma_0 A_F+A_F^t\gamma_0 A_1$ and a lower mixed block $\gamma_1A_F+\gamma_FA_1$. Here the subscript $F$ denotes the affine generator at the background. The symbols $w_i,w_i'$ are independent throughout this check. The tensor and coefficient routes share harmonic conventions and Fourier arithmetic. No radial derivative of a perturbation profile is required by the vertical, time-independent descriptor.

## Controls and output contract

[check_atlas.py](check_atlas.py) is prepared with the CLI `--config CONFIG --output OUTPUT`. The JSON config is hashed for provenance; the scientific column and row ranges are fixed in source. A future supervisor supplies a fresh output directory, authenticates the executed checker and dependencies, and controls runtime. No budget estimate has been measured for this full calculation.

- Each column includes all three Cartesian maps and the corresponding three exact rotation controls. Each rotation must give zero discarded coefficients. The independent control module receives copies of the previously computed column dictionaries through `mapped_column(Label(field,ell,m), axis, rotation=False)`. It checks the archived axisymmetric rows, mode selection, phase and reality conventions, and defect controls.
- The eight nonzero $m=0$, $i=z$ output rows in the [archived map](../../reports/oblate-null-map-check-2026-11-27/map-checks.json) are $X_4,C_{{\rm p},4},C_{{\rm a},3},C_{{\rm a},4},V_{{\rm a},3},V_{{\rm a},4},W_{{\rm a},3},W_{{\rm a},4}$. The independent highest-mode scalar anchors for input $(\ell,m)=(2,2)$ and output $(3,3)$ are $-i/5$ for $n_xK_z$ and $-1/5$ for $n_yK_z$. The odd coordinate leg alone supplies discarded test weights $(4/105)(6,-1)$ in the $(X_4,C_{{\rm p},4})$ rows per $\rho w f$; its $V_{{\rm a},2}$ term proportional to $f'$ gives zero.
- The `catalog` phase writes `atlas-columns/catalog.json` with ordered labels, config/source digests and runtime. The `affine_background` phase checks all six affine generators. The `atlas_columns` phase has 78 steps, one per input column. Its raw events time each generator and rotation; each atomically published `column-NNN.json` records the catalog digest, six durations, reconstruction status and sparse rows for axes `x,y,z`. Each entry gives a zero-based catalog row index and exact rational strings `real`, `imag`. Multiply by the physical $w_i(r)$ and coefficient amplitude. An absent row in a published column has coefficient zero. A missing column is untested.
- The `independent_controls` phase runs on the assembled maps. Only `final_export` publishes the complete atlas index, with hashes for the catalog and every column, the controls record and scope flags. Partial checkpoints are write-only diagnostics. The checker does not load or resume them.

The atlas supplies coefficient tests for a future pairing with $E_{\rm full}^{(1)}$. It computes neither that residual nor its contraction on an admissible first-order history. The prepared export would establish no additional nonintegrable tangent, nonlinear gauge direction, constraint rank, common boundary domain or physical degree count. The next dynamical obligations are specified in [stabilization preparation, §5](../../reports/oblate-constraint-stabilization-preparation-2026-11-28/README.md).
