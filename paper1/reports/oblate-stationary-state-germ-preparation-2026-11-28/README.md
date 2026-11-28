# Stationary oblate solution germs in the fourteen state coordinates

**Preparation only.** This route would certify a local analytic solution of the stationary linearized Paper I equations by solving the full state ODE and proving that it preserves the six constraints. It applies to the zero-remaining-control, gauge-fixed, axisymmetric reflection-even sector at $N=3$ and $\epsilon=1/100$. The physical-radius minor, regular interval and responding initial datum described below have yet to be certified. The separate 1800-second retry retains the reduced-coordinate calculation. [Runtime-only retry preparation](../oblate-stationary-extended-budget-preparation-2026-11-28/README.md)

## Data and physical radius

Use independently pinned copies of the original stationary operator export and the intermediate state snapshot. Verify their provenance, coefficient order and derivative convention. Write $z=(v,w)$, with seven profiles ordered as $(\lambda_0,\lambda_2,X_2,V_{\rm polar,2},C_{\rm polar,2},a_0,a_2)$. The proposed ODE is $z'=Kz$. Its data consist of the original $6\times14$ constraint matrix $C$, the $14\times14$ matrix $K$, the $14\times8$ kernel embedding $M$, and its ordered free-row selector $J$. Require $B$ to have shape $14\times0$. The snapshot supplies candidate matrices and labels; physical regularity requires separate evidence. [Intermediate snapshot](../oblate-stationary-field-dynamics-timeout-2026-11-28/outputs/terminal-state-before-preservation.json) · [Original stationary export](../oblate-stationary-system-recheck-2026-11-27/outputs/stationary-operators.json)

Evaluate formal coefficients in $\mathbb Q(r,T)$ on $T(r)=\arctan(1/(10r))$. Their total derivative is $D=\partial_r-10(100r^2+1)^{-1}\partial_T$. Choose a positive radius $r_0$. Certify nonvanishing there of every denominator used in the proof, including those in the operators and discarded response. Exhibit six ordered columns $P$ of the **original** $C$ and certify $\det C_{:,P}(r_0,T(r_0))\ne0$. One exact method at rational $r_0$ is to enclose the physical arctangent by rational bounds and separate the relevant polynomial values from zero. Analyticity and continuity give a common open interval about $r_0$ on which all coefficients are regular and $C$ has rank six.

## Invariance and existence

Verify the exact identities $CM=0$, $JM=I_8$ and $SM=0$, where $S=DC+CK$. On the regular interval, $JM=I_8$ makes the eight columns of $M$ independent, and rank-nullity gives $\operatorname{im}M=\ker C$. Hence every row of $S$ annihilates $\ker C$ and lies in the row space of $C$. The certified minor gives an analytic multiplier
$$
L=S_{:,P}(C_{:,P})^{-1},\qquad S=LC .
$$
The equality holds on all fourteen columns by row-space membership and uniqueness of the coefficients on $P$.

For each real $y_0\in\mathbb R^8$, solve the analytic ODE $z'=Kz$ with $z(r_0)=M(r_0)y_0$. Along its unique local solution, $h=Cz$ satisfies $h'=Sz=Lh$ and $h(r_0)=0$. Uniqueness for this homogeneous equation gives $Cz=0$ throughout the germ. This argument requires neither $DM$ nor construction of a reduced dynamics matrix.

## Lift to the original equations and test the response

Let $P_v$ and $P_w$ select the upper and lower seven coordinates. Verify $(P_vK-P_w)M=0$ so that $v'=w$ throughout the germ and $v''=P_wKz$. Restore $X_0=0$ in the original profile order. With $E_j$ the three gauge-fixed coefficient blocks of the original operator, require the full matrix identity
$$
\bigl(E_0P_v+E_1P_w+E_2P_wK\bigr)M=0
$$
on **all 22 retained Euler rows**, including the $X_0$ equation. Each operator coefficient multiplies the indicated profile derivative.

Apply the same lift independently to the discarded row $D_{\rm disc}=6\mathscr E_{X,4}-\mathscr E_{C_{\rm polar},4}$, obtaining a row $d$ on the fourteen coordinates. A responding germ requires $(dM)_j(r_0)\ne0$ for a specified column $j$, with numerator and denominators certified at the same regular radius. Taking $y_0=e_j$ then gives a real solution germ with $D_{\rm disc}(r_0)\ne0$.

A certified zero row $dM$ would establish silence for this forward ODE family. A universal statement about solutions of the original stationary equations additionally requires a converse showing that every solution in the claimed sector lifts to this regular constrained ODE, including the treatment of the radial gauge. The argument concerns an open radial interval and provides no global boundary realization or nonlinear Dirac classification. [Certificate contract](certificate-contract.json)
