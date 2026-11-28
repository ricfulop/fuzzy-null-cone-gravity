# Exact principal block of the stationary Paper I system

The single local invocation passed at fixed oblateness epsilon=1/100. After
checking the radial gauge identity against all 22 retained rows and the discarded
response, the highest-radial-derivative block has exact formal rank 4. Its pivot
profiles are X2, C_polar2, a0, a2; lambda0, lambda2, V_polar2 remain free acceleration
columns at this initial stage. The eight selected equations are not eight
independent constraints. This is not a Dirac rank or a physical DOF count.

Both manufactured controls passed: a hidden differentiated constraint chain
retains a free response control, and a variable-coefficient constraint forces
its designated response to zero. The input conversion was reversible and
numerator/denominator polynomials were checked over Q(r,T), where
T=atan(1/(10r)) and T'=-10/(100r^2+1).

The source is 28973f993fec42a76cbea3fed3ec767b88f40c17; tracked experiment
bbeaa7e4-0ca4-4f18-8d91-85cc9e1331f7, run
e04ffc4a-d9c3-4c8a-ad61-3bf8da079ff3. The upstream completed export is independently
pinned by its SHA-256 in the source and launch record. Source snapshots match
both the recorded commit and all before/after digests.

Exit 0; 5.371265583 seconds supervised elapsed, 4.281256 seconds inside the
reduction recorder. The largest sampled process-group RSS was 109608960 bytes
(104.53 MiB); no instantaneous peak or hard memory cap was measured.
One scientific child; no GPU or remote compute. Temporary user authorization
msg_4c4bcc8d covers this invocation.

## Evidence and limits

- [Raw scientific log](verification.log), [full engine log](engine.log)
- [Execution and hashes](execution.json), [resources](resource-samples.json)
- [Exact principal matrices and pivots](outputs/principal-block.json)
- [Manufactured controls](outputs/manufactured-controls.json)
- [Invocation identity](outputs/invocation.json), [progress](outputs/progress.jsonl)

No DAE closure, actual radial solution germ, responding non-gauge seed, or
nonlinear compatibility value was computed. A future closure must keep the
original 22 rows, all differentiated constraints and remaining free controls.
An actual germ additionally requires regular coefficients on an open interval.
