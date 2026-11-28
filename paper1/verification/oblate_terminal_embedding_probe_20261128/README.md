# Diagnostic terminal constraint kernel

This prepared probe measures the cost of constructing a right kernel of the six stored constraints on fourteen state coordinates. The source is the step-02 checkpoint of the fixed-$\epsilon=1/100$ stationary reduction. The interrupted invocation reached terminal embedding after writing this checkpoint. Its aggregate stage timer does not identify which operation occupied that time.

## Input and output

The entry point accepts three required arguments. The future supervisor supplies the independently recorded lowercase SHA-256, the original checkpoint or a byte-identical snapshot, and a fresh output filename. Paths below are shell variables supplied by that supervisor.

```sh
python verification/oblate_terminal_embedding_probe_20261128/probe_terminal_embedding.py \
  --checkpoint "$CHECKPOINT_PATH" \
  --checkpoint-sha256 "$PINNED_CHECKPOINT_SHA256" \
  --output "$CANDIDATE_OUTPUT"
```

The archived source is [dae-step-02.json](../../reports/oblate-stationary-closure-terminal-timeout-2026-11-28/outputs/dae-step-02.json). Its SHA-256 read during preparation was `ca630b313848bc4230f72f05642d0821297e50322ed2f0dc5ccbaa43eb104ccc`. The script checks the supplied pin before parsing those same bytes. It verifies the checkpoint type, step index, declared row count and matrix shape. Only the original `constraints` rows enter the algebra. The script uses SymPy 1.14's `DomainMatrix` API and records the actual Python and SymPy versions.

## Exact calculation

1. Parse the stored integer-coefficient rational expressions, then convert them to a sparse matrix $C$ over $F=\mathbb Q(r,T)$ with $T=\arctan(1/(10r))$. All subsequent matrix arithmetic takes place over the field or its associated polynomial ring.
2. Call `C.rref_den(method="CD", keep_domain=False)`. The returned numerator $U$ and common denominator $d$ belong to the same domain. Convert both to $F$ and form $R=U/d$. Check the unit pivot columns and trailing zero rows.
3. For pivot columns $p_i$ and ordered free columns $f_j$, set $M_{f_i,j}=\delta_{ij}$ and $M_{p_i,j}=-R_{i,f_j}$. Check the free-row identity. Multiply the preserved **original** $C$ by $M$ in the field and require an exact zero matrix. Together with the pivot count, these checks certify a basis of the formal right kernel of this checkpoint matrix.
4. Serialize the candidate $M$, its shape, zero-based pivot and free-column indices, input digest, and checks. The output is published atomically at a fresh filename. The phase log gives begin/end events for input verification, manufactured controls, expression parsing, field conversion, RREF, denominator normalization, embedding construction, the original-matrix product, and serialization. Serialization includes publication; its end timing is in the log. Earlier phase durations are also saved in the output.

The prepared manufactured controls cover a known rational kernel with a dependent row, a perturbed candidate with a nonzero original-matrix product, a zero matrix, empty matrix dimensions, full column rank, and rejection of floating-point literals even when simplification would erase them.

## Interpretation

The RREF derives coordinates for the kernel of the already stored $C$. The original $C$ remains available for the independent product check. The historical control pivots and the rows differentiated during discovery are preserved because the probe begins after the last recorded elimination step. Reducing the state block during control elimination could change the rows that subsequently enter differentiation.

A successful diagnostic establishes $C M=0$ and the rank and nullity of that formal matrix. Certification of a stationary solution family still requires constraint preservation, reconstruction of the reduced radial dynamics, direct substitution in all 22 original retained Euler equations (including the $X_0$ equation), the independent discarded-response check, and analytic regularity at a physical radius. The candidate is a diagnostic artifact. The original reduction remains interrupted, and runtime for this method is unmeasured at preparation.
