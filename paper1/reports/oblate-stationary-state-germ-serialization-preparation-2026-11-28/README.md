# Exact endpoint serialization repair

This preparation addresses the failure in run `9afcd73b-3499-4f08-8185-0f162c716cfb`, preserved in the [failure archive](../oblate-stationary-state-germ-serialization-failure-2026-11-28/README.md). The run selected a responding candidate at radius 1, arctangent term count 25, and zero-based response column 1. It then failed in `serialize_minor` because Python refused decimal conversion of an interval endpoint exceeding its 4,300-digit limit. No physical germ or response certificate was exported.

## Change and independent review

The new stdlib `exact_decimal.py` formats an already-computed integer by repeated division into base-10^9 chunks. Each chunk has at most nine digits; their ordered concatenation gives the exact decimal string. A rational endpoint uses its exact integer numerator and positive denominator. This does not alter Python's global digit-limit setting. All determinant, pivot, denominator, arctangent-control and response endpoint serialization paths use that formatter. Polynomial coefficient and provenance strings retain their previous handling.

The supervisor snapshots and authenticates the helper and its test file as well as the previous eight source/config/view files. Its process lifecycle, input pins and terminal-evidence protocol are unchanged. Independent static review `msg_5b7209fd` found no blocker and authenticated the five changed/new Python files. No original operator, state matrix, field arithmetic, interval operation, pivot choice, candidate order or success predicate changed.

## Validation and remaining work

[validation.json](validation.json) captures six passing stdlib formatter tests (0.008 seconds), 27 passing mocked supervisor tests (0.354 seconds), shell syntax success, and source hashes. The formatter tests check signs, chunk boundaries, repeated decimal patterns, values above 4,300 digits, rational JSON payloads, unchanged interpreter limits, and invalid input rejection. No scientific module was imported. An AST comparison found all 91 existing checker functions equivalent after normalizing only old/new endpoint formatting expressions and removing the single added serialization-control call. Eighty-four existing functions are byte-equivalent at the AST level without normalization.

The next tracked child will run 20 mathematical controls, including the newly prepared 5,001-digit rational determinant/pivot serialization control. That new scientific control remains unexecuted. The same two scientific input pins, 900-second wall cap, sampled 2 GiB process-group RSS threshold, 256-term bound, 17-stage view and 2026-11-28 16:29 UTC cutoff remain in force. The preceding invocation's 261-second duration is evidence for the unchanged calculation up to the export failure; successful final serialization, acceptance and runtime remain unmeasured.

An accepted physical result still requires the complete certificate, matching input and source hashes, clean process-group shutdown, and the terminal marker authenticating both execution and metrics. Preparation and the previous selected-candidate log do not themselves certify a physical minor, regular germ, nonzero response, global boundary domain, nonlinear consistency or a degree-of-freedom count. The separate transcendental-method preparation is excluded from this invocation.
