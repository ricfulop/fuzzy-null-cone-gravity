#!/bin/sh
set -eu
cd "$(dirname "$0")/../.."
exec python3 verification/oblate_terminal_embedding_probe_20261128/run_bounded.py
