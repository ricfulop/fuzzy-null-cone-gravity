#!/bin/sh
set -eu
cd "$(dirname "$0")/../.."
exec python3 verification/oblate_null_shift_atlas_20261128/run_bounded.py
