#!/bin/sh
set -eu
cd "$(dirname "$0")/../.."
exec python3 verification/oblate_stationary_germ_20261128/run_bounded.py
