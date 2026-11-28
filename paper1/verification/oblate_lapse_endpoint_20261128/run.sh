#!/bin/sh
set -eu
cd "$(dirname "$0")/../.."
exec python3 verification/oblate_lapse_endpoint_20261128/run_bounded.py
