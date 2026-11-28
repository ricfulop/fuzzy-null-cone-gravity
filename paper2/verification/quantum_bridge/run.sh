#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
printf '::orx:: {"phase":"setup","step":0,"of":1,"unit":"checks"}\n'
python3 --version
printf '::orx:: {"phase":"setup","step":1,"of":1,"unit":"checks"}\n'
exec python3 -u verification/quantum_bridge/run.py
