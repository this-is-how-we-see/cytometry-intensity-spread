#!/bin/sh
# One command for the whole pilot. Run from anywhere.
set -e
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(dirname "$HERE")
cd "$ROOT"
[ -d code/.venv ] || python3 -m venv code/.venv
. code/.venv/bin/activate
pip install -q -e code 2>/dev/null || pip install -q numpy scipy pandas scikit-learn matplotlib flowio flowsom scanpy leidenalg igraph rdata
python3 -I code/mirror_zzph.py data/FR-FCM-ZZPH || echo "FlowRepository mirror incomplete (login needed); see data/README.md"
for n in Nilsson_rare Mosmann_rare; do
  [ -f data/HDCytoData/$n.npz ] || python3 -I code/rda_extract.py data/HDCytoData/${n}_SE.rda data/HDCytoData/$n.npz
done
python3 -I code/pilot.py . 2>&1 | grep -v DEBUG | tee results/pilot.log
python3 -I code/figures.py .
