#!/bin/bash
#SBATCH --job-name=media_function
#SBATCH --output=%x-%j.out
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
# Functional partition (Fig 4E), its NVNN vs nitrate-replete test, functional
# trajectories and the Fig 4 composite, all on the shared genome reference.
# Runs from a copy of the repo's analysis/ and data/ next to the mapping outputs.
set -euo pipefail
cd /springbrook/share/lifesci/lfsnhf/metagenomics/runs/media_reads_to_MAGs/natmicro_analysis/analysis
export MPLBACKEND=Agg
PY=/springbrook/share/lifesci/lfsnhf/conda_envs/media_mapping/bin/python
for pkg in pymupdf:pymupdf skbio:scikit-bio; do $PY -c "import ${pkg%%:*}" 2>/dev/null || $PY -m pip install -q ${pkg#*:}; done
$PY -c "import media as m; m.run_shared_partition(); m.run_partition_test(); m.run_functional_trajectories()"
$PY -c "import figures as f; r = f.compose_fig4(); print(\"E_check\", r[\"partition\"]); f.compose_si_media_partition()"
echo DONE
