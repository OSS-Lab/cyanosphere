"""
Link each dereplicated MAG to its source bin by sequence, and count detected MAGs per sample.

The 88 MAGs are named taxonomically (Aestuariivirga_sp020252405), while the coverage
tables name them by bin (OSS1_Bin_119); matching by name resolves only 63 of 88. Here
each MAG's contig headers are compared with every binned FASTA, and a bin is the source
of a MAG when more than 90% of the contigs match in both directions. This is
unambiguous for all 88 (one bin each, no ties).

The per-sample MAG count (`N_MAGs`, a denominator in the robustness table) is then the
number of the 88 MAGs detected in each metagenome at a stated coverage threshold. The
dereplicated coverage table holds all eight enrichment end-points but only four of the
eight reservoir samples (W_06, W_10, W_13, W_15), so the count is complete for the
enrichments only.

Needs the project archive (the binned FASTAs under */FASTA/ and the dereplicated coverage
table), and writes data/mag_bin_map.tsv and data/mag_sample_coverage.tsv, so that
everything downstream needs only this repository.

    python pipeline/05_link_mags_to_bins.py /path/to/project_archive
"""
from __future__ import annotations

import collections
import glob
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import common as C

# A MAG counts as detected in a metagenome above this mean coverage. Reported
# across a sweep rather than asserted, because it is a free parameter.
DETECTION = 1.0
SWEEP = [0.1, 0.25, 0.5, 1.0, 2.0, 5.0]
OVERLAP = 0.9  # fraction of contigs that must match, in both directions


def contig_sets_from_mag_fasta(path: Path) -> dict[str, set[str]]:
    """MAG name -> contig headers, from the concatenated `NAME__NODE_...` export."""
    mags: dict[str, set[str]] = collections.defaultdict(set)
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name, node = line[1:].strip().split("__", 1)
                mags[name].add(node.split()[0])
    return mags


def contig_set(path: str) -> set[str]:
    with open(path) as fh:
        return {l[1:].strip().split()[0] for l in fh if l.startswith(">")}


def build_map(root: Path, mag_fasta: Path) -> pd.DataFrame:
    mags = contig_sets_from_mag_fasta(mag_fasta)
    bins = {}
    for p in glob.glob(str(root / "*" / "FASTA" / "*.fasta")):
        bins[os.path.basename(p)[:-6]] = contig_set(p)
    print(f"{len(mags)} MAGs against {len(bins)} binned FASTAs")

    rows, ambiguous, unmatched = [], [], []
    for name, nodes in mags.items():
        hits = [
            (b, len(nodes & bn))
            for b, bn in bins.items()
            if bn and len(nodes & bn) / len(bn) > OVERLAP
            and len(nodes & bn) / len(nodes) > OVERLAP
        ]
        if len(hits) == 1:
            rows.append({"MAG_name": name, "bin": hits[0][0], "n_contigs": len(nodes)})
        elif hits:
            ambiguous.append(name)
        else:
            unmatched.append(name)
    if ambiguous or unmatched:
        raise SystemExit(f"mapping not clean: {len(ambiguous)} ambiguous, "
                         f"{len(unmatched)} unmatched")
    print(f"all {len(rows)} MAGs mapped to exactly one bin")
    return pd.DataFrame(rows)


def main() -> None:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else C.REPO.parent)
    mag_fasta = Path(sys.argv[2]) if len(sys.argv) > 2 else None
    if mag_fasta is None or not mag_fasta.exists():
        raise SystemExit("pass the concatenated 88-MAG contigs FASTA as argument 2 "
                         "(anvi-export-contigs on each MAG db, headers NAME__NODE_...)")

    m = build_map(root, mag_fasta)
    m.to_csv(C.DATA / "mag_bin_map.tsv", sep="\t", index=False)

    drep = pd.read_csv(root / "drep_cyanosphere_samples_coverages_taxonomies_quality_merged.csv")
    samples = [c for c in drep.columns if c.startswith(("M_", "W_"))]
    cov = (drep.set_index("mag").loc[m.bin, samples]
           .set_axis(m.MAG_name.values).sort_index())
    cov.index.name = "MAG_name"
    cov.to_csv(C.DATA / "mag_sample_coverage.tsv", sep="\t")
    print(f"coverage table: {cov.shape[0]} MAGs x {cov.shape[1]} samples")
    have = [s for s in C.SOURCE + C.ENRICHMENT if s in cov.columns]
    print(f"of the paper's 16 metagenomes, {len(have)} are present: {have}")
    print(f"missing: {sorted(set(C.SOURCE + C.ENRICHMENT) - set(have))}")

    # --- the count, over a threshold sweep ------------------------------------
    sweep = pd.DataFrame(
        {f"cov>{t}": (cov[have] > t).sum(axis=0) for t in SWEEP}
    )
    sweep.index.name = "sample"

    old = (pd.read_csv(C.DATA / "kegg_modules_enrichments.csv")
           .groupby("genome_name")["N_MAGs"].first())
    sweep["published_N_MAGs"] = old.reindex(sweep.index)
    sweep.to_csv(C.OUT / "n_mags_rebuilt.csv")
    print("\n--- MAGs detected per metagenome, by coverage threshold ---")
    print(sweep.to_string())

    chosen = sweep[f"cov>{DETECTION}"]
    both = sweep.dropna(subset=["published_N_MAGs"])
    if len(both) > 1:
        r = np.corrcoef(both[f"cov>{DETECTION}"], both.published_N_MAGs)[0, 1]
        print(f"\nAt coverage > {DETECTION}: mean {chosen.mean():.1f} MAGs per sample "
              f"(range {chosen.min()} to {chosen.max()}).")
        print(f"Correlation with the published N_MAGs column: r = {r:.3f}")
        print("The published column is not reproduced by any threshold in the sweep;")
        print("use the rebuilt count and state the rule.")


if __name__ == "__main__":
    main()
