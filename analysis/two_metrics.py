"""
The two contig-based metrics and the two-metric classification, for any set of metagenomes.

Input is what anvi'o already produces, one table per sample:

    anvi-estimate-metabolism -c SAMPLE-contigs.db -p SAMPLE/PROFILE.db -O SAMPLE \
                             --add-coverage --add-copy-number

(SAMPLE_modules.txt, long format), plus the number of populations in each sample (here,
97% RpL19 EcoPhylo clusters with non-zero coverage; the number of genomes estimated from
single-copy core genes works too) and the group of each sample.

    python two_metrics.py --modules kegg/*_modules.txt \
                          --populations populations.tsv \
                          --groups groups.tsv --test ENRICHMENT --reference SOURCE \
                          --out two_metrics_results

populations.tsv: sample<TAB>n_populations      groups.tsv: sample<TAB>group

Writes to --out:
    coverage_metric.tsv, count_metric.tsv     module x sample matrices
    differential_coverage.tsv, differential_count.tsv
    two_metric_classes.tsv                    one row per module, with its class

A table in the format of data/kegg_modules_*.csv (columns module, genome_name, avg_coverage,
stepwise_copy_number, stepwise_module_completeness) can be given to --modules instead.
Thresholds, test and fold change are those of the paper (common.py).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import common as C


def read_modules(paths: list[str]) -> pd.DataFrame:
    """Long table of module x sample with avg_coverage and stepwise_copy_number."""
    rows = []
    for p in paths:
        sep = "," if p.endswith((".csv", ".csv.gz")) else "\t"
        d = pd.read_csv(p, sep=sep)
        if "avg_coverage" not in d.columns:      # anvi'o names it <sample>_avg_coverage
            col = [c for c in d.columns if c.endswith("_avg_coverage")]
            if len(col) != 1:
                raise SystemExit(f"{p}: expected one *_avg_coverage column (run with --add-coverage)")
            d = d.rename(columns={col[0]: "avg_coverage"})
        if "stepwise_copy_number" not in d.columns:
            raise SystemExit(f"{p}: no stepwise_copy_number column (run with --add-copy-number)")
        rows.append(d[["module", "genome_name", "stepwise_module_completeness",
                       "avg_coverage", "stepwise_copy_number"]])
    return pd.concat(rows)


def metrics(mods: pd.DataFrame, npop: pd.Series) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Coverage-based: module coverage. Count-based: stepwise copy number / populations.
    Both only for modules at stepwise completeness >= C.PRESENCE in that sample."""
    d = mods[mods.stepwise_module_completeness >= C.PRESENCE]
    cov = d.pivot_table(index="module", columns="genome_name", values="avg_coverage", aggfunc="max")
    cnt = d.pivot_table(index="module", columns="genome_name", values="stepwise_copy_number", aggfunc="max")
    missing = set(cnt.columns) - set(npop.index)
    if missing:
        raise SystemExit(f"no population count for: {sorted(missing)}")
    cnt = cnt.div(npop.reindex(cnt.columns).astype(float), axis=1)
    return cov.fillna(0.0), cnt.fillna(0.0)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--modules", nargs="+", required=True)
    ap.add_argument("--populations", required=True)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--test", required=True, help="group whose increases are classified")
    ap.add_argument("--reference", required=True, help="group compared against")
    ap.add_argument("--out", default="two_metrics_results")
    a = ap.parse_args()

    npop = pd.read_csv(a.populations, sep="\t", header=None, index_col=0).iloc[:, 0]
    grp = pd.read_csv(a.groups, sep="\t", header=None, index_col=0).iloc[:, 0]
    g1, g2 = list(grp.index[grp == a.test]), list(grp.index[grp == a.reference])
    cov, cnt = metrics(read_modules(a.modules), npop)
    dcov, dcnt = C.differential(cov, g1, g2), C.differential(cnt, g1, g2)
    cls = C.niche_class(dcov, dcnt)

    out = Path(a.out)
    out.mkdir(exist_ok=True)
    cov.to_csv(out / "coverage_metric.tsv", sep="\t")
    cnt.to_csv(out / "count_metric.tsv", sep="\t")
    dcov.to_csv(out / "differential_coverage.tsv", sep="\t")
    dcnt.to_csv(out / "differential_count.tsv", sep="\t")
    t = pd.concat([dcov.add_prefix("coverage_"), dcnt.add_prefix("count_"), cls], axis=1)
    t.to_csv(out / "two_metric_classes.tsv", sep="\t")
    print(cls.value_counts().to_string())


if __name__ == "__main__":
    main()
