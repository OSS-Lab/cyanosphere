"""Export the raw tables behind the EcoPhylo heatmaps and SI figures to data/.

Run once, on a machine that holds the project archive (CyanoSphere_NAS_20260929/):

    python pipeline/07_export_raw_tables.py /path/to/CyanoSphere_NAS_20260929

It reads anvi'o PROFILE.db files directly with sqlite3 (no anvi'o needed) and
writes plain CSVs, so every heatmap that was exported from anvi-interactive has
its underlying numbers in the repository.

Writes
    data/ecophylo_coverage.csv   run, cluster, sample, mean_coverage_Q2Q3, detection
    data/ecophylo_clusters.csv   run, cluster, cluster_size, GTDB taxonomy of the representative
    data/assembly_stats.csv      per-sample metaSPAdes assembly statistics (anvi-display-contigs-stats)
    data/media_samples.csv       media-experiment sample IDs with condition and passage
    data/raw_kegg/modules_unbinned_{enrichments,media}.csv.gz
                                 the full anvi-estimate-metabolism output on the unbinned contigs of
                                 every sample (enzyme hits, gene caller IDs, per-gene coverage and
                                 detection, copy numbers); kegg_modules_*.csv are column subsets of these
"""

import re
import sqlite3
import sys
from pathlib import Path

import pandas as pd

NAS = Path(sys.argv[1])
DATA = Path(__file__).resolve().parents[1] / "data"

# EcoPhylo runs: run name -> merged profile
RUNS = {
    "rps16_source_enrichment": "03_ecophylo/RpS16_source_vs_enrichment/06_MERGED/Bacteria_71_Ribosomal_S16/PROFILE.db",
    "rps16_media_endpoints": "03_ecophylo/RpS16_media_P10/06_MERGED/Bacteria_71_Ribosomal_S16/PROFILE.db",
    "rps16_media_allpassages": "03_ecophylo/RpS16_media_all_passages/06_MERGED/Bacteria_71_Ribosomal_S16/PROFILE.db",
    "rpl19_media_allpassages": "03_ecophylo/RpL19_media_all_passages/METAGENOMICS_WORKFLOW/06_MERGED/Bacteria_71_Ribosomal_L19/PROFILE.db",
    # earlier RpL19 run that also profiled the P0 and P3 passages of the enrichment lineages
    # (+ M_70 and M_76, mapped to the same references with the run's settings on 7 Oct 2026,
    # see old_runs/ECOPHYLO2_RpL19_plus_M70_M76/rerun_add_M70_M76.sh)
    "rpl19_enrichment_lineages": "03_ecophylo/old_runs/ECOPHYLO2_RpL19_plus_M70_M76/PROFILE.db",
}
# runs whose merged profile carries no EcoPhylo taxonomy or tree: take them from the original run
ANNOTATION_FROM = {"rpl19_enrichment_lineages": "03_ecophylo/old_runs/ECOPHYLO2_RpL19/METAGENOMICS_WORKFLOW/06_MERGED/Ribosomal_L19/PROFILE.db"}
# sample names that were mistyped in a run's samples.txt (checked against the read paths)
RENAME = {"M_869": "M_86", "WAMOS_01": "W_01"}
CLUSTER_KEYS = ["cluster_size", "genomic_seq_in_cluster", "percent_identity",
                "t_domain", "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]


def cluster_name(item: str) -> str:
    # 'V_P10_R1_Bacteria_71_Ribosomal_S16_000000000002_split_00001' -> drop the split suffix
    return re.sub(r"_split_\d+$", "", item)


def export_ecophylo() -> None:
    cov, clu = [], []
    for run, rel in RUNS.items():
        with sqlite3.connect(f"file:{NAS / rel}?mode=ro", uri=True) as db:
            q2q3 = pd.read_sql("select item, layer as sample, value as mean_coverage_Q2Q3 from mean_coverage_Q2Q3_splits", db)
            det = pd.read_sql("select item, layer as sample, value as detection from detection_splits", db)
        with sqlite3.connect(f"file:{NAS / ANNOTATION_FROM.get(run, rel)}?mode=ro", uri=True) as adb:
            extra = pd.read_sql("select item_name as item, data_key, data_value from item_additional_data", adb)
            tree = pd.read_sql("select data from item_orders where name like '%_tree:NA:NA'", adb).data
        c = q2q3.merge(det, on=["item", "sample"], how="outer")
        c.insert(0, "run", run)
        c = c.rename(columns={"item": "cluster"})
        c["sample"] = c["sample"].replace(RENAME)
        c["cluster"] = c["cluster"].map(cluster_name)
        cov.append(c)
        t = extra[extra.data_key.isin(CLUSTER_KEYS)].pivot(index="item", columns="data_key", values="data_value")
        t.index = t.index.map(cluster_name)
        t = t.reindex(columns=CLUSTER_KEYS).rename_axis("cluster").reset_index()
        t.insert(0, "run", run)
        clu.append(t)
        if len(tree):  # the EcoPhylo marker phylogeny, used to place clusters without SCG taxonomy
            (DATA / "ecophylo_trees").mkdir(exist_ok=True)
            (DATA / "ecophylo_trees" / f"{run}.nwk").write_text(re.sub(r"_split_\d+", "", tree.iloc[0]))
        print(f"{run}: {c['cluster'].nunique()} clusters x {c['sample'].nunique()} samples")
    pd.concat(cov).to_csv(DATA / "ecophylo_coverage.csv", index=False, float_format="%.6g")
    pd.concat(clu).to_csv(DATA / "ecophylo_clusters.csv", index=False)


def export_assembly_stats() -> None:
    rows = {}
    for f in sorted((NAS / "01_coassemblies_STRONG/per_sample_assembly_stats").glob("contigs_stats_*.txt")):
        s = pd.read_csv(f, sep="\t", header=None, index_col=0).iloc[:, 0]
        rows[f.stem.replace("contigs_stats_", "")] = s.drop("contigs_db", errors="ignore")
    out = pd.DataFrame(rows).T.rename_axis("sample")
    out.to_csv(DATA / "assembly_stats.csv")
    print(f"assembly_stats: {len(out)} samples x {out.shape[1]} statistics")


def export_media_samples() -> None:
    """Media-experiment sample key from the sequencing log (README_allsamples.txt)."""
    pat = re.compile(r"^(M_\d+)\s+(WAMOS nutrients?\s.*)$")
    rows = []
    for line in (NAS / "08_sample_metadata_ncbi/sequencing_log_all_samples.txt").read_text().splitlines():
        m = pat.match(line.strip())
        if not m:
            continue
        sid, desc = m.group(1), re.sub(r"\s+", " ", m.group(2)).strip()
        d = desc.lower()
        cond = ("NVNN" if ("no nitrate" in d or re.search(r"\bnvn", d)) else
                "NV" if ("no vit" in d or re.search(r"\bnv", d)) else "V")
        p = re.search(r"\bP(\d+)\b", desc)
        passage = f"P{p.group(1)}" if p else ("parent" if "parent" in d else "intermediate (25/1/2024)")
        rep = re.search(r"(?:Repeat|Rep|R)\s?(\d)\b", desc)
        rows.append((sid, cond, passage, rep.group(1) if rep else "", desc))
    out = pd.DataFrame(rows, columns=["sample", "condition", "passage", "replicate", "lab_description"])
    # P10 end-points carry their SRA / manuscript names (V_P10_R1 ...)
    p10 = out.passage.eq("P10")
    out.insert(1, "sra_name", "")
    out.loc[p10, "sra_name"] = out.loc[p10, "condition"] + "_P10_R" + out.loc[p10, "replicate"]
    out.to_csv(DATA / "media_samples.csv", index=False)
    print(f"media_samples: {len(out)} samples")


ANVIO_COLS = ["module", "sample", "module_name", "module_class", "module_category",
              "module_subcategory", "module_definition", "stepwise_module_completeness",
              "stepwise_module_is_complete", "pathwise_module_completeness",
              "pathwise_module_is_complete", "proportion_unique_enzymes_present",
              "enzymes_unique_to_module", "unique_enzymes_hit_counts", "enzyme_hits_in_module",
              "gene_caller_ids_in_module", "gene_coverages", "avg_coverage", "gene_detection",
              "avg_detection", "pathwise_copy_number", "stepwise_copy_number"]


def export_kegg_raw() -> None:
    """anvi-estimate-metabolism --add-coverage --add-copy-number, one row per module per sample."""
    out = DATA / "raw_kegg"
    out.mkdir(exist_ok=True)
    src = NAS / "05_metabolism_KEGG/unbinned_contigs_modules"
    e = pd.read_csv(src / "all_cyanosphere_contigs_KEGG_modules.csv", low_memory=False,
                    encoding="utf-8-sig").rename(columns={"genome_name": "sample"})
    e[ANVIO_COLS].to_csv(out / "modules_unbinned_enrichments.csv.gz", index=False)
    parts = []
    for f in sorted((NAS / "05_metabolism_KEGG/unbinned_contigs_taxonomy_and_modules/INPUTS_MEDIA").glob("*_modules.txt")):
        d = pd.read_csv(f, sep="\t", low_memory=False)
        smp = f.name.replace("_modules.txt", "")
        d = d.rename(columns={c: c[len(smp) + 1:] for c in d.columns if c.startswith(smp + "_")})
        d["sample"] = smp
        parts.append(d)
    m = pd.concat(parts)
    m[ANVIO_COLS].to_csv(out / "modules_unbinned_media.csv.gz", index=False)
    print(f"raw_kegg: {len(e)} enrichment rows, {len(m)} media rows")


if __name__ == "__main__":
    export_ecophylo()
    export_assembly_stats()
    export_media_samples()
    export_kegg_raw()
