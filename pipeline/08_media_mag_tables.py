"""MAG tables for the media experiment (V, NV, NVNN co-assemblies) for data/.

Run once, on a machine that holds the project archive and anvi'o (anvio-dev):

    python pipeline/08_media_mag_tables.py /path/to/CyanoSphere_NAS_20260929 /path/to/anvio/bin

Each medium was co-assembled from all of its passages with MEGAHIT in the anvi'o metagenomics
workflow (metaSPAdes, also configured, failed on all three), and binned separately
(MEDIA/<VIT|NOVIT|NONANO>).
For the bins of at least medium quality (completion >= 50%, redundancy <= 10%) this writes

    data/media_mag_quality.csv      bin, medium, completion, redundancy, size, GTDB taxonomy
    data/media_mag_coverage.csv     bin x sample mean coverage, every passage profiled (P0 parent,
                                    P4, P5 parent, P10), sample names as in media_samples.csv style
    data/media_mag_modules.csv      module x bin pathwise and stepwise completeness

Module completeness is computed the same way for all three media, with anvi-estimate-metabolism
on an enzymes-txt per bin: the bin's gene calls (anvi-summarize bin_by_bin) joined to the KOfam
hits in the medium's contigs database, opened read-only. Nothing in the archive is written.
"""

import re
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

NAS, ANVIO = Path(sys.argv[1]), Path(sys.argv[2])
# the media co-assemblies (moved into 04_media_experiment/ when the archive was reorganised)
MEDIA = next(p for p in (NAS / "04_media_experiment/coassemblies_V_NV_NVNN", NAS / "MEDIA") if p.exists())
MASTER = next(p for p in (NAS / "04_media_experiment/media_bins_master.xlsx", NAS / "MEDIA-master.xlsx") if p.exists())
# per-bin results are cached here, so an interrupted run resumes where it stopped
CACHE = Path(__file__).resolve().parent / ".cache_08"
DATA = Path(__file__).resolve().parents[1] / "data"
MEDIUM = {"VIT": "V", "NOVIT": "NV", "NONANO": "NVNN"}
MIN_COMPLETION, MAX_REDUNDANCY = 50, 10
RANKS = ["t_domain", "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species"]


def sample_name(medium: str, col: str) -> str:
    # 'P10_Rep_R1N' -> 'V_P10_R1'; 'P0_Ori_R4' -> 'V_P0_parent'; 'P5_Ori_R2' -> 'V_P5_parent'
    p, kind, r = re.match(r"P(\d+)_(Rep|Ori)_R(\d)", col).groups()
    return f"{medium}_P{p}_" + ("parent" if kind == "Ori" else f"R{r}")


def taxonomy(cond: str) -> pd.DataFrame:
    if cond == "VIT":  # per-bin SCG taxonomy kept in the MEDIA master sheet, under MetaBAT names
        t = pd.read_excel(MASTER, sheet_name="VIT").set_index("bins")[RANKS]
        ren = pd.read_csv(MEDIA / "VIT/anvi_rename.txt", sep="\t").set_index("old_bin_name").new_bin_name
        t.index = t.index.map(ren)
        return t[t.index.notna()]
    f = {"NOVIT": "NOVIT/NV-taxonomy.txt", "NONANO": "NONANO/NVNN-taxonomy.txt"}[cond]
    t = pd.read_csv(MEDIA / f, sep="\t")
    return t.set_index(t.columns[0])[RANKS].replace("None", pd.NA)


def module_table(cond: str, bins: list[str]) -> pd.DataFrame:
    db = next((MEDIA / cond / "03_CONTIGS").glob("*.db"))
    with sqlite3.connect(f"file:{db}?mode=ro", uri=True) as con:
        ko = pd.read_sql("select gene_callers_id, accession from gene_functions where source='KOfam'", con)
    ko = ko.groupby("gene_callers_id").accession.first().to_dict()  # plain dict: safe across threads
    from concurrent.futures import ThreadPoolExecutor

    def one(b, tmp):
        done = CACHE / f"{b}.csv"
        if done.exists():
            return pd.read_csv(done)
        genes = pd.read_csv(MEDIA / cond / "SUMMARY/bin_by_bin" / b / f"{b}-gene_calls.txt",
                            sep="\t", usecols=["gene_callers_id"]).gene_callers_id
        hits = [(g, ko[g]) for g in genes if g in ko]
        enz = pd.DataFrame(hits, columns=["gene_id", "enzyme_accession"]).assign(source="KOfam")
        ef, prefix = Path(tmp) / f"{b}.txt", Path(tmp) / b
        enz.to_csv(ef, sep="\t", index=False)
        subprocess.run([str(ANVIO / "anvi-estimate-metabolism"), "--enzymes-txt", str(ef),
                        "-O", str(prefix), "--output-modes", "modules", "--include-zeros", "--quiet"],
                       check=True, capture_output=True)
        m = pd.read_csv(f"{prefix}_modules.txt", sep="\t")
        m["bin"] = b
        print(f"  {b}: {len(enz)} KOfam hits", flush=True)
        out = m[["module", "bin", "pathwise_module_completeness", "stepwise_module_completeness"]]
        out.to_csv(done, index=False)
        return out

    with tempfile.TemporaryDirectory() as tmp, ThreadPoolExecutor(max_workers=8) as pool:
        rows = list(pool.map(lambda b: one(b, tmp), bins))
    return pd.concat(rows)


def main() -> None:
    CACHE.mkdir(exist_ok=True)
    qual, cov, mods = [], [], []
    for cond, medium in MEDIUM.items():
        s = pd.read_csv(MEDIA / cond / "SUMMARY/bins_summary.txt", sep="\t").set_index("bins")
        keep = s[(s.percent_completion >= MIN_COMPLETION) & (s.percent_redundancy <= MAX_REDUNDANCY)]
        q = keep[["total_length", "num_contigs", "percent_completion", "percent_redundancy"]].copy()
        q.insert(0, "medium", medium)
        q = q.join(taxonomy(cond))
        qual.append(q)
        c = pd.read_csv(MEDIA / cond / "SUMMARY/bins_across_samples/mean_coverage.txt", sep="\t").set_index("bins")
        c = c.loc[keep.index]
        c.columns = [sample_name(medium, x) for x in c.columns]
        cov.append(c)
        print(f"{cond}: {len(keep)} of {len(s)} bins pass")
        mods.append(module_table(cond, list(keep.index)))
    pd.concat(qual).rename_axis("bin").to_csv(DATA / "media_mag_quality.csv")
    # one row per bin; a bin has coverage only in its own medium's samples
    pd.concat(cov).rename_axis("bin").to_csv(DATA / "media_mag_coverage.csv", float_format="%.5g")
    pd.concat(mods).to_csv(DATA / "media_mag_modules.csv", index=False, float_format="%.4g")


if __name__ == "__main__":
    main()
