"""
Supplementary tables built from the derived data (S1 accessions, S4 carriers,
S5 MAGs) and the symclatron lifestyle predictions for the 88 MAGs.

Run:  python supplementary.py      (from analysis/, after the modules before it in
      the README order)
"""
from __future__ import annotations

from pathlib import Path
import pandas as pd
import common as C
import sys
import numpy as np
from scipy.stats import fisher_exact, mannwhitneyu

# ============================================================================
# Supplementary tables: run_tables()
# ----------------------------------------------------------------------------
#   S1  accession numbers per sample (data/sra_accessions.csv, from NCBI
#       BioProject PRJNA1432497, fetched 2026-08-17); samples not yet deposited
#       are listed as such
#   S4  per replicate, who carries the variably encoded modules: which
#       cyanobacterial MAGs are present in each enrichment, which of them encode
#       the module, and which heterotrophs encode it instead
#   S5  MAG accessions and quality
#
# Reads data/sra_accessions.csv, data/media_samples.csv, data/mag_quality.tsv,
# data/mag_biosample_accessions.tsv, data/mag_biosample_reference.tsv and the
# MAG tables of common.mag_tables().
# ============================================================================

HERE = C.OUT


DATA = C.DATA


THRESHOLD = C.ENCODING


SAMPLES = C.ENRICHMENT


# modules whose cyano/heterotroph split varies across replicates (Fig 3E)
VARIABLE = {
    "M00123": "Biotin biosynthesis",
    "M00122": "Cobalamin biosynthesis (aerobic branch)",
    "M00924": "Cobalamin biosynthesis (anaerobic)",
    "M00899": "Thiamine salvage",
    "M00531": "Assimilatory nitrate reduction",
    "M00970": "Proline degradation",
    "M00133": "Polyamine biosynthesis (via arginine)",
}


def tex_escape(s: str) -> str:
    return str(s).replace("_", r"\_").replace("&", r"\&").replace("%", r"\%")


# ---------------------------------------------------------------- S1 accessions
def s1_accessions() -> None:
    run = pd.read_csv(DATA / "sra_accessions.csv")
    run = run.set_index("SampleName")
    # media lineages sequenced before P10 (parent cultures and the intermediate
    # passage), named as in media_mag_coverage.csv
    med = pd.read_csv(DATA / "media_samples.csv")
    early = []
    for cond in ("V", "NV", "NVNN"):
        m = med[(med.condition == cond) & ~med.passage.isin(["P10", "P12"])]
        for _, r in m.iterrows():
            if r.passage == "parent":
                early.append((f"{cond}_P0_parent", r["sample"]))
            elif pd.notna(r.replicate):
                early.append((f"{cond}_P4_R{int(r.replicate)}", r["sample"]))
            else:
                early.append((f"{cond}_P5_parent", r["sample"]))
    order = [(s, s) for s in SAMPLES] + \
            [(f"W_{i:02d}", f"W_{i:02d}") for i in (1, 6, 7, 10, 13, 14, 15, 18)] + \
            sorted(early, key=lambda x: (["V", "NV", "NVNN"].index(x[0].split("_")[0]), x[0])) + \
            [(f"V_P10_R{i}", f"V_P10_R{i}") for i in range(1, 6)] + \
            [(f"NV_P10_R{i}", f"NV_P10_R{i}") for i in range(1, 6)] + \
            [(f"NVNN_P10_R{i}", f"NVNN_P10_R{i}") for i in range(1, 5)]
    rows = []
    for s, lab in order:
        if s in run.index or lab in run.index:
            k = s if s in run.index else lab
            rows.append(f"{tex_escape(s)} & {run.loc[k, 'BioSample']} & {run.loc[k, 'Run']} \\\\")
        else:
            rows.append(f"{tex_escape(s)} & \\multicolumn{{2}}{{l}}{{\\textit{{in submission}}}} \\\\")
    tex = r"""% Sequence accessions, BioProject PRJNA1432497. Regenerate after every
% deposit: add the runs to data/sra_accessions.csv (Run,BioSample,SampleName).
\begin{longtable}{@{}lll@{}}
\caption{\textbf{Sequence accessions for all metagenomic samples.} BioProject
PRJNA1432497. M, long-term enrichments; W, reservoir samples; V, NV and NVNN,
media experiment: parent culture of each medium (P0), intermediate passage
(replicate lineages R1--R5, P4, and the continuing parent line, P5) and end
point (P10).}
\label{stab:accessions}\\
\toprule
Sample & BioSample & SRA run \\
\midrule
\endfirsthead
\multicolumn{3}{@{}l}{\emph{Supplementary Table~\ref{stab:accessions}, continued}}\\
\toprule
Sample & BioSample & SRA run \\
\midrule
\endhead
\bottomrule
\endfoot
""" + "\n".join(rows) + r"""
\end{longtable}
"""
    (HERE / "s1_accessions.tex").write_text(tex)
    missing = [s for s, lab in order if s not in run.index and lab not in run.index]
    print(f"Accessions: {len(order) - len(missing)}/{len(order)} samples deposited; "
          f"not yet: {', '.join(missing)}")


# ------------------------------------------------------------------------- S4
def s4_carriers() -> None:
    data, cov, tax = C.mag_tables()
    cyano = tax.phylum == C.CYANO_PHYLUM

    def lineage(g):
        """Lowest named GTDB rank; an unnamed genus or family is labelled from the rank above."""
        for rank in ("genus", "family", "order"):
            if pd.notna(tax.loc[g, rank]):
                return str(tax.loc[g, rank]) if rank == "genus" else f"unclassified {tax.loc[g, rank]}"
        return "unclassified"

    long_rows, tex_rows = [], []
    for m, label in VARIABLE.items():
        row = data.loc[m]
        enc = set(row[row >= THRESHOLD].index)
        for s in SAMPLES:
            present = cov.index[cov[s] >= C.MAG_DETECT]
            cy_present = [x for x in present if cyano[x]]
            cy_enc = [x for x in cy_present if x in enc]
            het_enc = [x for x in present if not cyano[x] and x in enc]
            cy_cov = float(cov.loc[cy_enc, s].sum())
            het_cov = float(cov.loc[het_enc, s].sum())
            tot = cy_cov + het_cov
            long_rows.append(dict(
                module=m, module_label=label, sample=s,
                cyano_mags_present=len(cy_present),
                cyano_lineages_present="; ".join(sorted({lineage(x) for x in cy_present})),
                cyano_encoders=len(cy_enc),
                cyano_encoder_lineages="; ".join(sorted({lineage(x) for x in cy_enc})),
                other_encoders=len(het_enc),
                other_encoder_families="; ".join(sorted({
                    str(tax.loc[x, "family"]) for x in het_enc
                    if pd.notna(tax.loc[x, "family"])})[:12]),
                cyano_pct=round(100 * cy_cov / tot, 1) if tot > 0 else None,
            ))
    s4 = pd.DataFrame(long_rows)
    s4.to_csv(HERE / "s4_per_replicate_carriers.csv", index=False)

    # compact LaTeX: one row per (module, replicate), lineages abbreviated
    for m, label in VARIABLE.items():
        sub = s4[s4.module == m]
        tex_rows.append(r"\multicolumn{6}{@{}l}{\emph{" + f"{m}: {label}" + r"}} \\")
        for _, r in sub.iterrows():
            het = r.other_encoder_families
            het_short = "; ".join([x for x in str(het).split("; ")][:3])
            if str(het).count(";") >= 3:
                het_short += r"; \ldots"
            cy = r.cyano_encoder_lineages or "--"
            pct = "--" if pd.isna(r.cyano_pct) else f"{r.cyano_pct:.0f}"
            tex_rows.append(
                f"{tex_escape(r['sample'])} & {r.cyano_mags_present} & "
                f"{tex_escape(cy) if cy else '--'} & {r.cyano_encoders} & "
                f"{r.other_encoders} & {pct} \\\\")
        tex_rows.append(r"\addlinespace")

    # longtable: 7 modules x 8 replicates do not fit on one page
    cols = r"Replicate & Cyano.\ MAGs & Cyano.\ encoder lineage(s) & Cyano.\ encoders & Other encoders & Cyano.\ \% \\"
    tex = r"""{\scriptsize
\begin{longtable}{@{}lr>{\raggedright\arraybackslash}p{5.2cm}rrr@{}}
\caption{\textbf{Per-replicate carriers of the variably encoded modules.} For each
enrichment replicate: the number of cyanobacterial MAGs detected (mean read coverage
$\geq 1$), which of those lineages encode the module (pathwise module completeness
$\geq 0.75$), how many cyanobacterial and other bacterial MAGs encode it, and the
resulting cyanobacterial share of the coverage of all detected MAGs encoding it (the
quantity plotted in Figure~3E; -- where no detected MAG encodes the module).
MAG coverage is linked to each genome by sequence. Full lineage lists, including the
other bacterial families, are in the accompanying CSV.}
\label{stab:carriers}\\
\toprule
""" + cols + r"""
\midrule
\endfirsthead
\multicolumn{6}{@{}l}{\emph{Supplementary Table~\ref{stab:carriers}, continued}}\\
\toprule
""" + cols + r"""
\midrule
\endhead
\bottomrule
\endfoot
""" + "\n".join(tex_rows) + r"""
\end{longtable}
}
"""
    (HERE / "s4_per_replicate_carriers.tex").write_text(tex)
    print(f"S4: {len(s4)} rows ({s4.module.nunique()} modules x {len(SAMPLES)} replicates)")


# ------------------------------------------------------------------------- S5
def s5_mags() -> None:
    cw = pd.read_csv(DATA / "mag_quality.tsv", sep="\t")
    gb = pd.read_csv(DATA / "mag_biosample_accessions.tsv", sep="\t")
    gb["MAG_name"] = gb.sample_name.str.replace("^MAG_", "", regex=True)
    ref = pd.read_csv(DATA / "mag_biosample_reference.tsv", sep="\t")
    ref = ref.rename(columns={"sample_name": "MAG_name"})
    d = (cw.merge(gb[["MAG_name", "biosample_accession"]], on="MAG_name", how="left")
           .merge(ref[["MAG_name", "total_length_bp", "num_contigs"]], on="MAG_name", how="left"))
    d = d.sort_values(["phylum", "MAG_name"])
    d.to_csv(HERE / "s5_mag_accessions.csv", index=False)

    rows = []
    for _, r in d.iterrows():
        name = tex_escape(r.MAG_name).replace(r"\_", r"\_\allowbreak{}")
        tax = tex_escape(r.NCBI_organism_name if pd.notna(r.NCBI_organism_name) else r.MAG_name)
        size = "--" if pd.isna(r.total_length_bp) else f"{int(r.total_length_bp):,}"
        nc = "--" if pd.isna(r.num_contigs) else f"{int(r.num_contigs):,}"
        acc = r.biosample_accession if pd.notna(r.biosample_accession) else r"\textit{pending}"
        rows.append(f"{name} & {tex_escape(r.phylum) if pd.notna(r.phylum) else 'unclassified'} & {tax} & {r.completeness:.1f} & "
                    f"{r.contamination:.1f} & {size} & {nc} & {acc} \\\\")
    head = (r"MAG & Phylum & Taxonomy (NCBI) & Compl. & Red. & Length & Contigs & BioSample \\" "\n"
            r" & & & (\%) & (\%) & (bp) & & \\")
    cols = (r"@{}>{\raggedright\arraybackslash}p{3.3cm}>{\raggedright\arraybackslash}p{2.8cm}"
            r">{\raggedright\arraybackslash}p{3.2cm}rrrrl@{}")
    tex = r"""{\scriptsize\setlength{\tabcolsep}{3pt}\setlength{\LTcapwidth}{\textwidth}
\begin{longtable}{""" + cols + r"""}
\caption{\textbf{Metagenome-assembled genomes: quality and accessions.} The
dereplicated, non-redundant MAG set recovered from the eight enrichment
co-assemblies and deposited under BioProject PRJNA1432497. Completeness and
redundancy are anvi'o single-copy-core-gene estimates. Genome (GenBank assembly)
accessions are added at proof stage; BioSample accessions are given here.}
\label{stab:mags}\\
\toprule
""" + head + r"""
\midrule
\endfirsthead
\multicolumn{8}{@{}l}{\emph{Supplementary Table~\ref{stab:mags}, continued}}\\
\toprule
""" + head + r"""
\midrule
\endhead
\bottomrule
\endlastfoot
""" + "\n".join(rows) + r"""
\end{longtable}
}
"""
    (HERE / "s5_mag_accessions.tex").write_text(tex)
    hq = (d.quality == "High Quality Draft").sum()
    print(f"S5: {len(d)} MAGs ({hq} high-quality, {len(d) - hq} medium-quality); "
          f"{d.biosample_accession.notna().sum()} with BioSample accessions")


def run_tables() -> None:
    s1_accessions()
    s4_carriers()
    s5_mags()


# ============================================================================
# Predicted lifestyle of the 88 MAGs (symclatron): run_symclatron()
# ----------------------------------------------------------------------------
# symclatron (Villada et al., Nat Biotechnol 2026) classifies a genome as
# free-living, host-associated or obligate intracellular from orthogroup content.
# It was trained on symbionts of eukaryotes, so "host-associated" here means that
# a genome resembles genomes of dependent organisms, not that it has a host. The
# question is whether the heterotrophs that carry the count-only functions also
# carry that signature.
#
# Genome completeness is one of the classifier's seven input features, and the
# MAGs range from 81.7% to 95.8% complete, so predicted class is first tested
# against completeness; a lifestyle result is reported only if it does not
# follow completeness.
#
# Inputs:
#   data/symclatron_results.tsv                 from `symclatron classify`
#   data/mag_quality.tsv                        completeness and taxonomy
#   outputs/mag_control_per_module.csv          niche class per module
# Writes outputs/symclatron_*.
# ============================================================================

RESULTS = C.DATA / "symclatron_results.tsv"


CONFIDENCE = 0.725  # the authors' optimal threshold: 1.5% false positives at 0.93 coverage


def load() -> pd.DataFrame:
    if not RESULTS.exists():
        raise SystemExit(f"missing {RESULTS}; run `symclatron classify` first")
    r = pd.read_csv(RESULTS, sep="\t")
    # the column names shift between releases, so find them rather than assume
    r = r.rename(columns={"taxon_oid": "MAG_name",
                          "classification_thresholded": "predicted_class"})
    r["MAG_name"] = r.MAG_name.astype(str).str.replace(r"\.faa$", "", regex=True)

    q = pd.read_csv(C.DATA / "mag_quality.tsv", sep="\t")
    d = q.merge(r[["MAG_name", "predicted_class", "confidence",
                   "classification", "completeness_UNI56"]],
                on="MAG_name", how="left")
    d["symbiont"] = d.predicted_class.fillna("").str.startswith("Symbiont")
    d["confident"] = d.confidence >= CONFIDENCE
    d["is_cyano"] = d.phylum == C.CYANO_PHYLUM
    return d


def run_symclatron() -> None:
    d = load()
    d.to_csv(C.OUT / "symclatron_per_mag.csv", index=False)
    n = len(d)
    print(f"{n} MAGs, {d.predicted_class.notna().sum()} classified\n")

    print("--- predicted class ---")
    print(d.predicted_class.value_counts(dropna=False).to_string())
    print(f"\nhigh confidence (>= {CONFIDENCE}): {int(d.confident.sum())} of {n}")

    # --- the control, run first ------------------------------------------------
    print("\n=== CONTROL: is predicted class confounded with completeness? ===")
    a = d.loc[d.symbiont, "completeness"]
    b = d.loc[~d.symbiont, "completeness"]
    print("(also testing completeness_UNI56, the classifier's own input feature)")
    ua = d.loc[d.symbiont, "completeness_UNI56"]
    ub = d.loc[~d.symbiont, "completeness_UNI56"]
    if len(ua) and len(ub):
        uu, up = mannwhitneyu(ua, ub, alternative="two-sided")
        print(f"UNI56 completeness: symbiont {ua.mean():.1f} vs free-living {ub.mean():.1f}, p = {up:.4f}")
    if len(a) and len(b):
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        print(f"completeness, symbiont-called   mean {a.mean():.2f} (n={len(a)})")
        print(f"completeness, free-living-called mean {b.mean():.2f} (n={len(b)})")
        print(f"two-sided Mann-Whitney U = {u:.1f}, p = {p:.4f}")
        verdict = "PASSES" if p >= 0.05 else "FAILS"
        print(f"\ncompleteness control {verdict} at alpha = 0.05.")
        if p < 0.05:
            print("Predicted lifestyle tracks assembly completeness in this set.")
            print("Report the association and do not report a lifestyle gradient.")
    else:
        print("only one predicted class present; no completeness test possible")
        p = np.nan

    # stratified view, reported either way
    d["completeness_stratum"] = pd.cut(
        d.completeness, [0, 85, 90, 95, 100],
        labels=["<85", "85-90", "90-95", ">95"]
    )
    strat = (
        d.groupby("completeness_stratum", observed=True)
        .agg(n=("MAG_name", "size"), symbiont_frac=("symbiont", "mean"))
        .round(3)
    )
    strat.to_csv(C.OUT / "symclatron_completeness_strata.csv")
    print("\n--- symbiont fraction within completeness strata ---")
    print(strat.to_string())

    # --- cyanobacteria against the heterotrophic fraction ----------------------
    tab = pd.crosstab(d.is_cyano, d.symbiont)
    print("\n--- cyanobacteria vs heterotrophs ---")
    print(tab.to_string())
    if tab.shape == (2, 2):
        odds, pf = fisher_exact(tab)
        print(f"Fisher exact odds ratio {odds:.2f}, p = {pf:.4f}")

    by_phylum = (
        d.groupby("phylum")
        .agg(n=("MAG_name", "size"), symbiont_frac=("symbiont", "mean"),
             mean_completeness=("completeness", "mean"))
        .round(3)
        .sort_values("n", ascending=False)
    )
    by_phylum.to_csv(C.OUT / "symclatron_by_phylum.csv")
    print("\n--- by phylum ---")
    print(by_phylum.to_string())

    if not np.isnan(p) and p < 0.05:
        print("\nStopping before the niche-class link: the completeness control failed,")
        print("so any association with niche class would not be interpretable.")
        sys.exit(0)

    print("\n(Next step, once the control passes: link predicted lifestyle to which")
    print("MAGs carry the community-derived modules, using mag_control_per_module.csv.)")


# ============================================================================
# Supplementary figures
# ----------------------------------------------------------------------------
# The SI figures that are drawn from data (the rest are schematics, photographs
# and anvi'o EcoPhylo exports). Each is written to outputs/si/ with the table
# behind it.
#     S2   MAG recovery and quality per STRONG co-assembly
#     S4   unassigned RpL19 coverage by taxonomic rank
#     S6   both metrics: volcanoes and PCAs, enrichment vs source
#     S7   enriched modules across MAG phyla
#     S9   module abundance, between-replicate variation and enrichment
#     S11  richness and rank-abundance across media (A, B; biomass C needs input)
# ============================================================================
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

SI = C.OUT / "si"
SI.mkdir(exist_ok=True)
SRC_C, ENR_C = "#55558f", "#8fd0b0"
MEDIA_C = {"V": "#6aae6a", "NV": "#f28468", "NVNN": "#e377b8"}
RANKS = ["phylum", "class", "order", "family", "genus", "species"]


def si_save(fig, stem, table=None):
    C.polish(fig)
    fig.savefig(SI / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(SI / f"{stem}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    if table is not None:
        table.to_csv(SI / f"{stem}.csv")


def letter(ax, s):
    ax.text(-0.13, 1.06, s, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


def si_s2_mag_recovery():
    b = pd.read_csv(DATA / "strong_bins_all.csv")
    b["group"] = b.original_sample.replace({"OSS3": "LAKE"})
    order = b.group.value_counts().index.tolist()
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.6))
    ax = axes[0, 0]
    ax.bar(order, b.group.value_counts()[order], color="#3b6e8f")
    ax.set_ylabel("Number of MAGs"); ax.set_title("MAGs per co-assembly")
    for ax, col, lab, ttl in [(axes[0, 1], "total length", "Genome length (Mbp)", "Genome length"),
                              (axes[1, 0], "completion", "SCG completeness (%)", "Completeness"),
                              (axes[1, 1], "redundancy", "SCG redundancy (%)", "Redundancy")]:
        vals = [b.loc[b.group == g, col] / (1e6 if col == "total length" else 1) for g in order]
        parts = ax.violinplot(vals, showextrema=False, widths=0.8)
        for pc in parts["bodies"]:
            pc.set_facecolor("#3b6e8f"); pc.set_alpha(0.55)
        for i, v in enumerate(vals, 1):
            ax.scatter(np.random.default_rng(i).normal(i, 0.06, len(v)), v, s=4, color="#222222")
        ax.set_xticks(range(1, len(order) + 1), order)
        ax.set_ylabel(lab); ax.set_title(ttl)
    for ax, l in zip(axes.flat, "ABCD"):
        letter(ax, l)
        ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    si_save(fig, "S2_mag_recovery", b.groupby("group")[["completion", "redundancy", "total length"]].describe())


def si_s4_unassigned():
    """Fraction of RpL19 coverage without a GTDB name at each rank, over all clusters
    (clusters with no GTDB record count as unnamed at every rank). GTDB output only."""
    cov = pd.read_csv(DATA / "rpl19_cluster_coverage.csv", index_col=0)
    rel = cov / cov.sum()
    gtdb = pd.read_csv(DATA / "rpl19_cluster_taxonomy.tsv", sep="\t", index_col="split_name").reindex(rel.index)
    out = {grp: pd.DataFrame({r: rel.loc[gtdb[f"t_{r}"].isna(), members].sum() for r in RANKS})
           for grp, members in (("Source", C.SOURCE), ("Enrichment", C.ENRICHMENT))}
    fig, ax = plt.subplots(figsize=(4.6, 2.8))
    x = np.arange(len(RANKS))
    for k, (grp, col) in enumerate((("Source", SRC_C), ("Enrichment", ENR_C))):
        g = out[grp]
        xo = x + (k - 0.5) * 0.38
        ax.bar(xo, g.median(), 0.38, color=col, edgecolor="#333333", alpha=0.85, lw=0.6,
               label=f"{'Reservoir' if grp == 'Source' else 'Enrichments'} (n = {len(g)})")
        for j, r in enumerate(RANKS):
            ax.scatter(xo[j] + np.random.default_rng(j).normal(0, 0.04, len(g)), g[r], s=7,
                       color="black", lw=0, zorder=3)
            if grp == "Enrichment":
                last = -1.0
                for smp, v in g[r][g[r] > 0.2].sort_values().items():   # labels kept 0.04 apart
                    last = max(v, last + 0.04)
                    ax.text(xo[j] + 0.16, last, smp.replace("_", " "), fontsize=5.5, va="center")
    ax.set_xticks(x, [r.capitalize() for r in RANKS], fontsize=7)
    ax.set_ylim(0, 1); ax.tick_params(labelsize=7)
    ax.set_ylabel("Fraction of RpL19 coverage\nwithout a GTDB name", fontsize=7)
    ax.legend(frameon=False, fontsize=6.5, loc="upper left", title="bars, median; points, samples",
              title_fontsize=6.5, alignment="left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    tab = pd.concat(out, names=["group", "sample"])
    si_save(fig, "S4_unassigned_coverage", tab)
    return tab


def si_s6_both_metrics():
    import function as enr
    t, cov_m, _ = enr.two_metric("rpl19")
    cov_m, cnt_m, _ = C.metric_matrices("enrichments", "rpl19")
    fig, axes = plt.subplots(2, 2, figsize=(10, 8.2))
    for ax, m, lab in ((axes[0, 0], "cnt", "count-based"), (axes[0, 1], "cov", "coverage-based")):
        d = t.dropna(subset=[f"{m}_l2fc"])
        y = -np.log10(d[f"{m}_padj"])
        sig = (d[f"{m}_padj"] < C.FDR) & (d[f"{m}_l2fc"].abs() > C.MIN_L2FC)
        up, dn = sig & (d[f"{m}_l2fc"] > 0), sig & (d[f"{m}_l2fc"] < 0)
        ax.scatter(d[f"{m}_l2fc"][~sig], y[~sig], s=8, color="#bbbbbb", lw=0, label="not significant")
        ax.scatter(d[f"{m}_l2fc"][up], y[up], s=16, color=ENR_C, edgecolor="#333333", lw=0.3,
                   label=f"higher in enrichment ({int(up.sum())})")
        ax.scatter(d[f"{m}_l2fc"][dn], y[dn], s=16, color=SRC_C, edgecolor="#333333", lw=0.3,
                   label=f"higher in source ({int(dn.sum())})")
        ax.axhline(-np.log10(C.FDR), color="#e0605e", ls="--", lw=0.9)
        for v in (-C.MIN_L2FC, C.MIN_L2FC):
            ax.axvline(v, color="#aaaaaa", ls="--", lw=0.8)
        ax.set_xlabel(f"{lab.capitalize()} log$_2$FC (enrichment / source)")
        ax.set_ylabel("$-$log$_{10}$ $q$")
        ax.set_title(lab.upper())
        ax.legend(fontsize=7, frameon=False, loc="upper left")
    samples = C.SOURCE + C.ENRICHMENT
    for ax, mat, lab in ((axes[1, 0], cnt_m, "count-based"), (axes[1, 1], C.compositional(cov_m), "coverage-based")):
        X = StandardScaler().fit_transform(mat[samples].T.values)
        p = PCA(2).fit(X); sc = p.transform(X); v = p.explained_variance_ratio_ * 100
        F, pv = C.permanova(mat, {"S": C.SOURCE, "E": C.ENRICHMENT})
        for grp, col, idx in (("Source", SRC_C, range(8)), ("Enrichment", ENR_C, range(8, 16))):
            ax.scatter(sc[list(idx), 0], sc[list(idx), 1], s=40, color=col, edgecolor="#333333", lw=0.5, label=grp)
        ax.set_xlabel(f"PC1 ({v[0]:.1f}%)"); ax.set_ylabel(f"PC2 ({v[1]:.1f}%)")
        ax.set_title(f"{lab.upper()} (z-scored)")
        ax.text(0.98, 0.03, f"PERMANOVA F = {F:.1f}, p = {pv:.3f}", transform=ax.transAxes, ha="right", fontsize=7.5)
        ax.legend(fontsize=7, frameon=False)
    for ax, l in zip(axes.flat, "ABCD"):
        letter(ax, l)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    si_save(fig, "S6_both_metrics", t[["cov_l2fc", "cov_padj", "cnt_l2fc", "cnt_padj", "niche_class", "module_name"]])


def si_s7_mag_phyla():
    comp, cov, tax = C.mag_tables()
    cls = pd.read_csv(C.OUT / "module_niche_classification.csv", index_col=0)
    col = "niche_class" if "niche_class" in cls.columns else cls.columns[0]
    mods = [m for m in cls.index[cls[col].isin(["both metrics", "count-based only", "coverage-based only"])] if m in comp.index]
    enc = comp.loc[mods] >= C.ENCODING
    phyla = tax.phylum.value_counts()
    phyla = phyla[phyla >= 2].index
    frac = pd.DataFrame({ph: enc.loc[:, tax.phylum == ph].mean(axis=1) * 100 for ph in phyla}).T
    spread = (frac > 0).mean()
    order = spread.sort_values().index
    frac = frac[order].loc[frac[order].gt(0).sum(axis=1).sort_values(ascending=False).index]
    fig = plt.figure(figsize=(13, 5.2))
    ax = fig.add_axes([0.22, 0.30, 0.68, 0.50])
    im = ax.imshow(frac.values, aspect="auto", cmap="Greys", vmin=0, vmax=100)
    ax.set_yticks(range(len(frac)), [f"{p} ({(tax.phylum == p).sum()})" for p in frac.index], fontsize=7.5)
    names = cls.reindex(frac.columns).get("module_name", pd.Series(frac.columns, index=frac.columns))
    ax.set_xticks(range(frac.shape[1]), [str(n).split(",")[0][:38] for n in names], rotation=90, fontsize=4.6)
    top = fig.add_axes([0.22, 0.82, 0.68, 0.10], sharex=ax)
    top.bar(range(frac.shape[1]), spread[order], color="#777777", width=0.8)
    top.set_ylabel("fraction of\nphyla", fontsize=7); top.tick_params(labelbottom=False, labelsize=6)
    top.set_title("Enriched modules (either metric), ordered from sporadic to widespread", fontsize=9)
    side = fig.add_axes([0.03, 0.30, 0.06, 0.50], sharey=ax)
    side.barh(range(len(frac)), (frac > 0).sum(axis=1), color="#777777")
    side.invert_xaxis(); side.tick_params(labelleft=False, labelsize=6)
    side.set_xlabel("# modules\nin phylum", fontsize=7)
    cax = fig.add_axes([0.92, 0.30, 0.012, 0.50])
    fig.colorbar(im, cax=cax).set_label("% of MAGs in phylum encoding module", fontsize=7)
    for a in (top, side):
        a.spines[["top", "right"]].set_visible(False)
    si_save(fig, "S7_mag_phyla_modules", frac)


def si_s9_cv():
    import function as enr
    t, cov_m, _ = enr.two_metric("rpl19")
    e = cov_m[[s for s in C.ENRICHMENT if s in cov_m.columns]]
    d = t.assign(mean_cov=e.mean(axis=1)).dropna(subset=["cov_cv", "cov_l2fc"])
    d = d[d.mean_cov > 0]
    sig = (d.cov_padj < C.FDR) & (d.cov_l2fc.abs() > C.MIN_L2FC)
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    norm = mpl.colors.TwoSlopeNorm(vcenter=0, vmin=min(d.cov_l2fc.min(), -1), vmax=d.cov_l2fc.max())
    sc = ax.scatter(d.mean_cov[~sig], d.cov_cv[~sig], c=d.cov_l2fc[~sig], cmap="RdYlGn", norm=norm,
                    s=22, edgecolor="#555555", lw=0.3, label="not significant")
    ax.scatter(d.mean_cov[sig], d.cov_cv[sig], c=d.cov_l2fc[sig], cmap="RdYlGn", norm=norm,
               s=60, marker="X", edgecolor="#333333", lw=0.4, label="FDR < 0.05")
    ax.set_xscale("log")
    ax.set_xlabel("Mean module coverage across enrichments")
    ax.set_ylabel("Coefficient of variation across enrichments")
    fig.colorbar(sc, ax=ax).set_label("Coverage-based log$_2$FC (enrichment / source)")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    si_save(fig, "S9_cv_vs_abundance", d[["module_name", "mean_cov", "cov_cv", "cov_l2fc", "cov_padj"]])


def si_s11_media_richness():
    e = pd.read_csv(DATA / "ecophylo_coverage.csv")
    e = e[e.run == "rps16_media_endpoints"]
    m = e.pivot_table(index="cluster", columns="sample", values="mean_coverage_Q2Q3").fillna(0)
    groups = {g: [s for s in m.columns if s.startswith(g + "_P10")] for g in ("V", "NV", "NVNN")}
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))
    rich = (m > 0).sum()
    for i, (g, cols) in enumerate(groups.items()):
        axes[0].boxplot(rich[cols], positions=[i], widths=0.5, patch_artist=True,
                        boxprops=dict(facecolor=MEDIA_C[g], alpha=0.5), medianprops=dict(color="#333333"))
        axes[0].scatter(np.random.default_rng(i).normal(i, 0.05, len(cols)), rich[cols], color=MEDIA_C[g],
                        edgecolor="#333333", lw=0.4, zorder=3)
        rel = m[cols] / m[cols].sum()
        ranked = pd.DataFrame({c: np.sort(rel[c].values)[::-1] for c in cols})
        ranked = ranked[(ranked > 0).any(axis=1)]
        for c in cols:
            axes[1].plot(np.arange(1, len(ranked) + 1), ranked[c].replace(0, np.nan), color=MEDIA_C[g], lw=0.6, alpha=0.4)
        axes[1].plot(np.arange(1, len(ranked) + 1), ranked.mean(axis=1).replace(0, np.nan), color=MEDIA_C[g], lw=2, label=g)
    axes[0].set_xticks(range(3), list(groups)); axes[0].set_ylabel("RpS16 populations detected")
    axes[1].set_yscale("log"); axes[1].set_xlabel("Population rank"); axes[1].set_ylabel("Relative coverage")
    axes[1].legend(frameon=False)
    for ax, l in zip(axes, "AB"):
        letter(ax, l); ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    si_save(fig, "S11AB_media_richness_rank", rich.rename("rps16_populations").to_frame())


def run_si_figures() -> None:
    for f in (si_s2_mag_recovery, si_s4_unassigned, si_s6_both_metrics, si_s7_mag_phyla, si_s9_cv,
              si_s11_media_richness):
        f()
        print(f"SI: {f.__name__}")


# ============================================================================
# Ordination robustness (SI): alternative distances on every dataset
# ----------------------------------------------------------------------------
# The main figures ordinate with Bray-Curtis on relative coverage (common.bray_pcoa).
# Here the same samples are ordinated with the alternatives a reader might expect:
# PCA on z-scored or log relative coverage (a PCA is a PCoA on Euclidean distances),
# Bray-Curtis on raw coverage (depth and biomass included), weighted and unweighted
# UniFrac on the marker tree for composition, and the count-based metric for function.
# Each panel gives the PERMANOVA p and the mean distance of each group's samples to
# their centroid. Separation holds under every measure; group spread does not, and it
# is the reason the main figures use one measure throughout.
#     outputs/si/S_ordination_robustness.{pdf,csv}
# ============================================================================
def _media_raw_rps16() -> pd.DataFrame:
    e = pd.read_csv(C.DATA / "ecophylo_coverage.csv")
    e = e[e.run == "rps16_media_endpoints"]
    return e.pivot_table(index="cluster", columns="sample", values="mean_coverage_Q2Q3").fillna(0)


def _tree(path: Path, tips: list[str]):
    from skbio import TreeNode
    t = TreeNode.read([path.read_text()])
    for tip in t.tips():
        tip.name = tip.name.replace(" ", "_")
    t = t.shear([x for x in tips if x in {n.name for n in t.tips()}]).root_at_midpoint()
    for n in t.traverse():
        if n.length is None:
            n.length = 0.0
    return t


def _distances(raw: pd.DataFrame, samples: list[str], tree=None) -> dict:
    """name -> skbio DistanceMatrix for one dataset (features x samples, raw values)."""
    from skbio import DistanceMatrix
    from skbio.diversity import beta_diversity
    from scipy.spatial.distance import pdist, squareform
    from sklearn.preprocessing import StandardScaler
    raw = raw[samples].fillna(0.0)
    raw = raw[raw.sum(axis=1) > 0]
    rel = raw / raw.sum()
    euc = lambda X: DistanceMatrix(squareform(pdist(X)), ids=samples)
    out = {"PCA, z-scored": euc(StandardScaler().fit_transform(rel.T.values)),
           "PCA, log relative": euc(np.log10(1 + 1e4 * rel.T.values)),
           "Bray-Curtis, relative": beta_diversity("braycurtis", rel.T.values, ids=samples),
           "Bray-Curtis, raw": beta_diversity("braycurtis", raw.T.values, ids=samples)}
    if tree is not None:
        keep = [x for x in rel.index if x in {t.name for t in tree.tips()}]
        for name, metric, kw in (("weighted UniFrac", "weighted_unifrac", {"normalized": True}),
                                 ("unweighted UniFrac", "unweighted_unifrac", {})):
            out[name] = beta_diversity(metric, raw.loc[keep].T.values, ids=samples, taxa=keep,
                                       tree=tree, **kw)
    return out


def ordination_robustness():
    import matplotlib.pyplot as plt
    import media as med
    import taxonomy as tx
    from figures import ellipse, CI95, MEDIA
    enr = {"Source": C.SOURCE, "Enrichment": C.ENRICHMENT}
    cov_m, cnt_m, _ = C.metric_matrices("enrichments", "rpl19")
    _, mcov, mcnt, _ = med.contrasts("rps16")
    mraw = _media_raw_rps16()
    mgr = {g: [x for x in med.GROUPS[g] if x in mraw.columns] for g in ("V", "NV", "NVNN")}
    mgf = {g: [x for x in med.GROUPS[g] if x in mcov.columns] for g in ("V", "NV", "NVNN")}
    tree_e, cov_e = tx.load()
    s_e = C.SOURCE + C.ENRICHMENT
    s_m = sum(mgr.values(), [])
    tree_m = _tree(C.DATA / "ecophylo_trees" / "rps16_media_endpoints.nwk", list(mraw.index))
    rows_def = [
        ("Source vs enrichment: composition (RpL19)", _distances(cov_e, s_e, tree_e), enr),
        ("Source vs enrichment: function", {**_distances(cov_m, s_e),
             "count-based, Bray-Curtis": _distances(cnt_m, s_e)["Bray-Curtis, relative"]}, enr),
        ("Media: composition (RpS16)", _distances(mraw, s_m, tree_m), mgr),
        ("Media: function", {**_distances(mcov, sum(mgf.values(), [])),
             "count-based, Bray-Curtis": _distances(mcnt, sum(mgf.values(), []))["Bray-Curtis, relative"]}, mgf),
    ]
    color = {"Source": C.GROUP_COLOR["Source"], "Enrichment": C.GROUP_COLOR["Enrichment"], **MEDIA}
    short = {"Source": "S", "Enrichment": "E"}
    nrow = max(len(d) for _, d, _ in rows_def)
    # datasets as columns, measures as rows: a full portrait page with four wide columns
    fig, axes = plt.subplots(nrow, len(rows_def), figsize=(7.2, 9.4))
    fig.subplots_adjust(left=0.05, right=0.975, top=0.92, bottom=0.05, wspace=0.24, hspace=0.62)
    out = []
    for j, (title, dists, groups) in enumerate(rows_def):
        for i in range(nrow):
            ax = axes[i, j]
            if i >= len(dists):
                ax.axis("off")
                ax.text(0.5, 0.5, "UniFrac not defined\nfor functions", ha="center", va="center",
                        fontsize=6.5, color="#999999", transform=ax.transAxes)
                continue
            name, d = list(dists.items())[i]
            r = C.pcoa_stats(d, groups)
            for g, members in groups.items():
                xy = r["scores"].loc[[m for m in members if m in r["scores"].index]].values
                ellipse(ax, xy, color[g], n_std=CI95)
                ax.scatter(xy[:, 0], xy[:, 1], s=13, color=color[g], edgecolor="#333333", lw=0.3, zorder=3)
            main = "count-based, Bray-Curtis" if title == "Media: function" else "Bray-Curtis, relative"
            ax.set_title(name.replace("count-based, ", "count-based, ") + (" (main)" if name == main else ""),
                         fontsize=6.5,
                         fontweight="bold" if name == main else "normal")
            disp = ", ".join(f"{short.get(g, g)} {v:.1f}" if v >= 1 else f"{short.get(g, g)} {v:.2f}"
                             for g, v in r["centroid_dist"].items())
            pp = "< 0.001" if r["p"] < 0.001 else f"= {r['p']:.3f}"
            ax.text(0.02, 0.02, f"PERMANOVA p {pp}\nspread {disp}", transform=ax.transAxes, fontsize=5.6,
                    va="bottom", color="#444444")
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_xlabel(f"Axis 1 ({r['explained'][0]:.0f}%)", fontsize=6, labelpad=1)
            ax.set_ylabel(f"Axis 2 ({r['explained'][1]:.0f}%)", fontsize=6, labelpad=1)
            ax.margins(x=0.15, y=0.45)
            ax.spines[["top", "right"]].set_visible(False)
            if i == 0:
                head = title.replace(": ", ":\n")
                ax.text(0.5, 1.22, head, transform=ax.transAxes, ha="center", va="bottom",
                        fontsize=7.5, fontweight="bold", clip_on=False)
            out.append(dict(dataset=title, distance=name, PCo1=r["explained"][0], PCo2=r["explained"][1],
                            F=r["F"], p=r["p"], disp_p=r["disp_p"],
                            **{f"spread_{g}": v for g, v in r["centroid_dist"].items()}))
    (C.OUT / "si").mkdir(exist_ok=True)
    hs = [plt.Line2D([], [], ls="", marker="o", ms=5, color=c, mec="#333333", mew=0.3, label=g)
          for g, c in color.items()]
    fig.legend(handles=hs, loc="lower center", ncol=5, fontsize=7, frameon=False,
               bbox_to_anchor=(0.5, 0.0))
    C.polish(fig)
    fig.savefig(C.OUT / "si" / "S_ordination_robustness.pdf")
    fig.savefig(C.OUT / "si" / "S_ordination_robustness.png", dpi=200)
    plt.close(fig)
    t = pd.DataFrame(out)
    t.to_csv(C.OUT / "si" / "S_ordination_robustness.csv", index=False)
    return t


# ============================================================================
# Ordination at every taxonomic rank (SI)
# ----------------------------------------------------------------------------
# Bray-Curtis PCoA of reservoir and enrichment samples at phylum, class, order,
# family, genus and population (97% RpL19), with PERMANOVA, PERMDISP and each
# group's spread (mean distance to centroid). Same tables as Fig 2C.
#     outputs/si/S_pcoa_by_rank.{pdf,png,csv}
# ============================================================================
def pcoa_by_rank():
    import matplotlib.pyplot as plt
    import taxonomy as tx
    from figures import ellipse, CI95, fmt_p
    lt = tx.level_tables()
    levels = ["phylum", "class", "order", "family", "genus", "population", "function, coverage",
              "function, count"]
    titles = {"population": "Population (97% RpL19)", "function, coverage": "Function, coverage-based",
              "function, count": "Function, count-based"}
    groups = {"Source": C.SOURCE, "Enrichment": C.ENRICHMENT}
    fig, axes = plt.subplots(2, 4, figsize=(7.2, 4.3))
    fig.subplots_adjust(left=0.07, right=0.99, top=0.86, bottom=0.09, wspace=0.42, hspace=0.6)
    axes[0, 0].text(-0.3, 1.25, "C", transform=axes[0, 0].transAxes, fontsize=10, fontweight="bold")
    rows = []
    for ax, lev in zip(axes.flat, levels):
        r = C.bray_pcoa(lt[lev], groups)
        for g, members in groups.items():
            xy = r["scores"].loc[members].values
            ellipse(ax, xy, C.GROUP_COLOR[g], n_std=CI95)
            ax.scatter(xy[:, 0], xy[:, 1], s=14, marker="o" if g == "Enrichment" else "s",
                       color=C.GROUP_COLOR[g], edgecolor="#333333" if g == "Enrichment" else "white",
                       lw=0.4, zorder=3)
        cd = r["centroid_dist"]
        ax.set_title(titles.get(lev, lev.capitalize()), fontsize=7, fontweight="bold")
        ax.text(0.02, 0.02, f"spread S {cd['Source']:.2f}, E {cd['Enrichment']:.2f}"
                f"\nPERMDISP {fmt_p(r['disp_p'])}", transform=ax.transAxes, fontsize=6,
                color="#444444", va="bottom")
        ax.set_xlabel(f"PCo1 ({r['explained'][0]:.0f}%)", fontsize=6.5)
        ax.set_ylabel(f"PCo2 ({r['explained'][1]:.0f}%)", fontsize=6.5)
        ax.set_xticks([]); ax.set_yticks([])
        ax.margins(0.2)
        ax.spines[["top", "right"]].set_visible(False)
        rows.append(dict(level=lev, PCo1=r["explained"][0], PCo2=r["explained"][1], F=r["F"], p=r["p"],
                         disp_p=r["disp_p"], spread_source=cd["Source"], spread_enrichment=cd["Enrichment"]))
    hs = [plt.Line2D([], [], ls="", marker="o" if g == "Enrichment" else "s", ms=5, color=C.GROUP_COLOR[g],
                     mec="#333333" if g == "Enrichment" else "white", mew=0.4, label=g) for g in groups]
    fig.legend(handles=hs, loc="upper center", ncol=2, fontsize=7, frameon=False, bbox_to_anchor=(0.5, 1.0))
    (C.OUT / "si").mkdir(exist_ok=True)
    C.polish(fig)
    fig.savefig(C.OUT / "si" / "S_pcoa_by_rank.pdf")
    fig.savefig(C.OUT / "si" / "S_pcoa_by_rank.png", dpi=200)
    plt.close(fig)
    t = pd.DataFrame(rows)
    t.to_csv(C.OUT / "si" / "S_pcoa_by_rank.csv", index=False)
    return t


# ============================================================================
# Cyanobacterial morphotypes in every sample (SI)
# ----------------------------------------------------------------------------
# Long-term enrichments and reservoir (RpL19, final run) and media P10 cultures
# (RpS16). Cyanobacterial clusters are grouped as picocyanobacteria (family
# Cyanobiaceae), filamentous (families of filamentous genera), other cyanobacteria,
# or unresolved: placed inside the cyanobacteria on the marker tree but further than
# PLACE_MAX from any named tip, which cannot be told apart from a plastid homologue
# (in NV_P10_R2-R4 and V_P10_R3 reads map to no cyanobacterial genome, pipeline/09).
#     outputs/si/S_morphotypes.{pdf,png,csv}
# ============================================================================
FILAMENTOUS = {"Leptolyngbyaceae", "Pseudanabaenaceae", "Microcoleaceae", "Phormidesmiaceae",
               "Nostocaceae", "Elainellaceae", "JAAUUE01"}
MORPH_COL = {"Picocyanobacteria (Cyanobiaceae)": "#1b9e77", "Filamentous": "#d95f02",
             "Other cyanobacteria": "#7570b3", "Unresolved (cyanobacteria or plastid)": "#bdbdbd"}


def _morph(tax, place_max):
    cy = tax.t_phylum == C.CYANO_PHYLUM
    far = tax.placed & (tax.patristic > place_max)
    return pd.Series(np.select(
        [~cy, cy & far, cy & (tax.t_family == "Cyanobiaceae"), cy & tax.t_family.isin(FILAMENTOUS)],
        ["none", "Unresolved (cyanobacteria or plastid)", "Picocyanobacteria (Cyanobiaceae)", "Filamentous"],
        default="Other cyanobacteria"), index=tax.index)


def morphotypes():
    import matplotlib.pyplot as plt
    import taxonomy as tx
    from figures import PLACE_MAX, MEDIA
    import media as med
    rel = tx.load_clusters(relative=True)
    m1 = rel.groupby(_morph(C.rpl19_taxonomy().reindex(rel.index), PLACE_MAX)).sum()
    e = pd.read_csv(C.DATA / "ecophylo_coverage.csv"); e = e[e.run == "rps16_media_endpoints"]
    m = e.pivot_table(index="cluster", columns="sample", values="mean_coverage_Q2Q3").fillna(0)
    r2 = m / m.sum()
    m2 = r2.groupby(_morph(C.ecophylo_taxonomy("rps16_media_endpoints").reindex(r2.index), PLACE_MAX)).sum()
    blocks = [("Reservoir", C.SOURCE, m1), ("Long-term enrichments", C.ENRICHMENT, m1)]
    blocks += [(g, [f"{g}_P10_R{i}" for i in range(1, 6) if f"{g}_P10_R{i}" in m2.columns], m2)
               for g in ("V", "NV", "NVNN")]
    fig, ax = plt.subplots(figsize=(7.2, 2.8))
    x, ticks, rows = 0, [], []
    for name, samples, mat in blocks:
        x0 = x
        for smp in samples:
            bottom = 0.0
            for cat, col in MORPH_COL.items():
                v = mat.loc[cat, smp] if cat in mat.index else 0.0
                ax.bar(x, v, bottom=bottom, color=col, width=0.85, lw=0)
                bottom += v
            rows.append(dict(group=name, sample=smp, **{c: (mat.loc[c, smp] if c in mat.index else 0.0)
                                                         for c in MORPH_COL}))
            ticks.append((x, smp)); x += 1
        ax.text((x0 + x - 1) / 2, 1.02, name, ha="center", va="bottom", fontsize=7,
                color=MEDIA.get(name, "#333333"), fontweight="bold")
        x += 0.8
    ax.set_xticks([t[0] for t in ticks], [t[1].replace("_P10", "") for t in ticks], rotation=90, fontsize=5.5)
    ax.set_ylim(0, 1); ax.set_ylabel("Share of marker coverage")
    ax.spines[["top", "right"]].set_visible(False)
    hs = [plt.Rectangle((0, 0), 1, 1, color=c, label=k) for k, c in MORPH_COL.items()]
    ax.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=4, fontsize=6, frameon=False)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.88, bottom=0.36)
    (C.OUT / "si").mkdir(exist_ok=True)
    C.polish(fig)
    fig.savefig(C.OUT / "si" / "S_morphotypes.pdf"); fig.savefig(C.OUT / "si" / "S_morphotypes.png", dpi=200)
    plt.close(fig)
    t = pd.DataFrame(rows); t.to_csv(C.OUT / "si" / "S_morphotypes.csv", index=False)
    return t


# ============================================================================
# Media: who carries each functional category (SI)
# ----------------------------------------------------------------------------
#     outputs/si/S_media_category_partition.{pdf,png}; numbers in media.category_partition*
# ============================================================================
def media_category_partition_figure():
    """Arrow chart per passage (P4: all 15 lineages; P10: lineages that still carry
    cyanobacteria): for each functional category, log2 of the fraction of cyanobacterial over
    the fraction of heterotrophic genome coverage that encodes it (media.category_bias), with
    nitrate (V + NV) and without (NVNN)."""
    import matplotlib.pyplot as plt
    import media as med
    from figures import MEDIA
    from statsmodels.stats.multitest import multipletests
    cats = med.CATEGORY_ORDER
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True)
    rng = np.random.default_rng(0)
    out = {}
    for ax, ps, ttl in ((axes[0], "P4", "P4: all lineages (5 without, 10 with nitrate)"),
                        (axes[1], "P10", "P10: lineages with cyanobacteria (4 without, 6 with)")):
        long, t = med.category_bias(passage=ps)
        t["q"] = multipletests(t.p, method="fdr_bh")[1]
        out[ps] = t
        for i, c in enumerate(cats):
            y = len(cats) - 1 - i
            g = long[long.category == c]
            rep, nv = g[g.medium != "NVNN"].log2_bias, g[g.medium == "NVNN"].log2_bias
            ax.scatter(rep, y + 0.13 + rng.uniform(-0.04, 0.04, len(rep)), s=7, color="#7a8796", lw=0, zorder=2)
            ax.scatter(nv, y - 0.13 + rng.uniform(-0.04, 0.04, len(nv)), s=7, color=MEDIA["NVNN"], lw=0,
                       alpha=0.75, zorder=2)
            sig = abs(t.loc[c, "nvnn"] - t.loc[c, "replete"]) > 1      # black: shift > 1 log2 unit
            ax.annotate("", xy=(t.loc[c, "nvnn"], y), xytext=(t.loc[c, "replete"], y), zorder=4,
                        arrowprops=dict(arrowstyle="-|>", lw=2.0 if sig else 1.0, mutation_scale=9,
                                        color="black" if sig else "#9a9a9a"))
            q = t.loc[c, "q"]
            ax.text(1.01, y, f"{q:.3f}", transform=ax.get_yaxis_transform(), ha="left", va="center",
                    fontsize=6, fontweight="bold" if q < 0.05 else "normal")
        ax.text(1.01, len(cats) - 0.35, "$q$", transform=ax.get_yaxis_transform(), ha="left", fontsize=6)
        ax.axvline(0, color="#bbbbbb", lw=0.8, zorder=1)
        ax.set_xlim(-4.2, 4.2)
        ax.set_title(ttl, fontsize=7, fontweight="bold", pad=12)
        ax.text(0.01, 1.01, "\u2190 heterotrophic", transform=ax.transAxes, fontsize=6, ha="left")
        ax.text(0.99, 1.01, "cyanobacterial \u2192", transform=ax.transAxes, fontsize=6, ha="right")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
    axes[0].set_yticks(range(len(cats)), cats[::-1], fontsize=7)
    fig.supxlabel("log$_2$ (fraction of cyanobacterial / heterotrophic genome coverage encoding the category)",
                  fontsize=7, y=0.12)
    hs = [plt.Line2D([], [], ls="", marker="o", ms=4, color="#7a8796", label="with nitrate (V, NV)"),
          plt.Line2D([], [], ls="", marker="o", ms=4, color=MEDIA["NVNN"], label="without nitrate (NVNN)"),
          plt.Line2D([], [], color="black", lw=1.5, label="mean, with \u2192 without nitrate")]
    fig.legend(handles=hs, loc="lower center", ncol=3, fontsize=6, frameon=False, handletextpad=0.3)
    fig.subplots_adjust(left=0.17, right=0.95, top=0.88, bottom=0.25, wspace=0.22)
    C.polish(fig)
    fig.savefig(C.OUT / "si" / "S_media_category_partition.pdf")
    fig.savefig(C.OUT / "si" / "S_media_category_partition.png", dpi=220)
    plt.close(fig)
    return pd.concat(out, names=["passage"])


def media_pairwise_dissimilarity():
    """S19 (former Fig 4E): Bray-Curtis dissimilarity between every pair of P10 media
    cultures for composition (RpS16) and function (both metrics, relative), within each
    medium and between media, with the PERMANOVA R2 of the medium."""
    import matplotlib.pyplot as plt
    import media as med
    from figures import f4_convergence, rps16_matrix, small_style
    small_style()
    _, mcov, mcnt, _ = med.contrasts("rps16")
    fig, ax = plt.subplots(figsize=(7.2, 2.6))
    conv, stats = f4_convergence(ax, rps16_matrix(), mcov, mcnt)
    fig.subplots_adjust(left=0.09, right=0.99, top=0.80, bottom=0.17)
    C.polish(fig)
    fig.savefig(C.OUT / "si" / "S_media_pairwise.pdf")
    fig.savefig(C.OUT / "si" / "S_media_pairwise.png", dpi=220)
    plt.close(fig)
    conv.to_csv(C.OUT / "si" / "S_media_pairwise.csv")
    return conv, stats


# ============================================================================
# Convergence under every measure (SI; replaces the 24-panel ordination grid)
# ----------------------------------------------------------------------------
# Convergence is read as the similarity of replicate enrichments to one another,
# against the similarity of reservoir samples to one another (Fig 2C), not as the
# distance of the enrichments from the reservoir. Panel A repeats Fig 2C under three
# dissimilarities that weight populations differently (Bray-Curtis: abundance;
# Hellinger: square-root, damps the dominants; Jaccard: presence at >= 0.1%, no
# abundance), with the matching UniFrac on the marker tree (weighted with
# Bray-Curtis, unweighted with Jaccard), and with RpS16 in place of RpL19; panel B
# does the same for function. NMDS was tried and is degenerate on these data (stress
# ~0, each group collapsed to a point) because nearly every reservoir-enrichment
# distance exceeds every within-group distance, so it is not shown.
#     outputs/si/S_convergence_robustness.{pdf,png,csv}
#     outputs/si/S_ordination_robustness_table.tex
# ============================================================================
CONV_LEVELS = ["phylum", "class", "order", "family", "genus", "species", "population"]
DETECT = 0.001


def _rps16_level_tables() -> tuple[dict, pd.DataFrame]:
    e = pd.read_csv(C.DATA / "ecophylo_coverage.csv")
    e = e[e.run == "rps16_source_enrichment"]
    raw = e.pivot_table(index="cluster", columns="sample", values="mean_coverage_Q2Q3").fillna(0)
    raw = raw[C.SOURCE + C.ENRICHMENT]
    raw = raw[raw.sum(axis=1) > 0]
    rel = raw / raw.sum()
    tax = pd.read_csv(C.DATA / "ecophylo_clusters.csv")
    tax = tax[tax.run == "rps16_source_enrichment"].set_index("cluster").reindex(rel.index)
    own = pd.Series(rel.index, index=rel.index)
    out = {"phylum": rel.groupby(tax.t_phylum.fillna("unassigned")).sum()}
    for r in CONV_LEVELS[1:-1]:
        out[r] = rel.groupby(tax[f"t_{r}"].where(tax[f"t_{r}"].notna(), own)).sum()
    out["population"] = rel
    return out, raw


def _pair_d(m: pd.DataFrame, members: list, measure: str) -> np.ndarray:
    from scipy.spatial.distance import pdist
    X = m[members].T.values.astype(float)
    if measure == "braycurtis":
        return pdist(X, "braycurtis")
    if measure == "hellinger":
        return pdist(np.sqrt(X), "euclidean") / np.sqrt(2)
    if measure == "jaccard":
        return pdist(X >= DETECT, "jaccard")
    raise ValueError(measure)


def convergence_robustness():
    import itertools
    import matplotlib.pyplot as plt
    import taxonomy as tx
    from skbio.diversity import beta_diversity
    from figures import ellipse, CI95, F2_MARKER, small_style
    small_style()
    lt = tx.level_tables()
    lt16, raw16 = _rps16_level_tables()
    tree16 = _tree(C.DATA / "ecophylo_trees" / "rps16_source_enrichment.nwk", list(raw16.index))
    keep = [x for x in raw16.index if x in {t.name for t in tree16.tips()}]
    uf16 = {k: beta_diversity(k, raw16.loc[keep].T.values, ids=list(raw16.columns), taxa=keep,
                              tree=tree16, **kw).to_data_frame()
            for k, kw in (("weighted_unifrac", {"normalized": True}), ("unweighted_unifrac", {}))}
    uf19 = {k: pd.read_csv(C.OUT / f"phylo_{k}_distances.csv", index_col=0)
            for k in ("weighted_unifrac", "unweighted_unifrac")}
    panels = [("Bray-Curtis (abundance)", lt, "braycurtis", uf19["weighted_unifrac"], "weighted\nUniFrac"),
              ("Hellinger (damps dominants)", lt, "hellinger", None, None),
              ("Jaccard (presence only)", lt, "jaccard", uf19["unweighted_unifrac"], "unweighted\nUniFrac"),
              ("RpS16, Bray-Curtis", lt16, "braycurtis", uf16["weighted_unifrac"], "weighted\nUniFrac")]
    groups = {"Source": C.SOURCE, "Enrichment": C.ENRICHMENT}
    xlab = ["Phylum", "Class", "Order", "Family", "Genus", "Species", "Population"]

    fig = plt.figure(figsize=(7.2, 2.9))
    gs = fig.add_gridspec(1, 5, left=0.075, right=0.995, top=0.88, bottom=0.25, wspace=0.34,
                          width_ratios=[1, 1, 1, 1, 0.95])
    rows = []
    for j, (title, tabs, meas, uf, uflab) in enumerate(panels):
        ax = fig.add_subplot(gs[0, j])
        xs = np.arange(len(CONV_LEVELS))
        for g, members in groups.items():
            col = C.GROUP_COLOR[g]
            mk = dict(marker=F2_MARKER[g], ms=3.4, mec="#333333" if g == "Enrichment" else "white", mew=0.4)
            med = []
            for x, lev in zip(xs, CONV_LEVELS):
                s = 1 - _pair_d(tabs[lev], members, meas)
                med.append(np.median(s))
                rows.append(dict(panel=title, level=lev, group=g, median_similarity=np.median(s),
                                 q25=np.quantile(s, 0.25), q75=np.quantile(s, 0.75)))
                ax.plot([x, x], np.quantile(s, [0.25, 0.75]), color=col, lw=2.2, alpha=0.35,
                        solid_capstyle="butt")
            ax.plot(xs, med, color=col, lw=1.3, label=g, **mk)
            if uf is not None:
                s = np.array([1 - uf.loc[a, b] for a, b in itertools.combinations(members, 2)])
                rows.append(dict(panel=title, level=uflab.replace("\n", " "), group=g,
                                 median_similarity=np.median(s), q25=np.quantile(s, 0.25),
                                 q75=np.quantile(s, 0.75)))
                ax.plot([len(xs) + 0.4] * 2, np.quantile(s, [0.25, 0.75]), color=col, lw=2.2, alpha=0.35,
                        solid_capstyle="butt")
                ax.plot([len(xs) + 0.4], [np.median(s)], ls="", color=col, **mk)
        ticks, labs = list(xs), xlab[:]
        if uf is not None:
            ax.axvline(len(xs) - 0.3, color="#bbbbbb", lw=0.6, ls=(0, (3, 2)))
            ticks.append(len(xs) + 0.4)
            labs.append(uflab)
        ax.set_xticks(ticks, labs, rotation=90, fontsize=6)
        ax.set_xlim(-0.5, len(xs) + 0.9)
        ax.set_ylim(0, 1)
        ax.set_title(title, fontsize=7, fontweight="bold")
        ax.tick_params(axis="x", length=0, pad=2)
        ax.tick_params(axis="y", labelsize=6)
        ax.spines[["top", "right"]].set_visible(False)
        if j == 0:
            ax.set_ylabel("Similarity between two samples\nof the same group (1 - d)", fontsize=6.5)
            ax.text(-0.5, 1.08, "A", transform=ax.transAxes, fontsize=10, fontweight="bold")

    ax = fig.add_subplot(gs[0, 4])
    fx = {"braycurtis": 0, "hellinger": 1, "jaccard": 2}
    for k, (fname, off) in enumerate((("function, coverage", 0.0), ("function, count", 3.6))):
        for meas, x0 in fx.items():
            for g, members in groups.items():
                col = C.GROUP_COLOR[g]
                x = off + x0 + (-0.17 if g == "Source" else 0.17)
                s = 1 - _pair_d(lt[fname], members, meas)
                rows.append(dict(panel="Function", level=f"{fname}, {meas}", group=g,
                                 median_similarity=np.median(s), q25=np.quantile(s, 0.25),
                                 q75=np.quantile(s, 0.75)))
                ax.plot([x, x], np.quantile(s, [0.25, 0.75]), color=col, lw=2.2, alpha=0.35,
                        solid_capstyle="butt")
                ax.plot([x], [np.median(s)], ls="", marker=F2_MARKER[g], ms=3.4, color=col,
                        mec="#333333" if g == "Enrichment" else "white", mew=0.4)
    ax.axvline(3.0 - 0.2, color="#bbbbbb", lw=0.6, ls=(0, (3, 2)))
    ax.set_xticks([0, 1, 2, 3.6, 4.6, 5.6], ["Bray-Curtis", "Hellinger", "Jaccard"] * 2, rotation=90, fontsize=6)
    ax.text(1, 0.6, "coverage-\nbased", ha="center", va="top", fontsize=6, color="#555555")
    ax.text(4.6, 0.6, "count-\nbased", ha="center", va="top", fontsize=6, color="#555555")
    hs = [plt.Line2D([], [], color=C.GROUP_COLOR[g], lw=1.3, marker=F2_MARKER[g], ms=3.4,
                     mec="#333333" if g == "Enrichment" else "white", mew=0.4, label=g) for g in groups]
    ax.legend(handles=hs, loc="lower center", fontsize=6, frameon=False, handlelength=1.5)
    ax.set_xlim(-0.6, 6.2)
    ax.set_ylim(0, 1)
    ax.set_title("Function", fontsize=7, fontweight="bold")
    ax.tick_params(axis="x", length=0, pad=2)
    ax.tick_params(axis="y", labelsize=6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(-0.3, 1.08, "B", transform=ax.transAxes, fontsize=10, fontweight="bold")
    nm_rows = []
    (C.OUT / "si").mkdir(exist_ok=True)
    C.polish(fig)
    fig.savefig(C.OUT / "si" / "S_convergence_robustness.pdf")
    fig.savefig(C.OUT / "si" / "S_convergence_robustness.png", dpi=220)
    plt.close(fig)
    t = pd.DataFrame(rows)
    t.to_csv(C.OUT / "si" / "S_convergence_robustness.csv", index=False)
    return t


def ordination_robustness_table() -> str:
    """The PERMANOVA and PERMDISP results of ordination_robustness() as a LaTeX table,
    which replaces the 24-panel figure in the SI."""
    t = pd.read_csv(C.OUT / "si" / "S_ordination_robustness.csv")
    f = lambda p: "$<$0.001" if p < 0.001 else f"{p:.3f}"
    lines, last = [], None
    for _, r in t.iterrows():
        if r.dataset != last:
            if last is not None:
                lines.append(r"\addlinespace[2pt]")
            lines.append(rf"\multicolumn{{5}}{{@{{}}l}}{{\textit{{{r.dataset}}}}} \\")
            last = r.dataset
        sp = [f"{g} {r[f'spread_{g}']:.2f}" for g in ("Source", "Enrichment", "V", "NV", "NVNN")
              if pd.notna(r.get(f"spread_{g}"))]
        name = r.distance.replace("Bray-Curtis", "Bray--Curtis")
        func = r.dataset.endswith("function")
        if func and not r.distance.startswith("count"):
            name = "coverage-based, " + name
        main = {"Media: function": "count-based, Bray-Curtis"}.get(r.dataset, "Bray-Curtis, relative")
        if r.distance == main:
            name += " (main)"
        lines.append(f"\\quad {name} & {r.PCo1:.0f}/{r.PCo2:.0f} & {f(r.p)} & {f(r.disp_p)} & "
                     f"{', '.join(sp)} \\\\")
    tex = ("\\begin{tabular}{@{}lrrrl@{}}\n\\toprule\n"
           "Measure & Axes 1/2 (\\%) & PERMANOVA $p$ & PERMDISP $p$ & Spread (distance to centroid) \\\\\n"
           "\\midrule\n" + "\n".join(lines) + "\n\\bottomrule\n\\end{tabular}\n")
    (C.OUT / "si" / "S_ordination_robustness_table.tex").write_text(tex)
    return tex


def main() -> None:
    run_tables()
    run_symclatron()
    run_si_figures()
    print(ordination_robustness().round(3).to_string())
    print(pcoa_by_rank().round(3).to_string())
    print(convergence_robustness().round(3).to_string())
    ordination_robustness_table()
    print(morphotypes().round(3).to_string())
    print(media_category_partition_figure().round(3).to_string())
    print(*media_pairwise_dissimilarity(), sep='\n')


if __name__ == "__main__":
    main()
