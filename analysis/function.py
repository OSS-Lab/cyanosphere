"""
The two-metric functional analysis of the enrichments against the reservoir:
classes, robustness, Table 1, the MAG-level check of the count-based metric and
its calibration (Fig 3, Table 1, Supplementary Note S1, Table S6).

Run:  python function.py      (from analysis/, after the modules before it in
      the README order)
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import common as C
from scipy.stats import spearmanr, theilslopes

# ============================================================================
# Two-metric classification, enrichments against reservoir: run_two_metric()
# ----------------------------------------------------------------------------
# Produces
#   outputs/table1_two_metric.{csv,tex}      Table 1
#   outputs/module_niche_classification.csv  every module, both metrics, its class
#   outputs/denominator_robustness.csv       the classification with other ways of
#                                            counting populations (Table S6)
#   outputs/fig3.{pdf,svg}                   Figure 3, all five panels
#
# Run:  python function.py
# ============================================================================

# Modules drawn in panel E: every row of Table 1, plus the principal
# photosynthesis, nitrogen and cofactor modules, so the panel and the table
# describe the same functions.
PANEL_E_MODULES = [
    "M00161", "M00163", "M00145", "M00162", "M00165", "M00097", "M00112",
    "M00932", "M00615", "M00531", "M00529", "M00530", "M00899", "M00122",
    "M00924", "M00123", "M00118", "M00970",
    "M00597", "M00156", "M00985", "M00555", "M00134", "M00642",
    "M00116", "M00133",
]


SHORT = {
    "M00161": "Photosystem II", "M00163": "Photosystem I",
    "M00145": "NAD(P)H:quinone oxidoreductase", "M00162": "Cytochrome b6f complex",
    "M00165": "Calvin cycle", "M00097": "beta-Carotene biosynthesis",
    "M00112": "Tocopherol biosynthesis", "M00932": "Phylloquinone biosynthesis",
    "M00615": "Nitrate assimilation", "M00531": "Assimilatory nitrate reduction",
    "M00529": "Denitrification", "M00530": "Dissimilatory nitrate reduction",
    "M00899": "Thiamine salvage", "M00122": "Cobalamin biosynthesis (aerobic)",
    "M00924": "Cobalamin biosynthesis (anaerobic)", "M00123": "Biotin biosynthesis",
    "M00118": "Glutathione biosynthesis", "M00970": "Proline degradation",
    "M00597": "Anoxygenic photosystem II", "M00156": "Cytochrome c oxidase (cbb3)",
    "M00985": "Sulfide oxidation", "M00555": "Betaine biosynthesis",
    "M00134": "Polyamine biosynthesis (ornithine)",
    "M00642": "Multidrug efflux pump MexJK-OprM",
    "M00116": "Menaquinone biosynthesis", "M00133": "Polyamine biosynthesis (arginine)",
}


ENR_ORDER = ["M_63", "M_67", "M_71", "M_77", "M_132", "M_133", "M_134", "M_135"]


# ---------------------------------------------------------------- the analysis
def two_metric(marker: str = "rpl19"):
    cov_m, cnt_m, names = C.metric_matrices("enrichments", marker)
    cov = C.differential(cov_m, C.ENRICHMENT, C.SOURCE)
    cnt = C.differential(cnt_m, C.ENRICHMENT, C.SOURCE)
    cls = C.niche_class(cov, cnt)
    t = pd.DataFrame({
        "cov_l2fc": cov.l2fc, "cov_p": cov.p, "cov_padj": cov.padj,
        "cnt_l2fc": cnt.l2fc, "cnt_p": cnt.p, "cnt_padj": cnt.padj,
        "niche_class": cls,
    }).join(names)
    # between-replicate variability of each module across the eight enrichments
    e = cov_m[[s for s in C.ENRICHMENT if s in cov_m.columns]]
    t["cov_cv"] = (e.std(axis=1) / e.mean(axis=1)).replace([np.inf, -np.inf], np.nan)
    return t, cov_m, cnt_m


def denominator_robustness(reference: pd.DataFrame) -> pd.DataFrame:
    """
    Does the classification depend on the analytical choices behind it? Each row
    changes one choice and recomputes every class:

    denominator   RpL19 clusters (primary), RpS16 clusters, detected MAGs, or
                  the anvi'o assembly-level L19 estimate used until Sep 2026
    depth         coverage divided by the sample's summed RpL19 cluster coverage,
                  i.e. module coverage per genome equivalent, which removes
                  sequencing depth (a single-copy marker, so no compositional
                  closure over modules)
    test          the SciPy asymptotic rank-sum in place of the exact test
    """
    rows = []

    def add(label, t):
        agree = (t.niche_class == reference.niche_class.reindex(t.index)).mean()
        rows.append(dict(
            variant=label,
            n_both=(t.niche_class == "both metrics").sum(),
            n_count_only=(t.niche_class == "count-based only").sum(),
            n_coverage_only=(t.niche_class == "coverage-based only").sum(),
            n_neither=(t.niche_class == "neither").sum(),
            agreement_with_primary=round(100 * agree, 1),
            count_only_modules=";".join(sorted(t.index[t.niche_class == "count-based only"])),
        ))

    for marker in ("rpl19", "rps16", "mags", "anvio_l19"):
        add(f"denominator: {marker}", two_metric(marker)[0])

    cov_m, cnt_m, _ = C.metric_matrices("enrichments", "rpl19")
    cl = pd.read_csv(C.DATA / "rpl19_cluster_coverage.csv", index_col=0)
    ge = cl.sum().reindex(cov_m.columns)
    cov = C.differential(cov_m.div(ge, axis=1), C.ENRICHMENT, C.SOURCE)
    cnt = C.differential(cnt_m, C.ENRICHMENT, C.SOURCE)
    add("coverage per RpL19 genome equivalent", C.niche_class(cov, cnt).to_frame())

    cov = C.differential(cov_m, C.ENRICHMENT, C.SOURCE)
    for d in (cov, cnt):
        d["padj"] = d["padj_asymptotic"]
    add("asymptotic rank-sum test", C.niche_class(cov, cnt).to_frame())
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- panel E data
def division_of_labour(detect: float = C.MAG_DETECT) -> pd.DataFrame:
    """
    Per enrichment replicate, the share of the coverage of the MAGs encoding a
    module that belongs to Cyanobacteriota.

    A MAG encodes a module at pathwise completeness >= ENCODING and is present in
    a replicate at mean coverage >= `detect`. Coverage comes from
    `mag_sample_coverage.tsv`, which is linked to each MAG by sequence (contig
    headers, see pipeline/05_link_mags_to_bins.py), and phylum from `mag_quality.tsv`, so no
    name matching is involved. A cell is NA when no detected MAG encodes the
    module in that replicate: absence of an encoder is not 0% cyanobacterial.
    """
    mc = pd.read_csv(C.DATA / "mag_module_pathwise_completeness.tsv", sep="\t").set_index("module")
    cov = pd.read_csv(C.DATA / "mag_sample_coverage.tsv", sep="\t").set_index("MAG_name")
    phylum = pd.read_csv(C.DATA / "mag_quality.tsv", sep="\t").set_index("MAG_name").phylum
    mags = [g for g in mc.columns if g in cov.index]
    cyano = {g for g in mags if phylum.get(g) == C.CYANO_PHYLUM}
    rows = []
    for m in PANEL_E_MODULES:
        if m not in mc.index:
            continue
        enc = [g for g in mags if mc.loc[m, g] >= C.ENCODING]
        for smp in ENR_ORDER:
            present = [g for g in enc if cov.loc[g, smp] >= detect]
            c = float(cov.loc[[g for g in present if g in cyano], smp].sum())
            h = float(cov.loc[[g for g in present if g not in cyano], smp].sum())
            rows.append(dict(module=m, sample=smp,
                             n_cyano=sum(g in cyano for g in present),
                             n_other=sum(g not in cyano for g in present),
                             cyano_percent=100 * c / (c + h) if present else np.nan))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- the figure
def panel_B(ax, cov_m):
    samples = [s for s in C.ENRICHMENT + C.SOURCE if s in cov_m.columns]
    comp = C.compositional(cov_m)
    X = StandardScaler().fit_transform(comp[samples].T.values)
    p = PCA(n_components=2).fit(X)
    sc = p.transform(X)
    var = p.explained_variance_ratio_ * 100
    for grp, members in [("Source", C.SOURCE), ("Enrichment", C.ENRICHMENT)]:
        idx = [i for i, s in enumerate(samples) if s in members]
        ax.scatter(sc[idx, 0], sc[idx, 1], s=64, color=C.GROUP_COLOR[grp],
                   edgecolor="white", lw=0.8, label=grp, zorder=3)
    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    F, p = C.permanova(comp, {"Enrichment": C.ENRICHMENT, "Source": C.SOURCE})
    ax.set_title("B   Functional community structure", loc="left")
    ax.annotate(f"PERMANOVA $F$ = {F:.1f}, $p$ = {p:.3f}", xy=(0.03, 0.03),
                xycoords="axes fraction", fontsize=7.5, color="#555555")
    ax.legend(loc="best")
    return var, F, p


def panel_D(ax, t):
    x = t.cov_l2fc.values
    y = -np.log10(t.cov_padj.clip(1e-12).values)
    col = np.where(t.niche_class == "count-based only", C.CLASS_COLOR["count-based only"],
                   np.where(t.niche_class == "both metrics", C.CLASS_COLOR["both metrics"],
                            "#c8ccd0"))
    ax.scatter(x, y, s=18, c=col, edgecolor="none", alpha=0.85, zorder=3)
    ax.axhline(-np.log10(C.FDR), color="#b0392b", lw=0.8, ls="--")
    ax.axvline(C.MIN_L2FC, color="#999999", lw=0.7, ls=":")
    ax.axvline(-C.MIN_L2FC, color="#999999", lw=0.7, ls=":")
    offsets = {"M00163": (4, 6), "M00161": (4, -9), "M00530": (5, 2),
               "M00165": (-6, 7), "M00145": (-4, -10), "M00597": (5, 2),
               "M00156": (5, 2), "M00985": (5, -6)}
    for m, off in offsets.items():
        if m in t.index:
            ax.annotate(SHORT.get(m, m), (t.loc[m, "cov_l2fc"],
                        -np.log10(max(t.loc[m, "cov_padj"], 1e-12))),
                        fontsize=6.3, xytext=off, textcoords="offset points",
                        color="#333333")
    ax.set_xlabel("log$_2$ fold change, enrichment / source (coverage-based)")
    ax.set_ylabel("$-$log$_{10}$ FDR")
    ax.set_title("D   Differentially abundant functions", loc="left")
    ax.legend(handles=[mpatches.Patch(color=C.CLASS_COLOR["both metrics"], label="both metrics"),
                       mpatches.Patch(color=C.CLASS_COLOR["count-based only"], label="count-based only"),
                       mpatches.Patch(color="#c8ccd0", label="not enriched")],
              loc="upper left", fontsize=7)


def panel_E(ax, dol, t):
    """Split-triangle cells, as in the published panel: lower-left heterotroph,
    upper-right cyanobacterial, shaded by share, with the percentage printed.
    Rows are grouped by two-metric class and, within a class, ordered from
    entirely cyanobacterial to entirely heterotrophic."""
    piv = dol.pivot_table(index="module", columns="sample", values="cyano_percent")
    cls = t.niche_class.reindex(piv.index)
    order = []
    for block in ("both metrics", "count-based only", "coverage-based only"):
        members = [m for m in piv.index if cls.get(m) == block]
        members.sort(key=lambda m: -np.nanmean(piv.loc[m].values))
        order += [(block, m) for m in members]
    n = len(order)
    h = 0.5
    spans = {}
    for i, (block, m) in enumerate(order):
        yy = n - 1 - i
        spans.setdefault(block, []).append(yy)
        for j, s in enumerate(ENR_ORDER):
            v = piv.loc[m, s] if s in piv.columns else np.nan
            if np.isnan(v):
                ax.add_patch(mpatches.Rectangle((j - h, yy - h), 1, 1,
                             facecolor="#f0f0f0", edgecolor="white", lw=0.5))
                continue
            ax.add_patch(mpatches.Polygon([[j - h, yy - h], [j + h, yy - h], [j - h, yy + h]],
                         closed=True, facecolor=C.HETERO_C, edgecolor="white", lw=0.5,
                         alpha=0.12 + 0.88 * (100 - v) / 100))
            ax.add_patch(mpatches.Polygon([[j + h, yy - h], [j + h, yy + h], [j - h, yy + h]],
                         closed=True, facecolor=C.CYANO_C, edgecolor="white", lw=0.5,
                         alpha=0.12 + 0.88 * v / 100))
            ax.text(j, yy, f"{v:.0f}", ha="center", va="center", fontsize=5.6, zorder=5)
        ax.text(-0.72, yy, f"{m}  {SHORT.get(m, m)}", ha="right", va="center", fontsize=7.2)
    for block, ys in spans.items():
        ax.text(len(ENR_ORDER) - 0.30, np.mean(ys), block.replace(" ", "\n"),
                rotation=270, ha="left", va="center", fontsize=7,
                fontweight="bold", color=C.CLASS_COLOR[block])
        ax.axhline(max(ys) + 0.5, color="#888888", lw=0.8)
    ax.set_xlim(-8.6, len(ENR_ORDER) + 1.2)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xticks(range(len(ENR_ORDER)))
    ax.set_xticklabels(ENR_ORDER, fontsize=7.5, rotation=45, ha="left")
    ax.xaxis.tick_top()
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_aspect("equal")
    ax.set_title("E   Who encodes each enriched function (cyanobacterial % shown)",
                 loc="left", y=1.05)
    ax.legend(handles=[mpatches.Patch(facecolor=C.CYANO_C, label="Cyanobacteria"),
                       mpatches.Patch(facecolor=C.HETERO_C, label="Heterotrophs"),
                       mpatches.Patch(facecolor="#f0f0f0", label="no encoding MAG detected")],
              loc="upper center", bbox_to_anchor=(0.45, -0.02), ncol=3, fontsize=7.5)


def placeholder(ax, letter, title):
    ax.add_patch(mpatches.FancyBboxPatch((0.03, 0.03), 0.94, 0.94,
                 transform=ax.transAxes, boxstyle="round,pad=0.01",
                 facecolor="#f6f6f6", edgecolor="#c2c2c2", lw=1.0, ls=(0, (5, 4))))
    ax.text(0.5, 0.55, title, transform=ax.transAxes, ha="center", va="center",
            fontsize=9, color="#8a8a8a")
    ax.text(0.5, 0.40, "schematic, composited in Affinity", transform=ax.transAxes,
            ha="center", va="center", fontsize=7, color="#b0b0b0", style="italic")
    ax.set_title(letter, loc="left")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


# ---------------------------------------------------------------- table 1
def write_table1(t: pd.DataFrame) -> None:
    sig_cov = t[(t.niche_class.isin(["both metrics", "coverage-based only"]))]
    top_cov = sig_cov.sort_values("cov_l2fc", ascending=False).head(10)
    only_cnt = t[t.niche_class == "count-based only"].sort_values("cnt_l2fc", ascending=False).head(10)
    rows_df = pd.concat([top_cov, only_cnt]).sort_values("cov_l2fc", ascending=False)
    lines = []
    for m, r in rows_df.iterrows():
        def cell(v, p, up):
            s = f"{v:.2f}"
            ps = "$<$0.001" if p < 0.001 else f"{p:.3f}"
            return (rf"\textbf{{{s}}}", rf"\textbf{{{ps}}}") if up else (s, ps)
        c1, c2 = cell(r.cov_l2fc, r.cov_padj, r.cov_padj < C.FDR and r.cov_l2fc > C.MIN_L2FC)
        c3, c4 = cell(r.cnt_l2fc, r.cnt_padj, r.cnt_padj < C.FDR and r.cnt_l2fc > C.MIN_L2FC)
        name = C.tex_escape(r.module_name)
        if len(name) > 52:
            name = name[:50] + r"\ldots"
        lines.append(f"{m} & {name} & {c1} & {c2} & {c3} & {c4} \\\\")
    tex = r"""\begin{table}[h]
\centering
\caption{KEGG modules most strongly enriched in the end-point enrichment
communities relative to the reservoir communities. Fold changes
(log\textsubscript{2}FC, group means, per-module half-minimum pseudocount) and
Benjamini--Hochberg-adjusted $p$ values (two-sided Wilcoxon rank-sum, $n=8$ per
group) are given separately for the coverage-based and the count-based metric.
Bold marks a significant increase on that metric ($p_\mathrm{adj}<0.05$ and
log\textsubscript{2}FC $>0.5$). Rows are the ten modules with the highest
coverage-based fold change among those enriched on the coverage-based metric,
together with the modules enriched on the count-based metric alone; all are
ordered by coverage-based fold change. The two metrics are on different scales,
so fold changes are comparable down a column but not across the two
(see \methodsref).}
\label{tab:kegg_top10}
\footnotesize
\begin{tabular}{@{}llrrrr@{}}
\toprule
 &  & \multicolumn{2}{c}{\textbf{Coverage-based}} & \multicolumn{2}{c}{\textbf{Count-based}} \\
\cmidrule(l){3-4}\cmidrule(l){5-6}
\textbf{Module} & \textbf{Module name} & \textbf{log\textsubscript{2}FC} & \textbf{$p_\mathrm{adj}$} & \textbf{log\textsubscript{2}FC} & \textbf{$p_\mathrm{adj}$} \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    (C.OUT / "table1_two_metric.tex").write_text(tex)
    print(f"wrote outputs/table1_two_metric.tex ({len(lines)} rows)")


def run_two_metric() -> None:
    C.style()
    t, cov_m, cnt_m = two_metric("rpl19")
    t.to_csv(C.OUT / "module_niche_classification.csv")
    t.to_csv(C.OUT / "table1_two_metric.csv")
    write_table1(t)

    print("\nniche classes (RpL19 denominator):")
    print(t.niche_class.value_counts().to_string())
    print("\nmodules enriched on the count-based metric alone:")
    print(t[t.niche_class == "count-based only"][["module_name", "cov_l2fc", "cnt_l2fc"]]
          .to_string(float_format=lambda v: f"{v:6.2f}"))
    print("\nmedian coefficient of variation across the eight enrichments, by class:")
    print(t.groupby("niche_class").cov_cv.median().round(3).to_string())

    rob = denominator_robustness(t)
    rob.to_csv(C.OUT / "denominator_robustness.csv", index=False)
    print("\ndenominator robustness:")
    print(rob.drop(columns="count_only_modules").to_string(index=False))

    dol = division_of_labour()
    dol.to_csv(C.OUT / "fig3e_division_of_labour.csv", index=False)

    fig = plt.figure(figsize=(11.5, 15.5))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.40, 0.62, 1.15], hspace=0.30, wspace=0.24)
    placeholder(fig.add_subplot(gs[0, 0]), "A", "Functional annotation workflow\n(reads, contigs, KEGG modules)")
    placeholder(fig.add_subplot(gs[0, 1]), "C", "The imposed environment\n(carbon-free, nitrate-replete, 12:12 h light:dark)")
    var, F, p = panel_B(fig.add_subplot(gs[1, 0]), cov_m)
    print(f"\npanel B: PC1={var[0]:.1f}% PC2={var[1]:.1f}%, PERMANOVA F={F:.2f} p={p:.3f}")
    panel_D(fig.add_subplot(gs[1, 1]), t)
    panel_E(fig.add_subplot(gs[2, :]), dol, t)
    C.save(fig, "fig3")


# ============================================================================
# MAG-based check of the count-based metric: run_mag_control()
# ----------------------------------------------------------------------------
# The count-based metric is computed on raw contigs (stepwise copy number over the
# number of coexisting populations) and never looks at a genome. This compares it
# with the MAGs inside the eight enrichments. For every module and replicate: how
# many of the dereplicated MAGs detected in that replicate encode the module
# (pathwise completeness >= 0.75, the Fig 3E rule), divided by the number of MAGs
# detected there. If the count-based metric tracks how widely a function is held,
# this MAG-based breadth should follow it across modules.
#
# The check cannot reproduce the source-versus-enrichment fold change on MAGs:
# MAG coverage in data/ covers the eight enrichments only, so there is no MAG-based
# value for the reservoir. Functions carried by populations that were not binned
# show breadth on the count-based metric that the MAGs miss.
#
# Reads only data/. Writes outputs/mag_control_*.
# ============================================================================

DESC_COLS = ["module_name", "module_class", "module_category", "module_subcategory"]


def mag_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (module x MAG pathwise completeness, MAG x enrichment detection),
    from the sequence-linked MAG table shared with Figure 3E."""
    enc, cov, _ = C.mag_tables()
    detected = cov[C.ENRICHMENT] >= C.MAG_DETECT
    return enc, detected


def mag_breadth(enc: pd.DataFrame, detected: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per module and replicate: MAGs encoding it, and that as a fraction of MAGs present."""
    carries = (enc >= C.ENCODING)          # module x MAG
    n_enc = carries.astype(int) @ detected.astype(int)   # module x sample
    n_present = detected.sum(axis=0)                     # sample
    return n_enc, n_enc.div(n_present, axis=1)


def run_mag_control() -> None:
    enc, detected = mag_tables()
    n_enc, breadth = mag_breadth(enc, detected)
    n_present = detected.sum(axis=0)
    print(f"{enc.shape[1]} MAGs mapped, {enc.shape[0]} modules")
    print("MAGs detected per replicate:")
    print(n_present.to_string())

    cov, cnt, names = C.metric_matrices("enrichments")
    shared = cnt.index.intersection(breadth.index)
    cnt_e = cnt.loc[shared, C.ENRICHMENT]
    cov_e = cov.loc[shared, C.ENRICHMENT]
    br_e = breadth.loc[shared, C.ENRICHMENT]
    ne_e = n_enc.loc[shared, C.ENRICHMENT]
    print(f"\n{len(shared)} modules present in both the contig and MAG tables")

    # --- 1. per-replicate concordance across modules ---------------------------
    rows = []
    for s in C.ENRICHMENT:
        r_cnt = spearmanr(cnt_e[s], br_e[s])
        r_cov = spearmanr(cov_e[s], br_e[s])
        rows.append(
            {
                "sample": s,
                "n_MAGs_detected": int(n_present[s]),
                "rho_count_vs_MAG_breadth": round(r_cnt.statistic, 3),
                "p_count": r_cnt.pvalue,
                "rho_coverage_vs_MAG_breadth": round(r_cov.statistic, 3),
                "p_coverage": r_cov.pvalue,
            }
        )
    conc = pd.DataFrame(rows)
    conc.to_csv(C.OUT / "mag_control_concordance.csv", index=False)
    print("\n--- concordance per replicate (Spearman across modules) ---")
    print(conc.to_string(index=False))
    print(
        f"\nmean rho: count-based {conc.rho_count_vs_MAG_breadth.mean():.3f}, "
        f"coverage-based {conc.rho_coverage_vs_MAG_breadth.mean():.3f}"
    )

    # --- 2. by two-metric niche class -----------------------------------------
    cov_t = C.differential(cov, C.ENRICHMENT, C.SOURCE)
    cnt_t = C.differential(cnt, C.ENRICHMENT, C.SOURCE)
    cls = C.niche_class(cov_t, cnt_t).reindex(shared)

    per_module = pd.DataFrame(
        {
            "module": shared,
            "module_name": names.reindex(shared)["module_name"].values,
            "niche_class": cls.values,
            "l2fc_coverage": cov_t.reindex(shared)["l2fc"].values,
            "l2fc_count": cnt_t.reindex(shared)["l2fc"].values,
            "mean_count_metric": cnt_e.mean(axis=1).values,
            "mean_MAGs_encoding": ne_e.mean(axis=1).round(2).values,
            "mean_MAG_breadth": br_e.mean(axis=1).round(3).values,
            "min_MAG_breadth": br_e.min(axis=1).round(3).values,
            "max_MAG_breadth": br_e.max(axis=1).round(3).values,
        }
    ).set_index("module")
    per_module.to_csv(C.OUT / "mag_control_per_module.csv")

    summ = (
        per_module.groupby("niche_class")
        .agg(
            n_modules=("mean_MAG_breadth", "size"),
            median_MAGs_encoding=("mean_MAGs_encoding", "median"),
            median_MAG_breadth=("mean_MAG_breadth", "median"),
        )
        .round(3)
    )
    summ.to_csv(C.OUT / "mag_control_by_class.csv")
    print("\n--- MAG breadth by two-metric class ---")
    print(summ.to_string())

    # the six community-derived modules, individually
    six = per_module[per_module.niche_class == "count-based only"]
    print("\n--- the count-based-only (community-derived) modules on MAGs ---")
    print(
        six[["module_name", "l2fc_count", "mean_MAGs_encoding", "mean_MAG_breadth"]]
        .to_string()
    )

    # --- 3. where do the two disagree? ----------------------------------------
    # rank each module by the contig count metric and by MAG breadth; large
    # positive residuals are functions the contigs see as widely held but the
    # MAGs do not, which is what unbinned populations should produce
    rank_cnt = per_module.mean_count_metric.rank(pct=True)
    rank_mag = per_module.mean_MAG_breadth.rank(pct=True)
    per_module["rank_gap_contig_minus_MAG"] = (rank_cnt - rank_mag).round(3)
    per_module.sort_values("rank_gap_contig_minus_MAG", ascending=False).to_csv(
        C.OUT / "mag_control_per_module.csv"
    )
    top = per_module.nlargest(10, "rank_gap_contig_minus_MAG")
    print("\n--- widest contig-over-MAG rank gaps (functions the MAGs miss) ---")
    print(top[["module_name", "niche_class", "rank_gap_contig_minus_MAG"]].to_string())


# ============================================================================
# Table 1 in blocks: run_table1()
# ----------------------------------------------------------------------------
# One block per class. Each block header gives the class name, the class size and
# the rule used to select the rows shown. Unadjusted p and BH-adjusted q are printed
# for both metrics; bold marks the call, made on q (n = 8 against 8). The MAGs
# column is the mean number of dereplicated MAGs detected in an enrichment
# replicate that encode the module (from run_mag_control()).
#
# Writes outputs/table1_blocked.tex and a plain-text preview to stdout.
# Needs the outputs of run_mag_control().
# ============================================================================

TOP_N = 10


# blocks in the order they are argued in the text: what the medium imposes, what
# the community generates, then the class that is abundant but not widely held
BLOCKS = [
    ("both metrics", "Imposed (significant on both metrics)",
     f"top {TOP_N} by coverage-based log$_2$FC"),
    ("count-based only", "Community-derived (count-based metric only)",
     "all shown"),
    ("coverage-based only", "Coverage-based metric only",
     f"top {TOP_N} by coverage-based log$_2$FC"),
]


SHORTEN = {
    "Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate": "Glycolysis (Embden--Meyerhof)",
    "Dissimilatory nitrate reduction, nitrate => ammonia": "Dissimilatory nitrate reduction",
    "Polyamine biosynthesis, arginine => ornithine => putrescine": "Polyamine biosynthesis",
    "Betaine biosynthesis, choline => betaine": "Betaine biosynthesis",
    "Sulfide oxidation, sulfide => sulfur": "Sulfide oxidation",
    "Multidrug resistance, efflux pump MexJK-OprM": "Multidrug efflux pump MexJK-OprM",
    "Cytochrome c oxidase, cbb3-type": "Cytochrome $c$ oxidase, cbb$_3$-type",
}


def fmt_p(v: float) -> str:
    if pd.isna(v):
        return "--"
    return "$<$0.001" if v < 0.001 else f"{v:.3f}"


def build() -> pd.DataFrame:
    t = pd.read_csv(C.OUT / "table1_two_metric.csv").set_index("module")
    try:
        mag = pd.read_csv(C.OUT / "mag_control_per_module.csv").set_index("module")
        t["mags"] = mag["mean_MAGs_encoding"]
    except FileNotFoundError:
        raise SystemExit("run function.py (run_mag_control) first")
    t["name"] = t["module_name"].map(lambda s: SHORTEN.get(s, s))
    return t


def rows_for(t: pd.DataFrame, cls: str) -> pd.DataFrame:
    sub = t[t.niche_class == cls].sort_values("cov_l2fc", ascending=False)
    return sub if cls == "count-based only" else sub.head(TOP_N)


def run_table1() -> None:
    t = build()

    lines = [
        r"\begin{table}[h]",
        r"\centering",
        r"\caption{KEGG modules enriched in the end-point enrichment communities "
        r"relative to the reservoir communities, grouped by the two-metric "
        r"classification. Each block header gives the full size of the class and "
        r"the rule used to select the rows shown, so that a truncated block and a "
        r"complete one are not confused. Fold changes are group means with a "
        r"per-module half-minimum pseudocount; $p$ is a two-sided Wilcoxon "
        r"rank-sum ($n=8$ per group) and $q$ its Benjamini--Hochberg adjustment "
        r"across all 235 modules tested. Bold marks a significant increase on that "
        r"metric ($q<0.05$ and log\textsubscript{2}FC $>0.5$). MAGs is the mean "
        r"number of dereplicated MAGs detected in an enrichment replicate that "
        r"encode the module (pathwise completeness $\geq0.75$), an independent "
        r"genome-level check on the count-based metric. The two metrics are on "
        r"different scales, so fold changes are comparable down a column but not "
        r"across the two (see \methodsref).}",
        r"\label{tab:kegg_blocked}",
        r"\footnotesize",
        r"\begin{tabular}{@{}llrrrrrrr@{}}",
        r"\toprule",
        r" &  & \multicolumn{3}{c}{\textbf{Coverage-based}} & "
        r"\multicolumn{3}{c}{\textbf{Count-based}} & \\",
        r"\cmidrule(lr){3-5}\cmidrule(lr){6-8}",
        r"Module & Name & log\textsubscript{2}FC & $p$ & $q$ & "
        r"log\textsubscript{2}FC & $p$ & $q$ & MAGs \\",
        r"\midrule",
    ]

    preview = []
    for i, (cls, label, rule) in enumerate(BLOCKS):
        n_total = int((t.niche_class == cls).sum())
        sub = rows_for(t, cls)
        if i:
            lines.append(r"\addlinespace[2pt]\midrule\addlinespace[2pt]")
        lines.append(
            rf"\multicolumn{{9}}{{@{{}}l}}{{\textit{{{label}}}\quad"
            rf"$n={n_total}$\quad\small({rule})}} \\"
        )
        lines.append(r"\addlinespace[1pt]")
        preview.append(f"\n{label}   n={n_total}   ({rule})")
        preview.append(f"  {'module':8}{'name':38}{'covFC':>7}{'p':>9}{'q':>8}"
                       f"{'cntFC':>7}{'p':>9}{'q':>8}{'MAGs':>7}")

        for m, r in sub.iterrows():
            cov_b = (r.cov_padj < C.FDR) and (r.cov_l2fc > C.MIN_L2FC)
            cnt_b = (r.cnt_padj < C.FDR) and (r.cnt_l2fc > C.MIN_L2FC)
            cov_fc = rf"\textbf{{{r.cov_l2fc:.2f}}}" if cov_b else f"{r.cov_l2fc:.2f}"
            cnt_fc = rf"\textbf{{{r.cnt_l2fc:.2f}}}" if cnt_b else f"{r.cnt_l2fc:.2f}"
            lines.append(
                f"{m} & {C.tex_escape(r['name'])} & {cov_fc} & {fmt_p(r.cov_p)} & "
                f"{fmt_p(r.cov_padj)} & {cnt_fc} & {fmt_p(r.cnt_p)} & "
                f"{fmt_p(r.cnt_padj)} & {r.mags:.1f} \\\\"
            )
            nm = r["name"][:37]
            preview.append(
                f"  {m:8}{nm:38}{r.cov_l2fc:7.2f}{r.cov_p:9.4f}{r.cov_padj:8.4f}"
                f"{r.cnt_l2fc:7.2f}{r.cnt_p:9.4f}{r.cnt_padj:8.4f}{r.mags:7.1f}"
            )

    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    (C.OUT / "table1_blocked.tex").write_text("\n".join(lines) + "\n")
    print("\n".join(preview))
    print("\nwrote outputs/table1_blocked.tex")


# ============================================================================
# Calibration of the count-based metric in genomes: run_calibration()
# ----------------------------------------------------------------------------
# If a module goes from 0.2 to 0.4 on the count-based metric, how many more genomes
# carry it? For every module in every enrichment replicate:
#
#   count-based metric   copy number / RpL19 populations, from raw contigs
#   MAG encoders         number of detected MAGs encoding the module at
#                        pathwise completeness >= 0.75
#   MAG breadth          MAG encoders / MAGs detected in the replicate
#
# Two calibrations: against MAG breadth (a fraction against a fraction, since the
# denominators differ: 45 to 94 populations, 12 to 28 MAGs), and a lookup from
# metric value to the expected number of encoding MAGs, with its spread. Only 12
# to 28 MAGs are detected per replicate at mean coverage >= 1, so the integer
# axis is coarse: this calibrates direction and rough magnitude, not an exact
# genome count.
#
# Needs the outputs of run_mag_control(). Writes outputs/calibration_*.
# ============================================================================

# the metric values a reader will meet in the text and in Table 1
LOOKUP = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.00]


def paired_long() -> pd.DataFrame:
    """One row per module per replicate: metric, MAG encoders, MAG breadth."""
    enc, cov, _ = C.mag_tables()
    det = cov[C.ENRICHMENT] >= C.MAG_DETECT

    n_enc = (enc >= C.ENCODING).astype(int) @ det.astype(int)   # module x sample
    n_present = det.sum(axis=0)

    _, cnt, names = C.metric_matrices("enrichments")
    idx = cnt.index.intersection(n_enc.index)

    rows = []
    for s in C.ENRICHMENT:
        rows.append(pd.DataFrame({
            "module": idx, "sample": s,
            "count_metric": cnt.loc[idx, s].values,
            "mag_encoders": n_enc.loc[idx, s].values,
            "mags_detected": int(n_present[s]),
        }))
    long = pd.concat(rows, ignore_index=True)
    long["mag_breadth"] = long.mag_encoders / long.mags_detected
    long["module_name"] = long.module.map(names["module_name"])
    return long


def run_calibration() -> None:
    C.style()
    long = paired_long()
    long.to_csv(C.OUT / "calibration_long.csv", index=False)
    print(f"{long.module.nunique()} modules x {long['sample'].nunique()} replicates "
          f"= {len(long)} paired observations")
    print(f"MAGs detected per replicate: {long.groupby('sample').mags_detected.first().min()}"
          f" to {long.groupby('sample').mags_detected.first().max()}")

    rho = spearmanr(long.count_metric, long.mag_encoders)
    print(f"\nSpearman, metric vs integer MAG encoders: rho = {rho.statistic:.3f}, "
          f"p = {rho.pvalue:.2e}")

    # --- binned calibration ---------------------------------------------------
    edges = [0, .025, .05, .075, .10, .15, .20, .25, .30, .40, .60, 1.0, 10]
    long["bin"] = pd.cut(long.count_metric, edges, include_lowest=True)
    tab = (long.groupby("bin", observed=True)
           .agg(n_obs=("mag_encoders", "size"),
                median_metric=("count_metric", "median"),
                mean_MAG_encoders=("mag_encoders", "mean"),
                median_MAG_encoders=("mag_encoders", "median"),
                q25=("mag_encoders", lambda x: x.quantile(.25)),
                q75=("mag_encoders", lambda x: x.quantile(.75)),
                mean_MAG_breadth=("mag_breadth", "mean"))
           .round(3))
    tab.to_csv(C.OUT / "calibration_bins.csv")
    print("\n--- calibration by metric bin ---")
    print(tab.to_string())

    # --- the lookup a reader wants -------------------------------------------
    # robust straight line through the cloud (Theil-Sen, insensitive to the
    # long right tail of multi-copy modules)
    sl, ic, lo, hi = theilslopes(long.mag_encoders, long.count_metric)
    print(f"\nTheil-Sen: MAG encoders = {sl:.2f} x metric + {ic:.2f}  "
          f"(slope 95% CI {lo:.2f} to {hi:.2f})")

    look = []
    for v in LOOKUP:
        near = long[(long.count_metric - v).abs() < max(0.02, v * .25)]
        look.append({
            "count_metric": v,
            "n_obs_nearby": len(near),
            "median_MAG_encoders": near.mag_encoders.median() if len(near) else np.nan,
            "IQR": (f"{near.mag_encoders.quantile(.25):.0f}-{near.mag_encoders.quantile(.75):.0f}"
                    if len(near) else "--"),
            "mean_MAG_breadth": round(near.mag_breadth.mean(), 3) if len(near) else np.nan,
            "TheilSen_fit": round(sl * v + ic, 2),
        })
    look = pd.DataFrame(look)
    look.to_csv(C.OUT / "calibration_lookup.csv", index=False)
    print("\n--- lookup: metric value -> whole encoding MAGs ---")
    print(look.to_string(index=False))

    # The straight line is only usable up to about 0.5: MAG encoders are capped
    # at the MAGs detected, so the relationship saturates and the fit
    # overshoots above that. Quote the empirical medians instead.
    def near(v):
        w = long[(long.count_metric - v).abs() < max(0.02, v * .25)]
        return w.mag_encoders.median(), w.mag_encoders.quantile(.25), w.mag_encoders.quantile(.75)
    m20, l20, h20 = near(0.20)
    m40, l40, h40 = near(0.40)
    print(f"\nThe question asked: on the count-based metric, 0.20 corresponds to a "
          f"median of {m20:.0f} encoding MAGs (IQR {l20:.0f} to {h20:.0f}) and 0.40 to "
          f"{m40:.0f} (IQR {l40:.0f} to {h40:.0f}). So a move from 0.2 to 0.4 is about "
          f"{m40-m20:.0f} more genomes out of the {long.mags_detected.min()} to "
          f"{long.mags_detected.max()} detected per replicate.")
    print("The Theil-Sen line fits below about 0.5 and overshoots above it, because "
          "the number of MAGs is capped; use the binned medians beyond that.")

    # --- figure ---------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    ax = axes[0]
    jit = np.random.default_rng(0).normal(0, .10, len(long))
    ax.scatter(long.count_metric, long.mag_encoders + jit, s=9, alpha=.22,
               c=C.GROUP_COLOR["Enrichment"], lw=0)
    xs = np.linspace(0, long.count_metric.quantile(.995), 100)
    ax.plot(xs, sl * xs + ic, color=C.CLASS_COLOR["count-based only"], lw=2,
            label=f"Theil-Sen: {sl:.1f}x + {ic:.1f}")
    b = tab[tab.n_obs >= 20]
    ax.errorbar(b.median_metric, b.median_MAG_encoders,
                yerr=[b.median_MAG_encoders - b.q25, b.q75 - b.median_MAG_encoders],
                fmt="o", ms=5, lw=1.4, capsize=3, color="#37474f",
                label="binned median, IQR", zorder=4)
    ax.set_xlim(0, long.count_metric.quantile(.995))
    ax.set_xlabel("count-based metric (copy number / populations)")
    ax.set_ylabel("MAGs encoding the module (integer)")
    ax.set_title("A   Calibration against whole genomes", loc="left")
    ax.legend()

    ax = axes[1]
    ax.scatter(long.count_metric, long.mag_breadth, s=9, alpha=.22,
               c=C.GROUP_COLOR["Enrichment"], lw=0)
    ax.plot([0, 1], [0, 1], ls="--", lw=.9, color="#888888")
    ax.text(.62, .68, "1:1", color="#888888", fontsize=8, rotation=38)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("count-based metric")
    ax.set_ylabel("fraction of detected MAGs encoding it")
    ax.set_title("B   Fraction against fraction", loc="left")

    C.save(fig, "fig_calibration")


# ============================================================================
# Supplementary Figure S15: the count-based metric against encoding MAGs
# ----------------------------------------------------------------------------
# For every module in every enrichment replicate, the number of detected MAGs encoding it (pathwise
# completeness >= 0.75, mean coverage >= 1) is set against each contig metric.
# Within a replicate the number of detected MAGs is fixed, so Spearman rho
# against the MAG count equals rho against MAG breadth (mag_control_concordance).
# Panel D shows how much of the cyanobacterial community the MAGs capture,
# which bounds every MAG-based attribution (Fig 3E, Table S4).
# ============================================================================

def run_mag_figure() -> None:
    C.style()
    long = paired_long()
    cov, _, _ = C.metric_matrices("enrichments")
    long["coverage_metric"] = [cov.at[m, s] if m in cov.index else np.nan
                               for m, s in zip(long.module, long["sample"])]
    long = long.dropna(subset=["coverage_metric"])
    conc = pd.read_csv(C.OUT / "mag_control_concordance.csv").set_index("sample")
    rep = C.mag_representation()

    def binned(x, edges):
        b = pd.cut(x, edges, include_lowest=True)
        g = long.groupby(b, observed=True).mag_encoders
        out = pd.DataFrame({"x": x.groupby(b, observed=True).median(), "med": g.median(),
                            "q25": g.quantile(.25), "q75": g.quantile(.75), "n": g.size()})
        return out[out.n >= 20]

    fig, axes = plt.subplots(2, 2, figsize=(10, 8.2))
    jit = np.random.default_rng(0).normal(0, .12, len(long))
    col = C.GROUP_COLOR["Enrichment"]

    # A: count-based metric
    ax = axes[0, 0]
    xmax = long.count_metric.quantile(.995)
    ax.scatter(long.count_metric, long.mag_encoders + jit, s=8, alpha=.22, c=col, lw=0)
    b = binned(long.count_metric, [0, .025, .05, .075, .10, .15, .20, .25, .30, .40, .60, 1.0, 10])
    ax.errorbar(b.x, b.med, yerr=[b.med - b.q25, b.q75 - b.med], fmt="o-", ms=5, lw=1.4,
                capsize=3, color=C.CLASS_COLOR["count-based only"], label="binned median, IQR", zorder=4)
    ax.set_xlim(0, xmax)
    ax.set_xlabel("count-based metric (copy number per population)")
    ax.set_ylabel("MAGs encoding the module")
    r = conc.rho_count_vs_MAG_breadth
    ax.text(.97, .05, f"Spearman $\\rho$ = {r.min():.2f}–{r.max():.2f}\n(per replicate)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8)
    ax.legend(loc="upper left")
    ax.set_title("A   Count-based metric", loc="left")

    # B: coverage-based metric
    ax = axes[0, 1]
    x = long.coverage_metric
    ax.scatter(x, long.mag_encoders + jit[:len(long)], s=8, alpha=.22, c=col, lw=0)
    ax.set_xscale("symlog", linthresh=1)
    edges = [0, 1e-9] + list(np.logspace(0, np.log10(x.max() * 1.01), 9))
    b = binned(x, edges)
    ax.errorbar(b.x, b.med, yerr=[b.med - b.q25, b.q75 - b.med], fmt="o-", ms=5, lw=1.4,
                capsize=3, color=C.CLASS_COLOR["both metrics"], label="binned median, IQR", zorder=4)
    ax.set_xlabel("coverage-based metric (mean module coverage)")
    ax.set_ylabel("MAGs encoding the module")
    r = conc.rho_coverage_vs_MAG_breadth
    ax.text(.97, .95, f"Spearman $\\rho$ = {r.min():.2f}–{r.max():.2f}\n(per replicate)",
            transform=ax.transAxes, ha="right", va="top", fontsize=8)
    ax.legend(loc="upper left")
    ax.set_title("B   Coverage-based metric", loc="left")

    # C: per-replicate rho, paired
    ax = axes[1, 0]
    order = C.ENRICHMENT
    xs = np.arange(len(order))
    for i, s in enumerate(order):
        ax.plot([i, i], [conc.at[s, "rho_coverage_vs_MAG_breadth"], conc.at[s, "rho_count_vs_MAG_breadth"]],
                color="#bbbbbb", lw=1, zorder=1)
    ax.scatter(xs, conc.loc[order, "rho_count_vs_MAG_breadth"], s=40, zorder=3,
               color=C.CLASS_COLOR["count-based only"], label="count-based")
    ax.scatter(xs, conc.loc[order, "rho_coverage_vs_MAG_breadth"], s=40, zorder=3,
               color=C.CLASS_COLOR["both metrics"], label="coverage-based")
    for i, s in enumerate(order):
        ax.text(i, 0.03, f"{int(conc.at[s, 'n_MAGs_detected'])}", ha="center", fontsize=7, color="#555555")
    ax.text(-0.6, 0.03, "MAGs:", ha="right", fontsize=7, color="#555555")
    ax.set_xticks(xs, order, rotation=0, fontsize=8)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Spearman $\\rho$ with MAGs encoding")
    ax.legend(loc="upper right", ncol=2)
    ax.set_title("C   Agreement with genomes, per replicate", loc="left")

    # D: cyanobacterial capture by MAGs
    ax = axes[1, 1]
    v = rep.reindex(order)
    ax.bar(xs, v, color=[C.CYANO_C if a >= 25 else "#c8e6c9" for a in v], edgecolor="#333333", lw=.5)
    for i, a in enumerate(v):
        ax.text(i, a + 2, f"{a:.0f}", ha="center", fontsize=7.5)
    ax.axhline(25, ls="--", lw=.8, color="#888888")
    ax.set_xticks(xs, order, fontsize=8)
    ax.set_ylim(0, 110)
    ax.set_ylabel("cyanobacterial RpL19 coverage\ncaptured by MAGs (%)")
    ax.set_title("D   Cyanobacteria recovered as MAGs", loc="left")

    fig.tight_layout()
    C.save(fig, "fig_si_count_vs_mags")
    print("rho count", conc.rho_count_vs_MAG_breadth.round(2).tolist())
    print("rho cov", conc.rho_coverage_vs_MAG_breadth.round(2).tolist())
    print("capture", v.round(0).to_dict())


def main() -> None:
    run_two_metric()
    run_mag_control()
    run_table1()
    run_calibration()
    run_mag_figure()


if __name__ == "__main__":
    main()
