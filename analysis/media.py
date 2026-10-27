"""
The single-inoculum media experiment (V, NV, NVNN): two-metric contrasts,
set-level tests, Supplementary Table S7 and the Fig 4 analysis.

Run:  python media.py      (from analysis/, after the modules before it in
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

# ============================================================================
# Two-metric analysis of the media experiment: run_media()
# ----------------------------------------------------------------------------
# Three questions:
#   1. Which functions respond when nitrate is withdrawn, and which when vitamins
#      are withdrawn?
#   2. Are anoxic and microaerobic functions among them? Individual modules have
#      little power at n = 4, so this is asked of a pre-specified set of modules
#      (common.ANOXIC_MODULES) as a set.
#   3. Do the modules classified as community-derived in the long-term enrichments
#      respond to a change in the medium?
#
# Produces
#   outputs/table2_media_two_metric.{csv,tex}
#   outputs/media_set_tests.csv
#   outputs/fig4.{pdf,svg}
#
# Run:  python media.py
# ============================================================================



NITRATE_REPLETE = C.MEDIA_V + C.MEDIA_NV


GROUPS = {"V": C.MEDIA_V, "NV": C.MEDIA_NV, "NVNN": C.MEDIA_NVNN}


def load_niche_class() -> pd.Series:
    """The enrichment-experiment classification, written by function.py."""
    f = C.OUT / "module_niche_classification.csv"
    if not f.exists():
        raise SystemExit("run function.py first")
    return pd.read_csv(f, index_col=0)["niche_class"]


def contrasts(marker: str = "rps16"):
    """
    RpS16 is the denominator for the media experiment, the marker used for its
    composition. RpL19 gives the same cyanobacterial composition once one RpS16
    cluster, a eukaryotic plastid placed inside the cyanobacteria, is set aside.
    """
    cov_m, cnt_m, names = C.metric_matrices("media", marker)
    out = {}
    for tag, g1, g2 in [("nitrate", C.MEDIA_NVNN, NITRATE_REPLETE),
                        ("vitamin", C.MEDIA_V, C.MEDIA_NV),
                        ("nvnn_vs_nv", C.MEDIA_NVNN, C.MEDIA_NV)]:
        out[tag] = {"cov": C.differential(cov_m, g1, g2),
                    "cnt": C.differential(cnt_m, g1, g2)}
    return out, cov_m, cnt_m, names


def assemble(res, names, cls) -> pd.DataFrame:
    t = pd.DataFrame(index=names.index)
    for tag, d in res.items():
        for metric, tab in d.items():
            t[f"{tag}_{metric}_l2fc"] = tab.l2fc
            t[f"{tag}_{metric}_p"] = tab.p
            t[f"{tag}_{metric}_padj"] = tab.padj
            t[f"{tag}_{metric}_sep"] = tab.separated
    t = t.join(names)
    t["enrichment_niche_class"] = cls.reindex(t.index)
    t["anoxic_set"] = t.index.isin(C.ANOXIC_MODULES)
    t["oxygen_set"] = t.index.isin(C.OXYGEN_MODULES)
    return t


def write_table2(t: pd.DataFrame) -> None:
    nit = t.dropna(subset=["nitrate_cov_l2fc"])
    sig = lambda m: (nit[f"nitrate_{m}_padj"] < C.FDR) & (nit[f"nitrate_{m}_l2fc"].abs() > C.MIN_L2FC)
    responders = nit[sig("cov") | sig("cnt") | (nit.nitrate_cov_sep == True)
                     | (nit.nitrate_cnt_sep == True)]
    order = responders.nitrate_cnt_l2fc.fillna(responders.nitrate_cov_l2fc)
    responders = responders.reindex(order.sort_values(ascending=False).index)
    community = [m for m in t.index
                 if t.loc[m, "enrichment_niche_class"] == "count-based only"
                 and m not in set(responders.index)]

    def block(idx):
        lines = []
        for m in idx:
            r = t.loc[m]

            def cell(v, q, sep):
                if pd.isna(v):
                    return "--", "--", "--"
                s = f"{v:.2f}"
                qs = f"{q:.3f}"
                if q < C.FDR and abs(v) > C.MIN_L2FC:
                    s, qs = rf"\textbf{{{s}}}", rf"\textbf{{{qs}}}"
                return s, qs, ("yes" if bool(sep) else "--")
            a1, a2, a3 = cell(r.nitrate_cov_l2fc, r.nitrate_cov_padj, r.nitrate_cov_sep)
            b1, b2, b3 = cell(r.nitrate_cnt_l2fc, r.nitrate_cnt_padj, r.nitrate_cnt_sep)
            raw = str(r.module_name)
            if len(raw) > 44:
                raw = raw[:42].rsplit(" ", 1)[0].rstrip(",") + " ..."
            name = C.tex_escape(raw).replace(" ...", r"\ldots")
            klass = {"both metrics": "both", "count-based only": "count only",
                     "coverage-based only": "coverage only",
                     "neither": "neither"}.get(r.enrichment_niche_class, "n.d.")
            lines.append(f"{m} & {name} & {a1} & {a2} & {a3} & {b1} & {b2} & {b3} & {klass} \\\\")
        return lines

    tex = r"""\begin{table}[h]
\centering
\caption{\textbf{How KEGG modules responded to nitrate withdrawal in the single-inoculum
experiment, on both metrics.} Nitrate-free replicates (NVNN, $n=4$) are compared
with the pooled nitrate-replete conditions (V and NV, $n=9$); no module differed
between V and NV on either metric. Fold changes use a per-module half-minimum
pseudocount; $q$ is the Benjamini--Hochberg adjustment of a two-sided exact
Wilcoxon rank-sum test (all 715 group assignments enumerated) across the modules
tested on that metric (217 coverage-based, 182 count-based; count-based metric
normalised by RpS16 populations). Bold marks $q<0.05$ with $|\log_2\mathrm{FC}|>0.5$.
``Sep.'' marks modules for which every nitrate-free replicate lies outside the range
of every nitrate-replete replicate. Rows are the modules significant or completely
separated on either metric, ordered by count-based fold change. The final column
gives how the same module behaved across the eight long-term enrichments relative
to the reservoir; the lower block lists the remaining modules that were increased on
the count-based metric alone in those enrichments.}
\label{tab:media_two_metric}
\footnotesize
\begin{tabular}{@{}llrrcrrcl@{}}
\toprule
 &  & \multicolumn{3}{c}{\textbf{Coverage-based}} & \multicolumn{3}{c}{\textbf{Count-based}} & \textbf{Enrichment} \\
\cmidrule(l){3-5}\cmidrule(l){6-8}
\textbf{Module} & \textbf{Module name} & \textbf{log\textsubscript{2}FC} & \textbf{$q$} & \textbf{Sep.} & \textbf{log\textsubscript{2}FC} & \textbf{$q$} & \textbf{Sep.} & \textbf{class} \\
\midrule
""" + "\n".join(block(responders.index)) + r"""
\midrule
\multicolumn{9}{@{}l}{\emph{Other modules increased on the count-based metric alone in the long-term enrichments}} \\
""" + "\n".join(block(community)) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    (C.OUT / "table2_media_two_metric.tex").write_text(tex)
    print(f"wrote outputs/table2_media_two_metric.tex "
          f"({len(responders)} responders + {len(community)} community-derived)")


# ---------------------------------------------------------------- the figure
def panel_C(ax, cov_m):
    samples = [s for g in ("V", "NV", "NVNN") for s in GROUPS[g] if s in cov_m.columns]
    comp = C.compositional(cov_m)
    X = StandardScaler().fit_transform(comp[samples].T.values)
    p = PCA(n_components=2).fit(X)
    sc = p.transform(X)
    var = p.explained_variance_ratio_ * 100
    for g in ("V", "NV", "NVNN"):
        idx = [i for i, s in enumerate(samples) if s in GROUPS[g]]
        ax.scatter(sc[idx, 0], sc[idx, 1], s=70, color=C.GROUP_COLOR[g],
                   edgecolor="white", lw=0.8, label=g, zorder=3)
    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    F, p = C.permanova(comp, GROUPS)
    pairs = {}
    for a, b in (("V", "NV"), ("V", "NVNN"), ("NV", "NVNN")):
        pairs[f"{a} vs {b}"] = C.permanova(comp, {a: GROUPS[a], b: GROUPS[b]})
    ax.set_title("C   Functional community structure", loc="left")
    txt = (f"PERMANOVA $F$ = {F:.2f}, $p$ = {p:.3f}\n"
           + "\n".join(f"{k}: $p$ = {v[1]:.3f}" for k, v in pairs.items()))
    ax.annotate(txt, xy=(0.02, 0.98), xycoords="axes fraction", fontsize=6.8,
                color="#555555", va="top")
    ax.legend(title="Medium", loc="lower right")
    return var, F, p, pairs


def panel_D(ax, t):
    x = t.nitrate_cov_l2fc.values
    y = -np.log10(t.nitrate_cov_p.clip(1e-4).values)
    grey = "#c8ccd0"
    col = np.where(t.oxygen_set, "#b71c1c", np.where(t.index.isin(C.N_ANAEROBIC_MODULES), "#e07b39", grey))
    ax.scatter(x, y, s=np.where(t.anoxic_set, 34, 16), c=col,
               edgecolor="none", alpha=0.9, zorder=3)
    ax.axhline(-np.log10(0.05), color="#b0392b", lw=0.8, ls="--")
    ax.annotate("$p$ = 0.05 (unadjusted)", xy=(ax.get_xlim()[0], -np.log10(0.05)),
                fontsize=6.5, va="bottom", color="#b0392b")
    for m, off in [("M00175", (6, 2)), ("M00930", (-58, 3)), ("M00931", (-58, -7)),
                   ("M00595", (6, 3)), ("M00597", (6, -8)), ("M00156", (6, 2)),
                   ("M00124", (6, -8))]:
        if m in t.index and not pd.isna(t.loc[m, "nitrate_cov_l2fc"]):
            ax.annotate(str(t.loc[m, "module_name"]).split(",")[0],
                        (t.loc[m, "nitrate_cov_l2fc"],
                         -np.log10(max(t.loc[m, "nitrate_cov_p"], 1e-4))),
                        fontsize=6.3, xytext=off, textcoords="offset points")
    ax.set_xlabel("log$_2$ fold change, nitrate-free / nitrate-replete (coverage-based)")
    ax.set_ylabel("$-$log$_{10}$ $p$")
    ax.set_title("D   Response to nitrate withdrawal", loc="left")
    ax.legend(handles=[mpatches.Patch(color="#b71c1c", label="oxygen-related (anoxic or microaerobic)"),
                       mpatches.Patch(color="#e07b39", label="anaerobic nitrogen metabolism"),
                       mpatches.Patch(color=grey, label="all other modules")],
              loc="upper left", fontsize=7)


def panel_E(ax, t, tests):
    """The pre-specified anoxic set against everything else, on both metrics."""
    data, labels, colours = [], [], []
    for metric, mlabel in (("cov", "coverage-based"), ("cnt", "count-based")):
        col = f"nitrate_{metric}_l2fc"
        data.append(t.loc[t.oxygen_set, col].dropna().values)
        labels.append(f"oxygen-related set\n({mlabel})")
        colours.append("#b71c1c")
        data.append(t.loc[~t.oxygen_set, col].dropna().values)
        labels.append(f"all other\n({mlabel})")
        colours.append("#9aa0a6")
    pos = [1, 1.85, 3.4, 4.25]
    bp = ax.boxplot(data, positions=pos, widths=0.62, patch_artist=True,
                    showfliers=False, medianprops=dict(color="black", lw=1.2))
    for patch, c in zip(bp["boxes"], colours):
        patch.set_facecolor(c)
        patch.set_alpha(0.30)
        patch.set_edgecolor(c)
    rng = np.random.default_rng(1)
    for p_, d, c in zip(pos, data, colours):
        ax.scatter(p_ + rng.uniform(-0.16, 0.16, len(d)), d, s=10, color=c,
                   alpha=0.75, edgecolor="none", zorder=3)
    ax.axhline(0, color="#888888", lw=0.7, ls=":")
    ax.set_xticks(pos)
    ax.set_xticklabels([l.split("\n")[0] for l in labels], fontsize=7.5)
    for x0, lab in ((np.mean(pos[:2]), "coverage-based metric"),
                    (np.mean(pos[2:]), "count-based metric")):
        ax.annotate(lab, xy=(x0, -0.13), xycoords=("data", "axes fraction"),
                    ha="center", fontsize=8, fontweight="bold")
    ax.set_ylabel("log$_2$ fold change, nitrate-free / nitrate-replete")
    ax.set_title("E   Oxygen-related functions, tested as a pre-specified set", loc="left")
    for (x0, x1), row in zip([(pos[0], pos[1]), (pos[2], pos[3])], tests):
        yy = max(np.max(data[0]), np.max(data[2])) * 1.05
        ax.plot([x0, x0, x1, x1], [yy, yy * 1.04, yy * 1.04, yy], lw=0.8, color="#444444")
        ax.text((x0 + x1) / 2, yy * 1.06, f"$p$ = {row['p']:.3f}",
                ha="center", fontsize=7)


def placeholder(ax, letter, title, note="schematic, composited in Affinity"):
    ax.add_patch(mpatches.FancyBboxPatch((0.03, 0.03), 0.94, 0.94,
                 transform=ax.transAxes, boxstyle="round,pad=0.01",
                 facecolor="#f6f6f6", edgecolor="#c2c2c2", lw=1.0, ls=(0, (5, 4))))
    ax.text(0.5, 0.55, title, transform=ax.transAxes, ha="center", va="center",
            fontsize=9, color="#8a8a8a")
    ax.text(0.5, 0.40, note, transform=ax.transAxes, ha="center", va="center",
            fontsize=7, color="#b0b0b0", style="italic")
    ax.set_title(letter, loc="left")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def run_media() -> None:
    C.style()
    cls = load_niche_class()
    res, cov_m, cnt_m, names = contrasts("rps16")
    t = assemble(res, names, cls)
    t.to_csv(C.OUT / "table2_media_two_metric.csv")

    print("significant after FDR correction:")
    for tag in ("nitrate", "vitamin", "nvnn_vs_nv"):
        for metric in ("cov", "cnt"):
            n = int(((t[f"{tag}_{metric}_padj"] < C.FDR) &
                     (t[f"{tag}_{metric}_l2fc"].abs() > C.MIN_L2FC)).sum())
            print(f"  {tag:11s} {metric}: {n}")

    print("\ncompletely separated on both metrics, nitrate axis:")
    both_sep = t[(t.nitrate_cov_sep == True) & (t.nitrate_cnt_sep == True)]
    print(both_sep[["module_name", "nitrate_cov_l2fc", "nitrate_cov_p",
                    "nitrate_cnt_l2fc", "nitrate_cnt_p"]]
          .to_string(float_format=lambda v: f"{v:7.3f}"))

    # set-level tests
    tests = []
    SETS = [("oxygen-related", set(C.OXYGEN_MODULES)),
            ("nitrogen, anaerobic", set(C.N_ANAEROBIC_MODULES)),
            ("all anoxic/microaerobic", set(C.ANOXIC_MODULES))]
    for axis in ("nitrate", "vitamin"):
        for label, members in SETS:
            for metric in ("cov", "cnt"):
                tab = pd.DataFrame({"l2fc": t[f"{axis}_{metric}_l2fc"]}).dropna()
                r = C.set_test(tab, members, f"{label}, {axis}, {metric}")
                tests.append(r)
    st = pd.DataFrame(tests)
    st.to_csv(C.OUT / "media_set_tests.csv", index=False)
    print("\nset-level tests (pre-specified anoxic and microaerobic modules):")
    print(st.to_string(index=False, float_format=lambda v: f"{v:7.3f}"))

    print("\nper-module, the anoxic set on the nitrate axis:")
    an = t[t.anoxic_set][["module_name", "nitrate_cov_l2fc", "nitrate_cov_p",
                          "nitrate_cnt_l2fc", "nitrate_cnt_p", "enrichment_niche_class"]]
    print(an.sort_values("nitrate_cov_l2fc", ascending=False)
          .to_string(float_format=lambda v: f"{v:7.2f}"))

    # anvi'o assembly-level L19 denominator, as a check (no RpL19 clusters here)
    res19, _, _, _ = contrasts("anvio_l19")
    a = res["nitrate"]["cnt"].l2fc
    b = res19["nitrate"]["cnt"].l2fc
    common_idx = a.index.intersection(b.index)
    print(f"\ncount-based nitrate log2FC, RpS16 against RpL19 denominator: "
          f"Pearson r = {np.corrcoef(a[common_idx], b[common_idx])[0,1]:.3f} "
          f"over {len(common_idx)} modules")

    write_table2(t)

    fig = plt.figure(figsize=(11.5, 11))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.42, 0.75, 0.75], hspace=0.42, wspace=0.26)
    placeholder(fig.add_subplot(gs[0, 0]), "A",
                "Single-inoculum media experiment\n(V, NV, NVNN; 5 replicates; P10)")
    placeholder(fig.add_subplot(gs[0, 1]), "B",
                "EcoPhylo abundance heatmap, RpS16 clusters",
                "anvi-interactive, composited in Affinity")
    var, F, p, pairs = panel_C(fig.add_subplot(gs[1, 0]), cov_m)
    panel_D(fig.add_subplot(gs[1, 1]), t)
    ox = [r for r in tests if r["set"].startswith("oxygen-related, nitrate")]
    panel_E(fig.add_subplot(gs[2, :]), t, ox)
    C.save(fig, "fig4")
    print(f"\npanel C: PC1={var[0]:.1f}% PC2={var[1]:.1f}%, "
          f"three-group PERMANOVA F={F:.2f} p={p:.3f}")
    for k, v in pairs.items():
        print(f"  {k}: F={v[0]:.2f} p={v[1]:.3f}")


# ============================================================================
# Trajectories across passages (RpS16, all passages)
# ----------------------------------------------------------------------------
# The same lineages were sequenced at three points: the parent culture of each
# medium (16 Aug 2023, after two pre-culture passages), an intermediate passage
# (25 Jan 2024, replicates R1-R5 plus the continuing parent line) and P10.
# Unlike source vs enrichment, these are paired in time, so they show when the
# media diverge and when the filamentous cyanobacteria come in.
#     outputs/media_trajectories.pdf / .csv
# ============================================================================
from scipy.spatial.distance import pdist, squareform  # noqa: E402

MED_COL = {"V": "#6aae6a", "NV": "#f28468", "NVNN": "#e377b8"}
STAGES = ["parent", "mid", "P10"]
CYANO_GROUPS = {"Cyanobiaceae": "picocyanobacteria (Cyanobiaceae)",
                "Nostocaceae": "Calothrix (Nostocaceae)",
                "Leptolyngbyaceae": "Leptolyngbya (Leptolyngbyaceae)"}
CYANO_COL = {"picocyanobacteria (Cyanobiaceae)": "#4f9fd1", "Calothrix (Nostocaceae)": "#c2185b",
             "Leptolyngbya (Leptolyngbyaceae)": "#f2b134", "other Cyanobacteriota": "#9e9e9e"}


def passage_table():
    meta = pd.read_csv(C.DATA / "media_samples.csv")
    meta["stage"] = np.where(meta.passage.eq("P10"), "P10",
                             np.where(meta.passage.eq("parent"), "parent", "mid"))
    meta["lineage"] = np.where(meta.replicate.notna() & meta.stage.ne("parent"),
                               meta.condition + "_R" + meta.replicate.fillna(0).astype(int).astype(str),
                               meta.condition + "_line")
    e = pd.read_csv(C.DATA / "ecophylo_coverage.csv")
    e = e[e.run == "rps16_media_allpassages"]
    cov = e.pivot_table(index="cluster", columns="sample", values="mean_coverage_Q2Q3").fillna(0)
    rel = cov / cov.sum()
    meta = meta.set_index("sample").loc[rel.columns]
    return rel, meta


def run_trajectories():
    rel, meta = passage_table()
    D = squareform(pdist(rel.T.values, "braycurtis"))
    n = len(D)
    J = np.eye(n) - 1 / n
    B = -0.5 * J @ (D ** 2) @ J
    w, v = np.linalg.eigh(B)
    o = np.argsort(w)[::-1][:2]
    xy = v[:, o] * np.sqrt(w[o])
    expl = w[o] / w[w > 0].sum() * 100
    meta[["PCo1", "PCo2"]] = xy
    Dm = pd.DataFrame(D, index=rel.columns, columns=rel.columns)

    # between-media distance at each stage, and distance of each lineage to its parent
    rows = []
    for st in STAGES:
        for a, b in (("V", "NV"), ("V", "NVNN"), ("NV", "NVNN")):
            ia = meta.index[(meta.condition == a) & (meta.stage == st)]
            ib = meta.index[(meta.condition == b) & (meta.stage == st)]
            for x in ia:
                for y in ib:
                    rows.append(dict(stage=st, pair=f"{a}-{b}", distance=Dm.loc[x, y]))
    between = pd.DataFrame(rows)
    cy = pd.read_csv(C.DATA / "ecophylo_clusters.csv")
    cy = cy[(cy.run == "rps16_media_allpassages") & (cy.t_phylum == C.CYANO_PHYLUM)].drop_duplicates("cluster").set_index("cluster")
    grp = cy.t_family.map(CYANO_GROUPS).fillna("other Cyanobacteriota")
    cyano = rel.loc[rel.index.intersection(grp.index)].groupby(grp).sum().T

    # per lineage: distance to its parent, and within-medium spread, over time
    to_parent, spread = [], []
    for cond in MED_COL:
        par = meta.index[(meta.condition == cond) & (meta.stage == "parent")][0]
        for st in ("mid", "P10"):
            ix = meta.index[(meta.condition == cond) & (meta.stage == st) & meta.lineage.str.contains("_R")]
            to_parent += [dict(condition=cond, stage=st, lineage=meta.loc[x, "lineage"], distance=Dm.loc[par, x]) for x in ix]
            spread += [dict(condition=cond, stage=st, distance=Dm.loc[x, y]) for i, x in enumerate(ix) for y in ix[i + 1:]]
    to_parent, spread = pd.DataFrame(to_parent), pd.DataFrame(spread)

    fig, axes = plt.subplots(1, 4, figsize=(15, 3.9), gridspec_kw=dict(width_ratios=[1.15, 0.8, 0.8, 1.5], wspace=0.38))
    ax = axes[0]
    for lin, d in meta.groupby("lineage"):
        cond = d.condition.iloc[0]
        par = meta[(meta.condition == cond) & (meta.stage == "parent")]
        path = pd.concat([par, d]).drop_duplicates()
        path = path.assign(o=path.stage.map({s: i for i, s in enumerate(STAGES)})).sort_values("o")
        ax.plot(path.PCo1, path.PCo2, color=MED_COL[cond], lw=0.8, alpha=0.6, zorder=1)
    mk = {"parent": "*", "mid": "o", "P10": "s"}
    sz = {"parent": 160, "mid": 26, "P10": 34}
    for (cond, st), d in meta.groupby(["condition", "stage"]):
        ax.scatter(d.PCo1, d.PCo2, marker=mk[st], s=sz[st], color=MED_COL[cond], edgecolor="#333333", lw=0.5, zorder=3)
    hs = [plt.Line2D([], [], color=c, lw=2, label=g) for g, c in MED_COL.items()]
    hs += [plt.Line2D([], [], marker=mk[s], ls="", color="#888888", markeredgecolor="#333333",
                      markersize=7 if s != "parent" else 11,
                      label={"parent": "parent (Aug 2023)", "mid": "Jan 2024", "P10": "P10"}[s]) for s in STAGES]
    ax.legend(handles=hs, fontsize=6.3, frameon=False, loc="best")
    ax.set_xlabel(f"PCoA 1 ({expl[0]:.0f}%)"); ax.set_ylabel(f"PCoA 2 ({expl[1]:.0f}%)")
    ax.set_title("Lineage trajectories (RpS16)", fontsize=9, fontweight="bold")

    x = {"parent": 0, "mid": 1, "P10": 2}
    for ax, tab, ttl, yl in ((axes[1], to_parent, "Distance from the parent", "Bray-Curtis to parent culture"),
                             (axes[2], spread, "Spread among replicates", "Bray-Curtis between replicates")):
        for k, cond in enumerate(MED_COL):
            d = tab[tab.condition == cond]
            if "lineage" in d:
                for lin, dd in d.groupby("lineage"):
                    ax.plot([0] + [x[s] + (k - 1) * 0.06 for s in dd.stage], [0] + list(dd.distance),
                            color=MED_COL[cond], lw=0.6, alpha=0.45)
            m = d.groupby("stage").distance.mean().reindex(["mid", "P10"])
            xs = [x[s] + (k - 1) * 0.06 for s in m.index]
            ax.plot(([0] if "lineage" in d else []) + xs, ([0] if "lineage" in d else []) + list(m.values),
                    color=MED_COL[cond], lw=2.2, marker="o", ms=4, label=cond)
        ax.set_xticks([0, 1, 2], ["parent", "Jan 2024", "P10"])
        ax.set_ylim(0, 1.02); ax.set_ylabel(yl, fontsize=8)
        ax.set_title(ttl, fontsize=9, fontweight="bold")
    axes[1].legend(fontsize=7, frameon=False, loc="lower right")

    ax = axes[3]
    order = meta.assign(o=meta.stage.map({s: i for i, s in enumerate(STAGES)}),
                        c=meta.condition.map({"V": 0, "NV": 1, "NVNN": 2})).sort_values(["c", "o", "lineage"]).index
    bottom = np.zeros(len(order))
    for g, col in CYANO_COL.items():
        if g in cyano.columns:
            vals = cyano.reindex(order)[g].fillna(0).values
            ax.bar(range(len(order)), vals, bottom=bottom, color=col, width=0.85, label=g)
            bottom += vals
    ax.set_xticks(range(len(order)), [f"{meta.loc[s, 'condition']} {dict(parent='parent', mid='Jan', P10='P10')[meta.loc[s, 'stage']]}"
                                      for s in order], rotation=90, fontsize=5.5)
    ax.set_ylabel("Cyanobacterial share of RpS16 coverage", fontsize=8)
    ax.legend(fontsize=6.3, frameon=False, loc="upper left")
    ax.set_title("Cyanobacteria through the passages", fontsize=9, fontweight="bold")
    for a in fig.axes:
        a.spines[["top", "right"]].set_visible(False)
    for a, l in zip(fig.axes, "ABCD"):
        a.text(-0.12, 1.04, l, transform=a.transAxes, fontsize=12, fontweight="bold")
    fig.savefig(C.OUT / "media_trajectories.pdf", bbox_inches="tight")
    fig.savefig(C.OUT / "media_trajectories.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    meta.join(cyano).to_csv(C.OUT / "media_trajectories.csv")
    to_parent.to_csv(C.OUT / "media_trajectory_to_parent.csv", index=False)
    summ = between.groupby(["pair", "stage"]).distance.agg(["mean", "std", "count"])
    summ.to_csv(C.OUT / "media_trajectory_between_media.csv")
    print(summ.round(3).to_string())


# ============================================================================
# Functional trajectories (media MAGs, all passages)
# ----------------------------------------------------------------------------
# All 35 media samples (P0 parent, P4 R1-R5, P5 parent line, P10 R1-R5) mapped to
# one shared genome set (pipeline/09: enrichment + media MAGs at 95% ANI). Per sample,
#     composition  relative coverage of each genome
#     function     relative module abundance: summed coverage of the genomes that
#                  encode the module (pathwise completeness >= ENCODING)
# and ask whether function settles while composition keeps changing.
#     outputs/media_functional_trajectories.pdf / .csv
# ============================================================================
KEY_MODULES = [("M00175", "Nitrogen fixation"), ("M00161", "Photosystem II"),
               ("M00165", "Calvin cycle"), ("M00123", "Biotin"), ("M00122", "Cobalamin"),
               ("M00116", "Menaquinone")]
PASSAGE_X = {"P0": 0, "P4": 4, "P5": 5, "P10": 10}


def media_sample_names() -> pd.Series:
    """M_ number -> V_P10_R1 style name (P0 parent, P4 replicates, P5 parent line, P10)."""
    m = pd.read_csv(C.DATA / "media_samples.csv").set_index("sample")
    rep = m.replicate.fillna(0).astype(int).astype(str)
    pas = np.select([m.passage.eq("P10"), m.passage.eq("parent"), m.replicate.notna()],
                    ["P10_R" + rep, "P0_parent", "P4_R" + rep], "P5_parent")
    return pd.Series(m.condition + "_" + pas, index=m.index)


def shared_reference():
    """Media reads mapped to one shared genome set (pipeline/09).

    The 88 enrichment MAGs and the 181 media MAGs, dereplicated at 95% ANI (242
    representatives), with all 35 media samples mapped competitively (bowtie2,
    CoverM mean coverage, reads >= 95% identity, genome >= 10% covered). Unlike the
    per-medium co-assembly MAGs, it holds the V and NV picocyanobacteria.
    Returns coverage (genome x sample), encoding (module x genome, pathwise >=
    ENCODING), phylum per genome, and the share of reads that aligned at all
    (bowtie2; CoverM's "unmapped" row misses pairs dropped by --no-unal).
    """
    D = C.DATA / "media_mapping"
    cov = pd.read_csv(D / "coverm_mean.tsv", sep="\t", index_col=0)
    cov.columns = [c.split(" ")[0] for c in cov.columns]
    names = media_sample_names()
    cov = cov.rename(columns=names)
    aln = pd.read_csv(D / "bowtie2_overall_alignment.txt", sep=" ", header=None, index_col=0)[1]
    aln = aln.str.rstrip("%").astype(float).rename(index=names) / 100
    q = pd.read_csv(C.DATA / "media_mag_quality.csv", index_col=0)
    mq = pd.read_csv(C.DATA / "mag_quality.tsv", sep="\t", index_col=0)
    phylum = pd.concat([q.t_phylum, mq.phylum]).reindex(cov.index)
    mods = pd.read_csv(C.DATA / "media_mag_modules.csv")
    e_media = mods.pivot_table(index="module", columns="bin", values="pathwise_module_completeness")
    e_enr = pd.read_csv(C.DATA / "mag_module_pathwise_completeness.tsv", sep="\t", index_col=0)
    e_enr = e_enr[[c for c in e_enr.columns if c in cov.index]]
    enc = (pd.concat([e_media, e_enr], axis=1).reindex(columns=cov.index) >= C.ENCODING)
    kos = pd.read_csv(D / "signature_kos.tsv", sep="\t").groupby("genome").ko.agg(set)
    for m, need in C.SIGNATURE_KOS.items():
        if m in enc.index:
            enc.loc[m] = enc.loc[m] & pd.Series({g: need <= kos.get(g, set()) for g in enc.columns})
    return cov, enc, phylum, aln


def mag_profiles():
    """Composition and function per sample from the shared reference (pipeline/09)."""
    cov, enc, phylum, aln = shared_reference()
    info = pd.DataFrame([dict(sample=s, medium=s.split("_")[0], passage=s.split("_")[1], rep=s.split("_")[2])
                         for s in cov.columns]).set_index("sample")
    rel = cov / cov.sum()
    comp = {m: rel[info.index[info.medium == m]] for m in ("V", "NV", "NVNN")}
    f = enc.astype(float).fillna(0) @ cov             # module x sample: coverage of encoding genomes
    func = f / f.sum()
    return comp, func, info, enc, cov, phylum


# Functions whose carriers could change between media (same rows as Fig 4F)
PARTITION_MODULES = [("M00175", "Nitrogen fixation"), ("M00615", "Nitrate assimilation"),
                     ("M00531", "Assimilatory nitrate reduction"), ("M00161", "Photosystem II"),
                     ("M00165", "Calvin cycle"), ("M00123", "Biotin"), ("M00122", "Cobalamin"),
                     ("M00899", "Thiamine salvage"), ("M00116", "Menaquinone")]


def shared_partition():
    """Per sample: cyanobacterial share of encoding-genome coverage for each module,
    plus the cyanobacterial share of mapped coverage and of all reads (Fig 4E)."""
    cov, enc, phylum, aln = shared_reference()
    cyano = phylum.eq(C.CYANO_PHYLUM)
    det = cov >= C.MAG_DETECT
    rows = []
    for m, lab in PARTITION_MODULES:
        e = enc.loc[m] if m in enc.index else pd.Series(False, index=cov.index)
        for s in cov.columns:
            carriers = e & det[s]
            c, h = cov.loc[carriers & cyano, s].sum(), cov.loc[carriers & ~cyano, s].sum()
            rows.append(dict(module=m, function=lab, sample=s, n_cyano=int((carriers & cyano).sum()),
                             n_other=int((carriers & ~cyano).sum()),
                             cyano_percent=100 * c / (c + h) if c + h > 0 else np.nan))
    part = pd.DataFrame(rows)
    share = pd.DataFrame({"cyano_of_mapped": 100 * cov[cyano].sum() / cov.sum(), "reads_aligned": 100 * aln})
    share["cyano_of_all_reads"] = share.cyano_of_mapped * share.reads_aligned / 100
    return part, share


def partition_test(min_samples: int = 3):
    """Does nitrate removal change who carries each function? (Fig 4E)

    Per P10 replicate and module, the cyanobacterial share of the coverage of the
    genomes encoding it (as Fig 3E). NVNN (n = 4) against nitrate-replete V + NV
    (n = 9), exact two-sided Mann-Whitney, BH across the modules with an encoder
    in at least `min_samples` replicates of each group. Also per replicate, the
    cyanobacterial share of all encoded module coverage (community level).
    """
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests
    cov, enc, phylum, _ = shared_reference()
    nvnn, rep = C.MEDIA_NVNN, C.MEDIA_V + C.MEDIA_NV
    samples = rep + nvnn
    cy = phylum.eq(C.CYANO_PHYLUM).values
    x = cov[samples].where(cov[samples] >= C.MAG_DETECT, 0)
    E = enc.fillna(False).astype(float)
    c_cov = E.values[:, cy] @ x.values[cy]          # module x sample, cyanobacterial carriers
    h_cov = E.values[:, ~cy] @ x.values[~cy]
    share = pd.DataFrame(100 * c_cov / np.where(c_cov + h_cov > 0, c_cov + h_cov, np.nan),
                         index=enc.index, columns=samples)
    rows = []
    for m, r in share.iterrows():
        a, b = r[nvnn].dropna(), r[rep].dropna()
        if len(a) < min_samples or len(b) < min_samples:
            continue
        p = mannwhitneyu(a, b, alternative="two-sided", method="exact").pvalue if (a.nunique() + b.nunique() > 1) else 1.0
        rows.append(dict(module=m, nvnn_mean=a.mean(), replete_mean=b.mean(), diff=a.mean() - b.mean(),
                         n_nvnn=len(a), n_replete=len(b), p=p))
    t = pd.DataFrame(rows).set_index("module")
    t["q"] = multipletests(t.p, method="fdr_bh")[1]
    comm = pd.Series(100 * c_cov.sum(0) / (c_cov + h_cov).sum(0), index=samples)
    pc = mannwhitneyu(comm[nvnn], comm[rep], alternative="two-sided", method="exact").pvalue
    return share, t, comm, pc


def run_partition_test():
    share, t, comm, pc = partition_test()
    share.to_csv(C.OUT / "media_partition_P10.csv")
    t.to_csv(C.OUT / "media_partition_test.csv")
    print(f"modules tested: {len(t)}; p < 0.05: {(t.p < 0.05).sum()}; q < 0.05: {(t.q < 0.05).sum()}")
    print(f"community level, cyanobacterial share of encoded function: NVNN {comm[C.MEDIA_NVNN].mean():.1f}% "
          f"vs replete {comm[C.MEDIA_V + C.MEDIA_NV].mean():.1f}%, exact MWU p = {pc:.4f}")
    print(comm.round(1).to_string())
    print(t.sort_values("p").head(25).round(4).to_string())


# ============================================================================
# Coarse-grained partitioning of function between cyanobacteria and heterotrophs
# ----------------------------------------------------------------------------
# Modules are grouped into functional categories (KEGG category, with photosynthesis,
# carbon fixation, nitrogen and sulfur metabolism taken from the subcategory). Per P10
# culture and category: the cyanobacterial share of the summed coverage of the genomes
# that encode the category's modules (shared genome reference, genomes detected at
# MAG_DETECT, encoding at ENCODING with the SIGNATURE_KOS rule). Tests: NVNN against
# V + NV and V against NV, exact two-sided Mann-Whitney, BH across categories.
#     outputs/media_category_partition.csv, outputs/media_category_partition_test.csv
# ============================================================================
CATEGORY_ORDER = ["Photosynthesis", "Carbon fixation", "Nitrogen metabolism", "Sulfur metabolism",
                  "Cofactors and vitamins", "Amino acids", "Carbohydrates", "Lipids", "Nucleotides"]


def module_categories() -> pd.Series:
    t = pd.read_csv(C.DATA / "mag_module_pathwise_completeness.tsv", sep="\t",
                    usecols=["module", "module_category", "module_subcategory"]).set_index("module")
    sub = {"Photosynthesis": "Photosynthesis", "Carbon fixation": "Carbon fixation",
           "Nitrogen metabolism": "Nitrogen metabolism", "Sulfur metabolism": "Sulfur metabolism"}
    cat = {"Metabolism of cofactors and vitamins": "Cofactors and vitamins",
           "Amino acid metabolism": "Amino acids", "Carbohydrate metabolism": "Carbohydrates",
           "Lipid metabolism": "Lipids", "Nucleotide metabolism": "Nucleotides"}
    return pd.Series({m: sub.get(r.module_subcategory, cat.get(r.module_category))
                      for m, r in t.iterrows()}).dropna()


def category_partition():
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests
    cov, enc, phylum, _ = shared_reference()
    groups = {"V": C.MEDIA_V, "NV": C.MEDIA_NV, "NVNN": C.MEDIA_NVNN}
    samples = sum(groups.values(), [])
    cy = phylum.eq(C.CYANO_PHYLUM).values
    x = cov[samples].where(cov[samples] >= C.MAG_DETECT, 0)
    E = enc.fillna(False).astype(float)
    c_cov = pd.DataFrame(E.values[:, cy] @ x.values[cy], index=enc.index, columns=samples)
    h_cov = pd.DataFrame(E.values[:, ~cy] @ x.values[~cy], index=enc.index, columns=samples)
    cats = module_categories().reindex(enc.index).dropna()
    c_cat, h_cat = c_cov.loc[cats.index].groupby(cats).sum(), h_cov.loc[cats.index].groupby(cats).sum()
    share = (100 * c_cat / (c_cat + h_cat).replace(0, np.nan)).reindex(CATEGORY_ORDER)
    long = share.stack().rename("cyano_percent").reset_index().rename(columns={"level_0": "category", "level_1": "sample"})
    long["medium"] = long["sample"].str.split("_").str[0]
    rows = []
    for cat_, r in share.iterrows():
        a, b = r[groups["NVNN"]].dropna(), r[groups["V"] + groups["NV"]].dropna()
        v, nv = r[groups["V"]].dropna(), r[groups["NV"]].dropna()
        rows.append(dict(category=cat_, n_modules=int((cats == cat_).sum()),
                         V=v.mean(), NV=nv.mean(), NVNN=a.mean(),
                         p_nitrate=mannwhitneyu(a, b, method="exact").pvalue,
                         p_vitamin=mannwhitneyu(v, nv, method="exact").pvalue))
    t = pd.DataFrame(rows).set_index("category")
    t["q_nitrate"] = multipletests(t.p_nitrate, method="fdr_bh")[1]
    t["q_vitamin"] = multipletests(t.p_vitamin, method="fdr_bh")[1]
    long.to_csv(C.OUT / "media_category_partition.csv", index=False)
    t.to_csv(C.OUT / "media_category_partition_test.csv")
    return long, t


def category_partition_relative(exclude=("V_P10_R3", "NV_P10_R2", "NV_P10_R4")):
    """As category_partition, against each culture's overall cyanobacterial share of detected
    genome coverage: the excess (percentage points) separates a change in who carries a
    category from a change in how abundant the cyanobacteria are. Cultures with no
    cyanobacteria (exclude) are left out: their primary producer has no genome in the
    reference, so every category is 0% cyanobacterial."""
    from scipy.stats import mannwhitneyu
    cov, enc, phylum, _ = shared_reference()
    long, _ = category_partition()
    long = long[~long["sample"].isin(exclude)].copy()
    smp = long["sample"].unique().tolist()
    x = cov[smp].where(cov[smp] >= C.MAG_DETECT, 0)
    comm = 100 * x[phylum.eq(C.CYANO_PHYLUM)].sum() / x.sum()
    long["community_cyano"] = long["sample"].map(comm)
    long["excess"] = long.cyano_percent - long.community_cyano
    nvnn = long.medium == "NVNN"
    rows = []
    for cat_, g in long.groupby("category", sort=False):
        a, b = g[g.medium == "NVNN"], g[g.medium != "NVNN"]
        rows.append(dict(category=cat_, replete_share=b.cyano_percent.mean(), nvnn_share=a.cyano_percent.mean(),
                         p_share=mannwhitneyu(a.cyano_percent, b.cyano_percent, method="exact").pvalue,
                         replete_excess=b.excess.mean(), nvnn_excess=a.excess.mean(),
                         p_excess=mannwhitneyu(a.excess, b.excess, method="exact").pvalue))
    t = pd.DataFrame(rows).set_index("category").reindex(CATEGORY_ORDER)
    long.to_csv(C.OUT / "media_category_partition_relative.csv", index=False)
    t.to_csv(C.OUT / "media_category_partition_relative_test.csv")
    return long, t, comm


def category_bias(exclude=("V_P10_R3", "NV_P10_R2", "NV_P10_R4"), eps=0.01, passage="P10"):
    """Per culture and category, independent of how abundant each group is: the mean
    fraction of cyanobacterial genome coverage that encodes the category's modules,
    against the same fraction for the other bacteria, as log2((fc + eps) / (fh + eps)).
    Positive: per unit of genome coverage, cyanobacteria encode the category more often
    than heterotrophs. NVNN against V + NV, exact Mann-Whitney. Cultures without
    cyanobacteria are excluded (no cyanobacterial coverage to normalise by)."""
    from scipy.stats import mannwhitneyu
    cov, enc, phylum, _ = shared_reference()
    if passage == "P10":
        nvnn = C.MEDIA_NVNN
        rep = [s_ for s_ in C.MEDIA_V + C.MEDIA_NV if s_ not in exclude]
    else:   # every replicate lineage at that passage (all carry cyanobacteria at P4)
        nvnn = [c for c in cov.columns if c.startswith(f"NVNN_{passage}_R")]
        rep = [c for c in cov.columns if c.startswith((f"V_{passage}_R", f"NV_{passage}_R")) and c not in exclude]
    smp = rep + nvnn
    cy = phylum.eq(C.CYANO_PHYLUM).values
    x = cov[smp].where(cov[smp] >= C.MAG_DETECT, 0)
    E = enc.fillna(False).astype(float)
    cats = module_categories().reindex(enc.index).dropna()
    tot_c, tot_h = x.values[cy].sum(0), x.values[~cy].sum(0)
    rows, stats = [], []
    for cat_ in CATEGORY_ORDER:
        e = E.loc[cats.index[cats == cat_]]
        fc = (e.values[:, cy] @ x.values[cy]).mean(0) / tot_c
        fh = (e.values[:, ~cy] @ x.values[~cy]).mean(0) / tot_h
        r = np.log2((fc + eps) / (fh + eps))
        for s_, a_, b_, c_ in zip(smp, fc, fh, r):
            rows.append(dict(category=cat_, sample=s_, medium=s_.split("_")[0], frac_cyano=a_, frac_other=b_, log2_bias=c_))
        rr = pd.Series(r, index=smp)
        stats.append(dict(category=cat_, replete=rr[rep].mean(), nvnn=rr[nvnn].mean(),
                          frac_cyano_replete=fc[:len(rep)].mean(), frac_cyano_nvnn=fc[len(rep):].mean(),
                          frac_other_replete=fh[:len(rep)].mean(), frac_other_nvnn=fh[len(rep):].mean(),
                          p=mannwhitneyu(rr[nvnn], rr[rep], method="exact").pvalue))
    long, t = pd.DataFrame(rows), pd.DataFrame(stats).set_index("category")
    long.to_csv(C.OUT / f"media_category_bias_{passage}.csv", index=False)
    t.to_csv(C.OUT / f"media_category_bias_test_{passage}.csv")
    return long, t


def run_category_partition():
    long, t = category_partition()
    print(t.round(3).to_string())
    l2, t2, comm = category_partition_relative()
    print(t2.round(3).to_string())
    return long, t


def sample_order(samples):
    key = {"P0": 0, "P4": 1, "P5": 2, "P10": 3}
    return sorted(samples, key=lambda s: (["V", "NV", "NVNN"].index(s.split("_")[0]), key[s.split("_")[1]], s))


def run_shared_partition():
    import matplotlib as mpl
    part, share = shared_partition()
    order = sample_order(share.index)
    piv = part.pivot_table(index="module", columns="sample", values="cyano_percent", dropna=False)
    piv = piv.reindex(index=[m for m, _ in PARTITION_MODULES], columns=order)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 5.6), sharex=True,
                                 gridspec_kw=dict(height_ratios=[1, 2.4], hspace=0.08))
    x = np.arange(len(order))
    cols = [MED_COL[s.split("_")[0]] for s in order]
    a1.bar(x, share.loc[order, "cyano_of_mapped"], color=cols, width=0.8, label="share of mapped coverage")
    a1.scatter(x, share.loc[order, "reads_aligned"], marker="_", s=60, color="#222222", label="reads aligned to the reference")
    a1.set_ylim(0, 100); a1.set_ylabel("Cyanobacteria (%)", fontsize=8)
    a1.legend(fontsize=6.5, frameon=False, loc="upper right", ncol=2)
    cmap = mpl.colormaps["BuGn"].copy(); cmap.set_bad("#eeeeee")
    im = a2.pcolormesh(np.arange(len(order) + 1) - 0.5, np.arange(len(piv) + 1) - 0.5,
                       np.ma.masked_invalid(piv.values), cmap=cmap, vmin=0, vmax=100,
                       edgecolors="white", linewidth=0.4)
    a2.set_ylim(len(piv) - 0.5, -0.5)
    a2.set_yticks(range(len(PARTITION_MODULES)), [l for _, l in PARTITION_MODULES], fontsize=7)
    a2.set_xticks(x, [s.replace("_", " ") for s in order], rotation=90, fontsize=6)
    for b in [i for i in range(1, len(order)) if order[i].split("_")[0] != order[i - 1].split("_")[0]]:
        for a in (a1, a2):
            a.axvline(b - 0.5, color="#555555", lw=0.8)
    cb = fig.colorbar(im, ax=a2, fraction=0.02, pad=0.01)
    cb.set_label("cyanobacterial share of\nencoding-genome coverage (%)", fontsize=7)
    a2.text(0, -0.42, "grey: no detected genome encodes the module", transform=a2.transAxes, fontsize=6, color="#555555")
    for a in (a1, a2):
        a.spines[["top", "right"]].set_visible(False)
    a1.set_title("Who carries each function, shared genome reference (242 genomes, all media samples)",
                 fontsize=9, fontweight="bold")
    fig.savefig(C.OUT / "media_partition_shared.pdf", bbox_inches="tight")
    fig.savefig(C.OUT / "media_partition_shared.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    part.to_csv(C.OUT / "media_partition_shared.csv", index=False)
    share.loc[order].to_csv(C.OUT / "media_cyano_share_shared.csv")
    print(share.loc[order].round(1).to_string())
    p10 = part[part["sample"].str.contains("_P10_")].assign(medium=lambda d: d["sample"].str.split("_").str[0])
    print(p10.pivot_table(index="function", columns="medium", values="cyano_percent").round(0).to_string())


def run_functional_trajectories():
    from scipy.spatial.distance import braycurtis
    comp, func, info, enc, cov, phylum = mag_profiles()
    info = info.drop(index=["NV_P10_R3"], errors="ignore")   # excluded throughout, see common.MEDIA_NV
    rows = []
    for medium in ("V", "NV", "NVNN"):
        par = f"{medium}_P0_parent"
        reps = info[(info.medium == medium) & info.rep.str.startswith("R")]
        for smp, r in reps.iterrows():
            rows.append(dict(medium=medium, passage=r.passage, rep=r.rep,
                             composition=braycurtis(comp[medium][par], comp[medium][smp]),
                             function=braycurtis(func[par], func[smp])))
    to_parent = pd.DataFrame(rows)
    spread = []
    for (medium, pas), d in info[info.rep.str.startswith("R")].groupby(["medium", "passage"]):
        ix = list(d.index)
        for i, a in enumerate(ix):
            for b in ix[i + 1:]:
                spread.append(dict(medium=medium, passage=pas, composition=braycurtis(comp[medium][a], comp[medium][b]),
                                   function=braycurtis(func[a], func[b])))
    spread = pd.DataFrame(spread)
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.7))
    for ax, tab, ttl, yl in ((axes[0], to_parent, "Distance from the parent culture", "Bray-Curtis dissimilarity\nto the P0 parent"),
                             (axes[1], spread, "Dissimilarity among replicates", "Bray-Curtis dissimilarity\nbetween replicates")):
        for k, medium in enumerate(("V", "NV", "NVNN")):
            off = (k - 1) * 0.25                       # dodge the media so the bars do not overlap
            for kind, ls in (("composition", "--"), ("function", "-")):
                d = tab[tab.medium == medium]
                m = d.groupby("passage")[kind].agg(["mean", "std"])
                m = m.assign(x=m.index.map(PASSAGE_X)).sort_values("x")
                pre = ([0], [0.0]) if tab is to_parent else ([], [])
                ax.plot(pre[0] + list(m.x + off), pre[1] + list(m["mean"]), ls=ls, color=MED_COL[medium], lw=1.4)
                ax.errorbar(m.x + off, m["mean"], yerr=m["std"], ls="none", color=MED_COL[medium],
                            marker="o", ms=3, capsize=1.5, lw=0.8)
        ax.set_xticks([0, 4, 10], ["P0", "P4", "P10"], fontsize=7); ax.set_ylim(0, 1)
        ax.tick_params(labelsize=7)
        ax.set_ylabel(yl, fontsize=7); ax.set_title(ttl, fontsize=7.5, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
    for a, l in zip(axes, "AB"):
        a.text(-0.28, 1.04, l, transform=a.transAxes, fontsize=11, fontweight="bold")
    hs = [plt.Line2D([], [], color="black", ls="--", label="composition (genomes)"),
          plt.Line2D([], [], color="black", ls="-", label="function (modules)")]
    hs += [plt.Line2D([], [], color=c, lw=1.6, label=g) for g, c in MED_COL.items()]
    fig.legend(handles=hs, fontsize=6.5, frameon=False, loc="lower center", ncol=5, handlelength=2.2)
    fig.subplots_adjust(left=0.12, right=0.98, top=0.88, bottom=0.25, wspace=0.45)
    C.polish(fig)
    fig.savefig(C.OUT / "media_functional_trajectories.pdf")
    fig.savefig(C.OUT / "media_functional_trajectories.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    to_parent.to_csv(C.OUT / "media_functional_trajectory_to_parent.csv", index=False)
    summ = pd.concat({"to_parent": to_parent.groupby(["medium", "passage"])[["composition", "function"]].mean(),
                      "spread": spread.groupby(["medium", "passage"])[["composition", "function"]].mean()})
    summ.to_csv(C.OUT / "media_functional_trajectory_summary.csv")
    print(summ.round(3).to_string())
    return summ


def nitrogen_module_partition(exclude=("V_P10_R3", "NV_P10_R2", "NV_P10_R4")):
    """Module-level version of the nitrogen row of category_bias (Results, media section),
    which replaced the category figure: the category's shift towards cyanobacteria in NVNN
    is the nitrogenase of the heterocystous diazotrophs alone. Per culture, the percentage of
    cyanobacterial and of heterotrophic genome coverage that encodes nitrogen fixation, the
    other nitrogen modules together, and denitrification; NVNN against V + NV, exact
    Mann-Whitney, at P4 (all lineages with cyanobacteria) and P10.
        outputs/media_nitrogen_modules.csv"""
    from scipy.stats import mannwhitneyu
    cov, enc, phylum, _ = shared_reference()
    cats = module_categories().reindex(enc.index).dropna()
    nmods = list(cats.index[cats == "Nitrogen metabolism"])
    sets = {"nitrogen fixation (M00175)": ["M00175"],
            "other nitrogen modules": [m for m in nmods if m != "M00175"],
            "denitrification (M00529)": ["M00529"]}
    cy = phylum.eq(C.CYANO_PHYLUM).values
    E = enc.fillna(False).astype(float)
    rows = []
    for passage in ("P4", "P10"):
        if passage == "P10":
            nvnn = C.MEDIA_NVNN
            rep = [s_ for s_ in C.MEDIA_V + C.MEDIA_NV if s_ not in exclude]
        else:
            nvnn = [c for c in cov.columns if c.startswith("NVNN_P4_R")]
            rep = [c for c in cov.columns if c.startswith(("V_P4_R", "NV_P4_R"))]
        smp = rep + nvnn
        x = cov[smp].where(cov[smp] >= C.MAG_DETECT, 0)
        for name, mods in sets.items():
            e = E.loc[mods].values
            for group, mask in (("cyanobacteria", cy), ("heterotrophs", ~cy)):
                f = pd.Series(100 * (e[:, mask] @ x.values[mask]).mean(0) / x.values[mask].sum(0), index=smp)
                rows.append(dict(passage=passage, modules=name, carriers=group,
                                 replete_mean=f[rep].mean(), replete_min=f[rep].min(), replete_max=f[rep].max(),
                                 nvnn_mean=f[nvnn].mean(), nvnn_min=f[nvnn].min(), nvnn_max=f[nvnn].max(),
                                 p=mannwhitneyu(f[nvnn], f[rep], method="exact").pvalue))
    t = pd.DataFrame(rows)
    t.to_csv(C.OUT / "media_nitrogen_modules.csv", index=False)
    return t


def main() -> None:
    run_media()
    run_trajectories()
    run_shared_partition()
    run_functional_trajectories()
    print(nitrogen_module_partition().round(3).to_string(index=False))


if __name__ == "__main__":
    main()
