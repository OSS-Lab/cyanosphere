"""
Taxonomic composition of the enrichments and the reservoir: diversity (Fig 1),
composition and PCA (Fig 2), UniFrac and the phylogenetic core, divergence
between replicates, and the core-tree figure (SI).

Run:  python taxonomy.py      (from analysis/, after the modules before it in
      the README order)
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import mannwhitneyu
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import common as C
import itertools
from skbio import TreeNode
from skbio.diversity import beta_diversity
from skbio.stats.distance import permanova
from skbio.stats.ordination import pcoa
from matplotlib.patches import Polygon
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
from pathlib import Path
from matplotlib.lines import Line2D

# ============================================================================
# Diversity and taxonomic composition: run_diversity_and_taxonomy()
# ----------------------------------------------------------------------------
# Figures 1 and 2: how the enrichment end-points differ from the reservoir in
# diversity (Fig 1) and in taxonomic composition (Fig 2).
#
# Both figures use EcoPhylo profiles of a single-copy ribosomal protein. The unit
# of counting is a distinct ribosomal-protein sequence cluster (97% nucleotide
# identity), called a population. RpL19 is used here; function.py repeats the
# functional analysis with RpS16.
#
# Two counts exist in the source data and are kept apart:
#   clusters  distinct 97% RpL19 sequence clusters (659 across all samples)
#   species   those clusters collapsed by GTDB species assignment (314)
# Richness, evenness and the rank-abundance curves are all computed on clusters,
# the quantity "population" refers to.
#
# Produces
#   outputs/diversity_metrics.csv
#   outputs/fig1.{pdf,svg}
#   outputs/fig2.{pdf,svg}
#
# Run:  python taxonomy.py
# ============================================================================

ORDER = C.SOURCE + C.ENRICHMENT


def load_clusters(relative: bool = False) -> pd.DataFrame:
    """97% RpL19 sequence clusters x samples. `relative` divides each sample by
    its own total, which is what the published Fig 2D ordination used."""
    d = pd.read_csv(C.DATA / "rpl19_cluster_coverage.csv")
    d = d.rename(columns={d.columns[0]: "cluster"}).set_index("cluster")
    d.columns = [c.replace("WAMOS_01", "W_01") for c in d.columns]
    d = d[[s for s in ORDER if s in d.columns]]
    return d / d.sum() if relative else d


def diversity(cov: pd.DataFrame) -> pd.DataFrame:
    rel = cov / cov.sum()
    H = -(rel * np.log(rel.where(rel > 0))).sum()
    m = pd.DataFrame({"shannon": H, "richness": (cov > 0).sum()})
    m["group"] = ["Source" if s in C.SOURCE else "Enrichment" for s in m.index]
    m["top3_percent"] = [100 * cov[s].nlargest(3).sum() / cov[s].sum() for s in m.index]
    return m


def panel_box(ax, m, col, title, ylabel):
    data = [m.loc[m.group == g, col].values for g in ("Source", "Enrichment")]
    bp = ax.boxplot(data, positions=[0, 1], widths=0.55, patch_artist=True,
                    showfliers=False, medianprops=dict(color="black", lw=1.2))
    for patch, g in zip(bp["boxes"], ("Source", "Enrichment")):
        patch.set_facecolor(C.GROUP_COLOR[g])
        patch.set_alpha(0.28)
        patch.set_edgecolor(C.GROUP_COLOR[g])
    rng = np.random.default_rng(0)
    for x, g in enumerate(("Source", "Enrichment")):
        v = m.loc[m.group == g, col].values
        ax.scatter(x + rng.uniform(-0.13, 0.13, len(v)), v, s=22,
                   color=C.GROUP_COLOR[g], edgecolor="white", lw=0.5, zorder=3)
    p = mannwhitneyu(data[0], data[1], alternative="two-sided").pvalue
    y = max(np.max(data[0]), np.max(data[1]))
    ax.plot([0, 0, 1, 1], [y * 1.04, y * 1.08, y * 1.08, y * 1.04], lw=0.8, color="#444444")
    ax.text(0.5, y * 1.10, f"$p$ = {p:.1e}", ha="center", fontsize=7.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Source", "Enrichment"])
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left")
    return p


def panel_rank(ax, cov):
    for g, members in (("Source", C.SOURCE), ("Enrichment", C.ENRICHMENT)):
        curves = []
        for s in members:
            if s not in cov.columns:
                continue
            v = np.sort(cov[s].values)[::-1]
            v = 100 * v / v.sum()
            curves.append(v[:60])
            ax.plot(range(1, len(v[:60]) + 1), v[:60], color=C.GROUP_COLOR[g],
                    alpha=0.22, lw=0.8)
        A = np.vstack(curves)
        ax.plot(range(1, A.shape[1] + 1), A.mean(0), color=C.GROUP_COLOR[g],
                lw=2.0, label=g)
    ax.set_yscale("log")
    ax.set_xlabel("Population rank")
    ax.set_ylabel("Relative coverage (%)")
    ax.set_title("D   Rank abundance", loc="left")
    ax.legend()


def phylum_table() -> pd.DataFrame:
    """Relative RpL19 coverage by phylum: the 97% clusters of load_clusters() summed by
    GTDB phylum, unnamed clusters placed on the RpL19 tree (common.rpl19_taxonomy).
    Replaces an earlier table from per-sample anvi-estimate-scg-taxonomy, which dropped
    unclassified RpL19 genes and so counted different populations from Fig 1 and 2A/D."""
    cov = load_clusters()
    tax = C.rpl19_taxonomy().reindex(cov.index)
    d = cov.groupby(tax.t_phylum.fillna("Unknown")).sum()
    return d / d.sum()


def panel_phylum(ax, ph):
    top = ph.mean(axis=1).sort_values(ascending=False)
    keep = top[top >= 0.01].index.tolist()
    other = ph.loc[[i for i in ph.index if i not in keep]].sum()
    m = pd.concat([ph.loc[keep], pd.DataFrame({"Other (<1% mean)": other}).T])
    cmap = plt.get_cmap("tab20")
    bottom = np.zeros(m.shape[1])
    for i, ph_name in enumerate(m.index):
        ax.bar(range(m.shape[1]), m.loc[ph_name].values, bottom=bottom, width=0.82,
               color=cmap(i % 20), edgecolor="white", lw=0.4, label=ph_name)
        bottom += m.loc[ph_name].values
    ax.set_xticks(range(m.shape[1]))
    ax.set_xticklabels(m.columns, rotation=90, fontsize=7)
    ax.set_ylabel("Relative coverage")
    ax.set_ylim(0, 1)
    ax.set_title("B   Phylum composition (stacked by decreasing mean coverage)", loc="left")
    ax.legend(fontsize=6.4, ncol=2, loc="center left", bbox_to_anchor=(1.005, 0.5))
    return m


def panel_pca(ax, mat, title, zscore=False):
    """
    Ordination of a composition matrix.

    zscore=False reproduces the published panels: the PCA is run on the matrix
    as it stands, so a phylum or a cluster with large absolute variance carries
    the first component. That is a deliberate choice for a compositional bar
    chart's companion ordination, but it has to be stated, because z-scoring the
    same matrix gives a much lower PC1 (36% instead of 97% at phylum level) by
    giving every phylum equal weight.
    """
    samples = [s for s in ORDER if s in mat.columns]
    X = mat[samples].T.values.astype(float)
    if zscore:
        X = StandardScaler().fit_transform(X)
    p = PCA(n_components=2).fit(X)
    sc = p.transform(X)
    var = p.explained_variance_ratio_ * 100
    for g, members in (("Source", C.SOURCE), ("Enrichment", C.ENRICHMENT)):
        idx = [i for i, s in enumerate(samples) if s in members]
        ax.scatter(sc[idx, 0], sc[idx, 1], s=64, color=C.GROUP_COLOR[g],
                   edgecolor="white", lw=0.8, label=g, zorder=3)
    F, pv = C.permanova(mat, {"Source": C.SOURCE, "Enrichment": C.ENRICHMENT})
    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    ax.set_title(title, loc="left")
    ax.annotate(f"PERMANOVA $F$ = {F:.1f}, $p$ = {pv:.3f}", xy=(0.03, 0.03),
                xycoords="axes fraction", fontsize=7.5, color="#555555")
    ax.legend(loc="best")
    return var, F, pv


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


def run_diversity_and_taxonomy() -> None:
    C.style()
    cov = load_clusters()
    m = diversity(cov)
    m.to_csv(C.OUT / "diversity_metrics.csv")
    print(m.groupby("group")[["shannon", "richness", "top3_percent"]]
          .agg(["mean", "min", "max"]).round(2).to_string())

    fig = plt.figure(figsize=(11, 7.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0], hspace=0.42, wspace=0.34)
    placeholder(fig.add_subplot(gs[0, :]), "A",
                "Sampling and enrichment strategy (eight independent lineages, >1 year)")
    p_sh = panel_box(fig.add_subplot(gs[1, 0]), m, "shannon", "B   Shannon diversity", "Entropy (H)")
    p_ri = panel_box(fig.add_subplot(gs[1, 1]), m, "richness", "C   Population richness",
                     "RpL19 clusters detected")
    panel_rank(fig.add_subplot(gs[1, 2]), cov)
    C.save(fig, "fig1")
    print(f"Fig 1: Shannon p={p_sh:.2e}, richness p={p_ri:.2e}")

    ph = phylum_table()
    fig = plt.figure(figsize=(11.5, 9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.95], hspace=0.5, wspace=0.3)
    placeholder(fig.add_subplot(gs[0, :]), "A",
                "EcoPhylo abundance heatmap of RpL19 clusters",
                "anvi-interactive, composited in Affinity")
    panel_phylum(fig.add_subplot(gs[1, 0]), ph)
    var, F, pv = panel_pca(fig.add_subplot(gs[1, 1]), load_clusters(relative=True),
                           "D   Ordination of RpL19 cluster composition")
    C.save(fig, "fig2")
    print(f"Fig 2: cluster PCA PC1={var[0]:.1f}% PC2={var[1]:.1f}%, PERMANOVA F={F:.1f} p={pv:.3f}")
    varp, Fp, pvp = panel_pca(plt.figure(figsize=(4, 4)).add_subplot(111), ph,
                              "phylum-level ordination")
    plt.close("all")
    print(f"       phylum PCA PC1={varp[0]:.1f}% PC2={varp[1]:.1f}%, "
          f"PERMANOVA F={Fp:.1f} p={pvp:.3f}")


# ============================================================================
# Phylogenetic beta-diversity and the core community phylogeny: run_phylogenetic_core()
# ----------------------------------------------------------------------------
# Both are computed on the EcoPhylo RpL19 tree and the per-sample coverage of each
# 97% sequence cluster.
#
# 1. UniFrac. Unweighted and weighted-normalised UniFrac between all 16 samples,
#    with PCoA and PERMANOVA. Being phylogenetic, it does not depend on taxonomic
#    names, which matters because much of the tree is unnamed at species level.
#
# 2. The core community phylogeny, branch-based, after Bewick and Camper
#    (Microbiome 2026, doi 10.1186/s40168-026-02520-8). A branch is core to a
#    habitat if the clade below it is detected in at least a threshold fraction of
#    that habitat's samples. Summing core branch length gives core Faith's PD;
#    comparing core branch sets between habitats gives core UniFrac and a
#    phylogenetic Venn partition, i.e. convergence measured as shared evolutionary
#    history rather than as shared taxon names.
#
# Habitats are the reservoir (8 samples) and the year-long enrichments (8
# independent replicates). Outputs go to outputs/ with the prefix `phylo_`.
#
# Reads only data/.
# ============================================================================

# --- parameters ---------------------------------------------------------------
# A cluster counts as detected in a sample above this mean coverage. EcoPhylo
# coverages already count only reads mapping the ribosomal protein, so zero means
# no read recruited; the threshold guards against single stray reads.
DETECTION = 0.1


# Fraction of a habitat's samples in which a branch must be detected to be core.
# 0.5 is the default of Bewick and Camper; the sweep shows the result does not
# depend on it.
CORE_THRESHOLD = 0.5


CORE_SWEEP = [0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0]


HABITAT = {"Source": C.SOURCE, "Enrichment": C.ENRICHMENT}


def load() -> tuple[TreeNode, pd.DataFrame]:
    """Return the midpoint-rooted RpL19 tree and the cluster-by-sample coverage."""
    cov = pd.read_csv(C.DATA / "rpl19_cluster_coverage.csv", encoding="utf-8-sig")
    cov = cov.set_index(cov.columns[0])
    cov.index.name = "cluster"
    cov = cov[C.SOURCE + C.ENRICHMENT]

    tree = TreeNode.read(str(C.DATA / "rpl19_cluster_tree.nwk"))
    # The Newick standard reads an unquoted underscore as a space, so skbio
    # returns "M 132 Ribosomal L19 ..." for tips written "M_132_Ribosomal_L19_...".
    # Put the underscores back so tip names match the coverage table.
    for tip in tree.tips():
        tip.name = tip.name.replace(" ", "_")
    # The EcoPhylo tree is unrooted (FastTree). The core-phylogeny measures need
    # a rooted tree, because "the clade below a branch" needs a direction.
    # Midpoint rooting is the neutral choice and is applied once, here, so both
    # analyses use the same tree.
    tree = tree.root_at_midpoint()
    for node in tree.traverse():
        if node.length is None:
            node.length = 0.0

    tips = {t.name for t in tree.tips()}
    missing = tips.symmetric_difference(set(cov.index))
    if missing:
        raise SystemExit(f"tree and coverage table disagree on {len(missing)} tips")
    return tree, cov.loc[[t.name for t in tree.tips()]]


# --- 1. UniFrac ---------------------------------------------------------------
def unifrac_block(tree: TreeNode, cov: pd.DataFrame) -> pd.DataFrame:
    """UniFrac distances, ordination and PERMANOVA for both weightings."""
    samples = list(cov.columns)
    groups = ["Source" if s in C.SOURCE else "Enrichment" for s in samples]
    counts = cov.T.values
    taxa = list(cov.index)

    rows = []
    for metric, kwargs in [
        ("unweighted_unifrac", {}),
        ("weighted_unifrac", {"normalized": True}),
    ]:
        dm = beta_diversity(metric, counts, ids=samples, taxa=taxa, tree=tree, **kwargs)
        res = permanova(dm, grouping=groups, permutations=9999, seed=0)
        ord_ = pcoa(dm, number_of_dimensions=3)
        pv = ord_.proportion_explained

        dm.to_data_frame().to_csv(C.OUT / f"phylo_{metric}_distances.csv")
        coords = ord_.samples.iloc[:, :2].copy()
        coords.columns = ["PCo1", "PCo2"]
        coords["group"] = groups
        coords.to_csv(C.OUT / f"phylo_{metric}_pcoa.csv")

        # within-group distance: how tightly the eight independent enrichments
        # cluster relative to the eight reservoir samples
        d = dm.to_data_frame()
        within = {
            g: float(np.mean([d.loc[a, b] for a, b in itertools.combinations(m, 2)]))
            for g, m in HABITAT.items()
        }
        between = float(
            np.mean([d.loc[a, b] for a in HABITAT["Source"] for b in HABITAT["Enrichment"]])
        )

        rows.append(
            {
                "metric": metric,
                "PCo1_pct": round(100 * pv.iloc[0], 1),
                "PCo2_pct": round(100 * pv.iloc[1], 1),
                "permanova_F": round(float(res["test statistic"]), 3),
                "permanova_p": float(res["p-value"]),
                "mean_within_source": round(within["Source"], 4),
                "mean_within_enrichment": round(within["Enrichment"], 4),
                "mean_between": round(between, 4),
            }
        )
    return pd.DataFrame(rows)


# --- 2. core community phylogeny ----------------------------------------------
def branch_presence(tree: TreeNode, cov: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
    """Per-branch detection matrix and branch lengths.

    A branch is detected in a sample when at least one tip of the clade below it
    is detected there. This is the branch-based core of Bewick and Camper, which
    they prefer to the tip-based version because it retains the deeper shared
    lineages when strains differ between samples.
    """
    detected = (cov > DETECTION).values  # tips x samples, tips in tree order
    tip_index = {name: i for i, name in enumerate(cov.index)}

    names, lengths, pres = [], [], []
    for node in tree.postorder(include_self=False):
        if node.is_tip():
            mask = detected[tip_index[node.name]]
        else:
            mask = np.zeros(cov.shape[1], dtype=bool)
            for child in node.children:
                mask |= child._mask
        node._mask = mask
        names.append(node.name or f"internal_{len(names)}")
        lengths.append(float(node.length))
        pres.append(mask)

    presence = pd.DataFrame(np.array(pres), index=names, columns=cov.columns)
    return presence, np.array(lengths)


def core_mask(presence: pd.DataFrame, samples: list[str], threshold: float) -> np.ndarray:
    """Branches detected in at least `threshold` of a habitat's samples."""
    need = int(np.ceil(threshold * len(samples)))
    return (presence[samples].sum(axis=1) >= need).values


def core_block(presence: pd.DataFrame, lengths: np.ndarray):
    """Core Faith's PD, core UniFrac and the phylogenetic Venn, over a sweep."""
    total = lengths.sum()
    sweep, venn = [], []

    for tau in CORE_SWEEP:
        masks = {g: core_mask(presence, s, tau) for g, s in HABITAT.items()}
        pds = {g: lengths[m].sum() for g, m in masks.items()}

        src, enr = masks["Source"], masks["Enrichment"]
        shared = lengths[src & enr].sum()
        union = lengths[src | enr].sum()
        # unweighted core UniFrac: branch length unique to one core, over the
        # branch length core to either
        core_uf = (union - shared) / union if union > 0 else np.nan

        sweep.append(
            {
                "threshold": tau,
                "core_PD_source": round(pds["Source"], 3),
                "core_PD_enrichment": round(pds["Enrichment"], 3),
                "core_PD_source_pct_of_tree": round(100 * pds["Source"] / total, 1),
                "core_PD_enrichment_pct_of_tree": round(100 * pds["Enrichment"] / total, 1),
                "core_unifrac": round(core_uf, 4),
                "shared_PD": round(shared, 3),
                "shared_pct_of_union": round(100 * shared / union, 1) if union > 0 else np.nan,
                "n_core_branches_source": int(src.sum()),
                "n_core_branches_enrichment": int(enr.sum()),
            }
        )
        if tau == CORE_THRESHOLD:
            venn = [
                {"compartment": "source core only", "branch_length": round(lengths[src & ~enr].sum(), 3)},
                {"compartment": "shared by both cores", "branch_length": round(shared, 3)},
                {"compartment": "enrichment core only", "branch_length": round(lengths[enr & ~src].sum(), 3)},
                {"compartment": "core to neither", "branch_length": round(total - union, 3)},
                {"compartment": "whole tree", "branch_length": round(total, 3)},
            ]
    return pd.DataFrame(sweep), pd.DataFrame(venn)


def faiths_pd_per_sample(presence: pd.DataFrame, lengths: np.ndarray, cols) -> pd.Series:
    """Faith's PD of each individual sample, for the alpha-diversity comparison."""
    return pd.Series({s: lengths[presence[s].values].sum() for s in cols}, name="faiths_pd")


def run_phylogenetic_core() -> None:
    tree, cov = load()
    print(f"tree: {cov.shape[0]} tips, coverage: {cov.shape[1]} samples")

    uf = unifrac_block(tree, cov)
    uf.to_csv(C.OUT / "phylo_unifrac_summary.csv", index=False)
    print("\n--- UniFrac ---")
    print(uf.to_string(index=False))

    presence, lengths = branch_presence(tree, cov)
    print(f"\nbranches: {len(lengths)}, total branch length {lengths.sum():.2f}")

    sweep, venn = core_block(presence, lengths)
    sweep.to_csv(C.OUT / "phylo_core_sweep.csv", index=False)
    venn.to_csv(C.OUT / "phylo_core_venn.csv", index=False)
    print("\n--- core community phylogeny, threshold sweep ---")
    print(sweep.to_string(index=False))
    print(f"\n--- phylogenetic Venn at threshold {CORE_THRESHOLD} ---")
    print(venn.to_string(index=False))

    fpd = faiths_pd_per_sample(presence, lengths, list(cov.columns)).to_frame()
    fpd["group"] = ["Source" if s in C.SOURCE else "Enrichment" for s in fpd.index]
    fpd.to_csv(C.OUT / "phylo_faiths_pd.csv")
    s, e = fpd.loc[C.SOURCE, "faiths_pd"], fpd.loc[C.ENRICHMENT, "faiths_pd"]
    u, p = mannwhitneyu(s, e, alternative="two-sided")
    print("\n--- Faith's PD per sample ---")
    print(f"source     mean {s.mean():.2f} (range {s.min():.2f}-{s.max():.2f})")
    print(f"enrichment mean {e.mean():.2f} (range {e.min():.2f}-{e.max():.2f})")
    print(f"two-sided Mann-Whitney U = {u:.1f}, p = {p:.2e}")


# ============================================================================
# Divergence between replicates: run_divergence()
# ----------------------------------------------------------------------------
# How far the eight independent enrichments diverged from one another, in
# phylogeny and in function, against the same reservoir baseline.
#
# The enrichments are more variable than the reservoir samples on every measure,
# phylogenetic and functional. Relative to the same baseline, phylogenetic
# composition diverged about half again as much as functional composition, so
# functional divergence is buffered relative to taxonomic divergence.
#
# The absolute values are confounded by the design: the reservoir samples come from
# one lake community, while the enrichments are independent year-long lineages. The
# comparison between metrics is not, because every metric sees the same two groups;
# that ratio is what panel C shows.
#
# Distances:
#   phylogeny  unweighted and weighted UniFrac on the RpL19 EcoPhylo tree
#   function   Bray-Curtis on per-sample-normalised coverage-based and
#              count-based KEGG module profiles
#
# Needs the outputs of run_phylogenetic_core(). Writes
# outputs/fig_divergence.{pdf,svg} and outputs/divergence_summary.csv.
# ============================================================================

COLS = C.SOURCE + C.ENRICHMENT


RNG = np.random.default_rng(0)


N_BOOT = 5000


LABELS = {
    "unweighted_unifrac": "Phylogeny\nunweighted UniFrac",
    "weighted_unifrac": "Phylogeny\nweighted UniFrac",
    "coverage": "Function\ncoverage-based",
    "count": "Function\ncount-based",
}


KIND_C = {"unweighted_unifrac": "#3B4A6B", "weighted_unifrac": "#5A6E96",
          "coverage": "#8D6E63", "count": "#B0301F"}


def distances() -> dict[str, pd.DataFrame]:
    d = {}
    for m in ("unweighted_unifrac", "weighted_unifrac"):
        d[m] = pd.read_csv(C.OUT / f"phylo_{m}_distances.csv", index_col=0).loc[COLS, COLS]
    cov, cnt, _ = C.metric_matrices("enrichments")
    for key, mat in (("coverage", cov), ("count", cnt)):
        x = C.compositional(mat[COLS]).T
        d[key] = pd.DataFrame(squareform(pdist(x.values, metric="braycurtis")),
                              index=COLS, columns=COLS)
    return d


def within(d: pd.DataFrame, g: list[str]) -> np.ndarray:
    return np.array([d.loc[a, b] for a, b in itertools.combinations(g, 2)])


def boot_ratio(d: pd.DataFrame) -> tuple[float, float, float]:
    """Ratio of mean within-enrichment to mean within-source distance, with a
    percentile CI bootstrapped over samples rather than over pairs, because the
    pairs are not independent."""
    obs = within(d, C.ENRICHMENT).mean() / within(d, C.SOURCE).mean()
    vals = []
    for _ in range(N_BOOT):
        s = list(RNG.choice(C.SOURCE, len(C.SOURCE), replace=True))
        e = list(RNG.choice(C.ENRICHMENT, len(C.ENRICHMENT), replace=True))
        ws, we = within(d, s), within(d, e)
        if len(ws) and len(we) and ws.mean() > 0:
            vals.append(we.mean() / ws.mean())
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return obs, lo, hi


def hull(ax, xy, color):
    if len(xy) < 3:
        return
    p = xy[np.argsort(np.arctan2(xy[:, 1] - xy[:, 1].mean(), xy[:, 0] - xy[:, 0].mean()))]
    ax.add_patch(Polygon(p, closed=True, facecolor=color, alpha=0.12,
                         edgecolor=color, lw=1.0, zorder=1))


def ordinate(ax, d: pd.DataFrame, title: str):
    ord_ = pcoa(DistanceMatrix(d.values, ids=list(d.index)), number_of_dimensions=2)
    pv = ord_.proportion_explained
    xy = ord_.samples.values[:, :2]
    for g, members in (("Source", C.SOURCE), ("Enrichment", C.ENRICHMENT)):
        idx = [COLS.index(s) for s in members]
        c = C.GROUP_COLOR[g]
        hull(ax, xy[idx], c)
        ax.scatter(xy[idx, 0], xy[idx, 1], s=46, c=c, edgecolor="white", lw=.8,
                   zorder=3, label=f"{g} (n=8)")
    ax.set_xlabel(f"PCo1 ({100*pv.iloc[0]:.1f}%)")
    ax.set_ylabel(f"PCo2 ({100*pv.iloc[1]:.1f}%)")
    ax.set_title(title, loc="left")
    ax.axhline(0, lw=.5, color="#cccccc", zorder=0)
    ax.axvline(0, lw=.5, color="#cccccc", zorder=0)


def run_divergence() -> None:
    C.style()
    d = distances()

    rows = []
    for k, dm in d.items():
        s, e = within(dm, C.SOURCE), within(dm, C.ENRICHMENT)
        obs, lo, hi = boot_ratio(dm)
        u, p = mannwhitneyu(e, s, alternative="two-sided")
        rows.append(dict(distance=k, mean_within_source=round(s.mean(), 4),
                         mean_within_enrichment=round(e.mean(), 4),
                         ratio=round(obs, 3), ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                         p=p))
    summ = pd.DataFrame(rows)
    summ.to_csv(C.OUT / "divergence_summary.csv", index=False)
    print(summ.to_string(index=False))

    fig = plt.figure(figsize=(11.5, 8.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=.42, wspace=.26)

    ax = fig.add_subplot(gs[0, 0])
    ordinate(ax, d["weighted_unifrac"], "A   Phylogeny: weighted UniFrac")
    ax.legend(loc="best")
    ordinate(fig.add_subplot(gs[0, 1]), d["count"],
             "B   Function: count-based module profile")

    ax = fig.add_subplot(gs[1, :])
    order = ["unweighted_unifrac", "weighted_unifrac", "coverage", "count"]
    for i, k in enumerate(order):
        s, e = within(d[k], C.SOURCE), within(d[k], C.ENRICHMENT)
        for off, vals, col, lab in ((-.16, s, C.GROUP_COLOR["Source"], "Source"),
                                    (.16, e, C.GROUP_COLOR["Enrichment"], "Enrichment")):
            x = i + off + RNG.normal(0, .028, len(vals))
            # each metric is on its own scale, so plot every distance relative to
            # that metric's own within-source mean; 1.0 is then "as variable as
            # the reservoir is" on all four
            ax.scatter(x, vals / s.mean(), s=13, c=col, alpha=.55, lw=0, zorder=2)
            ax.hlines(vals.mean() / s.mean(), i + off - .1, i + off + .1,
                      color=col, lw=2.2, zorder=3)
        r = summ.loc[summ.distance == k].iloc[0]
        ax.text(i, ax.get_ylim()[1] * .02 + 2.62,
                f"{r.ratio:.2f}\n[{r.ci_lo:.2f}, {r.ci_hi:.2f}]",
                ha="center", va="top", fontsize=8.5, color=KIND_C[k], fontweight="bold")
    ax.axhline(1, color="#888888", lw=.9, ls="--", zorder=1)
    ax.set_xticks(range(4))
    ax.set_xticklabels([LABELS[k] for k in order])
    ax.set_ylabel("pairwise distance\n(relative to mean within-reservoir distance)")
    ax.set_title("C   Divergence among independent replicates, on a common scale",
                 loc="left")
    ax.set_xlim(-.55, 3.55)
    ax.text(3.45, 1.03, "reservoir baseline", ha="right", va="bottom",
            fontsize=8, color="#888888")

    C.save(fig, "fig_divergence")


# ============================================================================
# The core RpL19 phylogeny as a figure: run_core_tree_figure()
# ----------------------------------------------------------------------------
# The branch-based core of run_phylogenetic_core(), drawn on the full RpL19 tree.
# Of the 659 clusters, 36 are detected in at least half of the eight enrichments,
# 190 in at least half of the reservoir samples, and one in both. The 36
# enrichment-core clusters are labelled at family level; reservoir-core branches
# are drawn in a lighter colour.
#
# Needs the outputs of run_phylogenetic_core(). Writes outputs/fig_core_tree.{pdf,svg}.
# ============================================================================



TAXONOMY = C.DATA / "rpl19_cluster_taxonomy.tsv"   # EcoPhylo 06_MISC_DATA/Ribosomal_L19_scg_taxonomy_data.tsv


ENR_C = "#2E8B63"      # enrichment core, the subject of the figure


SHARED_C = "#12171D"   # core to both


SRC_C = "#A9B6CC"      # reservoir core, deliberately held back


NONE_C = "#EAEDF1"


STYLE = {"shared": (SHARED_C, 1.30, 5),
         "enrichment": (ENR_C, 1.15, 4),
         "source": (SRC_C, .55, 3),
         "none": (NONE_C, .32, 1)}


PHYLUM_C = {
    "Cyanobacteriota": C.CYANO_C,
    "Pseudomonadota": "#6E86AE",
    "Actinomycetota": "#C2A25A",
    "Bacteroidota": "#B08A7E",
    "Verrucomicrobiota": "#8E9BA6",
    "Planctomycetota": "#7E8C74",
    "Chloroflexota": "#A79BB0",
    "unclassified": "#D8DCE1",
}


OTHER_C = "#C4CAD1"


def run_core_tree_figure() -> None:
    C.style()
    tree, cov = load()
    presence, lengths = branch_presence(tree, cov)
    src = core_mask(presence, C.SOURCE, CORE_THRESHOLD)
    enr = core_mask(presence, C.ENRICHMENT, CORE_THRESHOLD)

    nodes = list(tree.postorder(include_self=False))
    kind, tip_kind = {}, {}
    for i, nd in enumerate(nodes):
        s, e = bool(src[i]), bool(enr[i])
        k = "shared" if (s and e) else "enrichment" if e else "source" if s else "none"
        kind[id(nd)] = k
        if nd.is_tip():
            tip_kind[nd.name] = k

    tips = list(tree.tips())
    for i, t in enumerate(tips):
        t._y = float(i)
    for nd in tree.postorder(include_self=True):
        if not nd.is_tip():
            nd._y = float(np.mean([c._y for c in nd.children]))
    tree._x = 0.0
    for nd in tree.preorder(include_self=False):
        nd._x = nd.parent._x + max(float(nd.length or 0), 0)
    xmax = max(nd._x for nd in nodes)

    tax = pd.read_csv(TAXONOMY, sep="\t")
    phy = dict(zip(tax.split_name, tax.t_phylum))
    fam = dict(zip(tax.split_name, tax.t_family))
    gen = dict(zip(tax.split_name, tax.t_genus))
    names = [t.name for t in tips]
    plab = [phy.get(n) if isinstance(phy.get(n), str) else "unclassified" for n in names]

    fig = plt.figure(figsize=(6.6, 8.5))
    ax = fig.add_axes([0.035, 0.105, 0.485, 0.795])
    ax.axis("off")
    ax.set_ylim(len(tips) + 3, -4)
    ax.set_xlim(-xmax * .012, xmax * 1.02)

    for nd in nodes:
        col, lw, z = STYLE[kind[id(nd)]]
        ax.plot([nd.parent._x, nd._x], [nd._y, nd._y], color=col, lw=lw, zorder=z,
                solid_capstyle="butt")
    for nd in tree.postorder(include_self=True):
        if nd.is_tip() or not nd.children:
            continue
        ys = [c._y for c in nd.children]
        best = max(nd.children, key=lambda c: STYLE[kind[id(c)]][2])
        col, lw, z = STYLE[kind[id(best)]]
        ax.plot([nd._x, nd._x], [min(ys), max(ys)], color=col, lw=lw, zorder=z,
                solid_capstyle="butt")

    step = min([.5, .25, .2, .1], key=lambda v: abs(v - xmax / 5))
    y0 = len(tips) + 1.5
    ax.plot([0, step], [y0, y0], color="#5A6570", lw=1.0)
    ax.text(step / 2, y0 + 2.6, f"{step:g} substitutions per site",
            ha="center", va="top", fontsize=6.8, color="#5A6570")

    # phylum strip, for context only
    sx = fig.add_axes([0.535, 0.105, 0.020, 0.795])
    sx.axis("off"); sx.set_ylim(len(tips) + 3, -4); sx.set_xlim(0, 1)
    for i, lab in enumerate(plab):
        sx.add_patch(plt.Rectangle((0, i - .5), 1, 1, lw=0,
                                   color=PHYLUM_C.get(lab, OTHER_C)))

    # right margin: name the enrichment-core populations
    tx = fig.add_axes([0.560, 0.105, 0.430, 0.795])
    tx.axis("off"); tx.set_ylim(len(tips) + 3, -4); tx.set_xlim(0, 1)

    core_idx = [i for i, n in enumerate(names)
                if tip_kind.get(n) in ("enrichment", "shared")]
    runs = []
    for _, grp in itertools.groupby(enumerate(core_idx), lambda p: p[1] - p[0]):
        idxs = [i for _, i in grp]
        labs = []
        for i in idxs:
            n = names[i]
            f, g = fam.get(n), gen.get(n)
            if isinstance(g, str) and phy.get(n) == C.CYANO_PHYLUM:
                labs.append(g)
            elif isinstance(f, str):
                labs.append(f)
            else:
                labs.append("unclassified")
        top = pd.Series(labs).value_counts()
        txt = top.index[0] + (f" +{len(top) - 1}" if len(top) > 1 else "")
        runs.append([float(np.mean(idxs)), txt, len(idxs),
                     any(phy.get(names[i]) == C.CYANO_PHYLUM for i in idxs)])

    gap = len(tips) * 0.0235
    runs.sort()
    for i in range(1, len(runs)):
        if runs[i][0] - runs[i - 1][0] < gap:
            runs[i][0] = runs[i - 1][0] + gap
    for (y, txt, n, cy), i0 in zip(runs, [r[0] for r in runs]):
        pass
    for r, idxs in zip(runs, [None] * len(runs)):
        y, txt, n, cy = r
        tx.text(0.085, y, f"{txt}" + (f"  ×{n}" if n > 1 else ""),
                va="center", ha="left", fontsize=6.6,
                color=C.CYANO_C if cy else "#2F3844",
                fontweight="semibold" if cy else "normal")
    for i in core_idx:
        tx.plot([0.0, 0.045], [i, i], color=ENR_C, lw=.6, alpha=.8)

    handles = [
        Line2D([], [], color=ENR_C, lw=2.2, label="core to the enrichments"),
        Line2D([], [], color=SHARED_C, lw=2.2, label="core to both"),
        Line2D([], [], color=SRC_C, lw=2.0, label="core to the reservoir"),
        Line2D([], [], color=NONE_C, lw=2.0, label="core to neither"),
    ]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.035, 0.012),
               ncol=2, frameon=False, fontsize=7.0, handlelength=1.9,
               columnspacing=1.3, handletextpad=.6,
               title="branch colour", title_fontsize=7.0, alignment="left")

    total = lengths.sum()
    stats = {k: lengths[[kind[id(nd)] == k for nd in nodes]].sum()
             for k in ("shared", "enrichment", "source", "none")}
    fig.text(0.560, 0.030,
             f"branch length: both {stats['shared']:.1f} ({100*stats['shared']/total:.0f}%), "
             f"enrichment {stats['enrichment']:.1f} ({100*stats['enrichment']/total:.0f}%),\n"
             f"reservoir {stats['source']:.1f} ({100*stats['source']/total:.0f}%). "
             f"Core UniFrac between the two cores 0.777.",
             fontsize=6.6, color="#5A6570", va="bottom", linespacing=1.5)

    fig.text(0.035, 0.982, "The core CyanoSphere phylogeny",
             fontsize=11.5, family="serif", va="top")
    fig.text(0.035, 0.958,
             f"RpL19 tree, {len(tips)} sequence clusters. A branch is core to a habitat when the clade below it is\n"
             f"detected in at least {CORE_THRESHOLD:.0%} of that habitat's eight samples. Enrichment-core tips are named at\n"
             f"family level, or genus for cyanobacteria; $\\times n$ marks a run of adjacent tips. Of the {len(tips)} tips,\n"
             f"36 are core to the enrichments, 190 to the reservoir, and one to both.",
             fontsize=6.9, color="#69737E", va="top", linespacing=1.7)

    C.save(fig, "fig_core_tree")
    print(pd.Series(stats).round(3).to_string())
    print(f"\nenrichment-core tips {len(core_idx)}, in {len(runs)} runs")


# ============================================================================
# Enrichment lineages through time (Fig 1B-E)
# ----------------------------------------------------------------------------
# Four lineages (DW S4-S7) were sequenced at P0 (first culture, 49 days after
# inoculation from the lake), P3 and P12 (the end-points M_132-M_135); lineage 2
# at P3 (M_62) and P6 (M_63); M_71 and M_77 at P0 (M_70, M_76) and P3. M_67 was
# only sequenced at its end-point. The 11 earlier-passage samples were profiled
# against the cluster references of the final RpL19 run with that run's settings
# (pipeline/10_fig1_lineages), so passages, end-points and reservoir are counted in
# one set of populations. The control (M_132 re-profiled from SRA) reproduces the
# end-point column of rpl19_cluster_coverage.csv with mean coverage (r = 1.0, max
# difference 0.42), which is therefore the coverage used here.
#     outputs/fig1_samples.csv, outputs/fig1_lineages.csv
# ============================================================================
from scipy.stats import wilcoxon  # noqa: E402

# lineage (named by its end-point) -> [(sample, passage)], from the sequencing log and SI S1.
LINEAGES = {"M_132": [("M_79", 0), ("M_83", 3), ("M_132", 12)],
            "M_133": [("M_80", 0), ("M_84", 3), ("M_133", 12)],
            "M_134": [("M_81", 0), ("M_85", 3), ("M_134", 12)],
            "M_135": [("M_82", 0), ("M_86", 3), ("M_135", 12)],
            "M_63": [("M_62", 3), ("M_63", 6)],
            "M_71": [("M_70", 0), ("M_71", 3)],
            "M_77": [("M_76", 0), ("M_77", 3)],
            "M_67": [("M_67", None)]}
# lineages are told apart by marker, not colour: everything enrichment is green
LINEAGE_MARKER = {"M_132": "o", "M_133": "s", "M_134": "^", "M_135": "D",
                  "M_63": "v", "M_71": "P", "M_77": "X", "M_67": "h"}
ENRICH_DARK = "#2f7d5b"
STAGE_X = {0: 1.2, 3: 2.2}      # intermediate passages; every end-point sits at X_END
X_SOURCE, X_END = 0.0, 3.4


def lineage_coverage() -> pd.DataFrame:
    """Reservoir, end-point and earlier-passage samples on the final RpL19 clusters."""
    cov = load_clusters()
    new = pd.read_csv(C.DATA / "fig1_lineages" / "coverage_mean-COVs.txt", sep="\t", index_col=0)
    new = new[[c for c in new.columns if c not in cov.columns]]   # M_132 was the control
    return cov.join(new.reindex(cov.index).fillna(0.0))


def sample_metrics(cov: pd.DataFrame) -> pd.DataFrame:
    """Richness, Shannon and cyanobacterial share per sample, with its lineage and
    passage. Cyanobacteria include unnamed clusters placed in the phylum on the tree."""
    m = diversity(cov)
    tax = C.rpl19_taxonomy()
    cyano = tax.index[tax.t_phylum == C.CYANO_PHYLUM]
    m["cyano"] = cov.loc[cyano].sum() / cov.sum()
    m["lineage"], m["passage"] = None, np.nan
    for lin, points in LINEAGES.items():
        for s, pas in points:
            m.loc[s, ["lineage", "passage"]] = [lin, pas]
    m["stage"] = ["Source" if s in C.SOURCE else "End-point" if s in C.ENRICHMENT else "Passage"
                  for s in m.index]
    return m


def panel_course(ax, m, col, ylabel, label_ends=True):
    """Source box | lineages through their earlier sequenced passages | end-point box. Each
    sample is drawn once: a lineage's line ends at its own point in the end-point box,
    whatever passage the end-point was sampled at (P3 to P12). p, end-point vs source
    Mann-Whitney."""
    src, end = m.loc[C.SOURCE, col].values, m.loc[C.ENRICHMENT, col].values
    for x, v, g in ((X_SOURCE, src, "Source"), (X_END, end, "Enrichment")):
        bp = ax.boxplot([v], positions=[x], widths=0.62, patch_artist=True, showfliers=False,
                        medianprops=dict(color="black", lw=1.0),
                        whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7))
        bp["boxes"][0].set(facecolor=C.GROUP_COLOR[g], alpha=0.28, edgecolor=C.GROUP_COLOR[g])
    rng = np.random.default_rng(0)
    ax.scatter(X_SOURCE + rng.uniform(-0.15, 0.15, len(src)), src, s=12,
               color=C.GROUP_COLOR["Source"], edgecolor="white", lw=0.4, zorder=3)
    jit = dict(zip(C.ENRICHMENT, np.linspace(-0.2, 0.2, len(C.ENRICHMENT))))
    ends = []
    for lin, points in LINEAGES.items():
        xs = [X_END + jit[s] if s in C.ENRICHMENT else STAGE_X[p] for s, p in points]
        ys = [m.loc[s, col] for s, _ in points]
        if len(xs) > 1:
            ax.plot(xs, ys, color=ENRICH_DARK, lw=0.8, alpha=0.75, zorder=2)
        ax.scatter(xs, ys, s=13, marker=LINEAGE_MARKER[lin], color=C.GROUP_COLOR["Enrichment"],
                   edgecolor=ENRICH_DARK, lw=0.5, zorder=3)
        ends.append((xs[-1], ys[-1], lin))
    ax.axvspan(STAGE_X[0] - 0.45, STAGE_X[3] + 0.45, color="#f2f5f3", zorder=0, lw=0)
    p = mannwhitneyu(src, end, alternative="two-sided").pvalue
    lo, hi = ax.get_ylim()
    y = max(src.max(), end.max(), m.loc[m.passage.notna(), col].max()) + 0.10 * (hi - lo)
    ax.plot([X_SOURCE, X_SOURCE, X_END, X_END], [y - 0.03 * (hi - lo), y, y, y - 0.03 * (hi - lo)],
            lw=0.6, color="#444444")
    ax.text((X_SOURCE + X_END) / 2, y + 0.015 * (hi - lo), f"$p$ = {p:.1e}", ha="center",
            va="bottom", fontsize=6.3)
    ax.set_ylim(lo, y + 0.12 * (hi - lo))
    if label_ends:
        from figures import label_column
        label_column(ax, ends, X_END + 0.5, "left", fontsize=5.2, min_gap=6.0)
    ax.set_xticks([X_SOURCE, *STAGE_X.values(), X_END], ["Source", "P0", "P3", "End-\npoint"])
    ax.set_xlim(X_SOURCE - 0.5, X_END + 1.25)
    ax.spines["bottom"].set_bounds(X_SOURCE - 0.5, X_END + 0.5)
    ax.set_ylabel(ylabel)
    ax.spines[["top", "right"]].set_visible(False)
    return p


def panel_stage_rank(ax, cov, m, top=60):
    """Mean rank-abundance curve per stage: source, P0, P3, end-points; each sample in one
    stage only (M_71 and M_77, end-points at P3, are end-points here)."""
    stages = [("Source", C.SOURCE, C.GROUP_COLOR["Source"]),
              ("P0", m.index[m.passage == 0].tolist(), "#b4dfc6"),
              ("P3", m.index[(m.passage == 3) & ~m.index.isin(C.ENRICHMENT)].tolist(), "#6fbf9b"),
              ("End-points", C.ENRICHMENT, ENRICH_DARK)]
    ends = []
    for name, members, col in stages:
        A = []
        for s in members:
            v = np.sort(cov[s].values)[::-1]
            A.append(100 * v[:top] / v.sum())
        A = np.vstack(A)
        if name in ("Source", "End-points"):
            for a in A:
                ax.plot(range(1, top + 1), a, color=col, alpha=0.18, lw=0.5)
        ax.plot(range(1, top + 1), A.mean(0), color=col, lw=1.6)
        ends.append((top, A.mean(0)[-1], f"{name} ($n$ = {len(A)})"))
    ax.set_yscale("log")
    ax.set_xlim(-0.5, top * 1.5)
    ax.set_xticks([1, 20, 40, 60])
    ax.spines["bottom"].set_bounds(1, top)
    ax.set_yticks([0.01, 0.1, 1, 10, 100], ["0.01", "0.1", "1", "10", "100"])
    ax.set_xlabel("Population rank")
    ax.set_ylabel("Relative coverage (%)")
    from figures import label_column
    label_column(ax, ends, top + 1.5, "left", fontsize=6.0, min_gap=6.5)
    ax.spines[["top", "right"]].set_visible(False)


def lineage_tests(m: pd.DataFrame) -> pd.DataFrame:
    """Paired change from P0 to P3 and to the end-point in the six lineages sequenced at P0
    (two-sided Wilcoxon signed-rank; n = 6, so the smallest attainable p is 0.031)."""
    rows = []
    p0 = m[m.passage == 0].set_index("lineage")
    for later, sel in (("P3", m.passage == 3), ("end-point", m.index.isin(C.ENRICHMENT))):
        q = m[sel].set_index("lineage").reindex(p0.index)
        for col in ("richness", "shannon", "cyano"):
            a, b = p0[col], q[col]
            rows.append(dict(comparison=f"P0 -> {later}", metric=col, n=int(b.notna().sum()),
                             median_P0=a.median(), median_later=b.median(),
                             fell=int((b < a).sum()), p=wilcoxon(a, b).pvalue))
    return pd.DataFrame(rows)


def run_lineage_trajectories():
    cov = lineage_coverage()
    m = sample_metrics(cov)
    m.to_csv(C.OUT / "fig1_samples.csv")
    t = lineage_tests(m)
    t.to_csv(C.OUT / "fig1_lineages.csv", index=False)
    print(m.round(3).to_string())
    print(t.round(4).to_string())
    return cov, m, t


# ============================================================================
# Convergence by level of organisation
# ----------------------------------------------------------------------------
# How similar are two enrichments, and two reservoir samples, at each level from
# phylum to 97% RpL19 population, and in function? Median pairwise Bray-Curtis of
# relative coverage at each level (unnamed clusters keep their own identity below
# phylum, and are placed by the tree at phylum). Plus: how many populations the
# eight samples of each group share, and how rare the populations that dominate the
# enrichments were in the reservoir. Note: reservoir samples are not the inocula,
# so "not detected" means below detection in these eight reservoir metagenomes.
#     outputs/convergence_by_level.csv, outputs/population_sharing.csv,
#     outputs/dominant_populations.csv
# ============================================================================
RANKS = ("t_phylum", "t_class", "t_order", "t_family", "t_genus")


def level_tables() -> dict:
    """level name -> features x samples relative coverage (taxonomy and function)."""
    rel = load_clusters(relative=True)
    tax = C.rpl19_taxonomy().reindex(rel.index)
    out = {}
    for rank in RANKS:
        key = tax[rank]
        if rank != "t_phylum":        # an unnamed cluster is its own unit below phylum
            key = key.where(key.notna(), pd.Series(rel.index, index=rel.index))
        out[rank[2:]] = rel.groupby(key).sum()
    sp = pd.read_csv(C.DATA / "rpl19_cluster_taxonomy.tsv", sep="\t",
                     index_col="split_name").reindex(rel.index).t_species
    out["species"] = rel.groupby(sp.where(sp.notna(), pd.Series(rel.index, index=rel.index))).sum()
    out["population"] = rel
    cov_m, cnt_m, _ = C.metric_matrices("enrichments", "rpl19")
    out["function, coverage"] = C.compositional(cov_m)
    out["function, count"] = C.compositional(cnt_m)
    return out


def convergence_by_level() -> pd.DataFrame:
    rows = []
    for level, m in level_tables().items():
        for g, members in (("Source", C.SOURCE), ("Enrichment", C.ENRICHMENT)):
            d = pdist(m[members].T.values, "braycurtis")
            rows.append(dict(level=level, group=g, median=np.median(d), min=d.min(), max=d.max()))
    return pd.DataFrame(rows)


def population_sharing(detect: float = 0.001) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Populations (>= `detect` relative coverage) shared by k of the 8 samples of each
    group and the share of each sample they make up; and the dominant enrichment
    populations (> 20% in some enrichment) with their reservoir detection."""
    rel = load_clusters(relative=True)
    tax = C.rpl19_taxonomy().reindex(rel.index)
    rows = []
    for g, members in (("Source", C.SOURCE), ("Enrichment", C.ENRICHMENT)):
        k = (rel[members] >= detect).sum(axis=1)
        for m in (8, 6, 4):
            core = k.index[k >= m]
            rows.append(dict(group=g, in_at_least=m, populations=len(core),
                             mean_share=rel.loc[core, members].sum().mean()))
    dom = sorted({c for s in C.ENRICHMENT for c in rel[s].index[rel[s] > 0.2]})
    d = pd.DataFrame([dict(population=c, phylum=tax.loc[c, "t_phylum"], genus=tax.loc[c, "t_genus"],
                           dominant_in=",".join(s for s in C.ENRICHMENT if rel.loc[c, s] > 0.2),
                           max_enrichment=rel.loc[c, C.ENRICHMENT].max(),
                           reservoir_detected=int((rel.loc[c, C.SOURCE] > 0).sum()),
                           reservoir_max=rel.loc[c, C.SOURCE].max()) for c in dom])
    return pd.DataFrame(rows), d


def unnamed_sensitivity(far: float = 0.6) -> pd.DataFrame:
    """Median shared coverage at phylum, class and order under four treatments of
    clusters without a GTDB phylum: placed on the tree (as used), placements further
    than `far` from a named tip lumped as 'unassigned', those kept as their own units,
    and no placement at all (every unnamed cluster its own unit)."""
    rel = load_clusters(relative=True)
    t = C.rpl19_taxonomy().reindex(rel.index)
    own = pd.Series(rel.index, index=rel.index)
    is_far = t.placed & (t.patristic > far)
    rows = []
    for name in ("placed on tree", "far as one unassigned bin", "far as own units", "no placement"):
        for rank in ("t_phylum", "t_class", "t_order"):
            k = t[rank].copy()
            if name == "far as one unassigned bin":
                k[is_far] = "unassigned"
            elif name == "far as own units":
                k[is_far] = own[is_far]
            elif name == "no placement":
                k[t.placed] = own[t.placed]
            m = rel.groupby(k.where(k.notna(), own)).sum()
            for g, members in (("Source", C.SOURCE), ("Enrichment", C.ENRICHMENT)):
                rows.append(dict(treatment=name, level=rank[2:], group=g,
                                 median_shared=np.median(1 - pdist(m[members].T.values, "braycurtis"))))
    out = pd.DataFrame(rows)
    out.to_csv(C.OUT / "unnamed_sensitivity.csv", index=False)
    return out


def run_convergence_by_level():
    lv = convergence_by_level()
    sh, dom = population_sharing()
    lv.to_csv(C.OUT / "convergence_by_level.csv", index=False)
    sh.to_csv(C.OUT / "population_sharing.csv", index=False)
    dom.to_csv(C.OUT / "dominant_populations.csv", index=False)
    print(unnamed_sensitivity().pivot_table(index=["treatment", "level"], columns="group",
                                           values="median_shared").round(2).to_string())
    print(lv.pivot(index="level", columns="group", values="median").round(2).to_string())
    print(sh.round(2).to_string(index=False))
    print(dom.round(4).to_string(index=False))
    return lv, sh, dom


def main() -> None:
    run_diversity_and_taxonomy()
    run_lineage_trajectories()
    run_convergence_by_level()
    run_phylogenetic_core()
    run_divergence()
    run_core_tree_figure()


if __name__ == "__main__":
    main()
