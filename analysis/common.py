"""
Shared definitions for the CyanoSphere analysis.

Everything that more than one script needs lives here: where the data is, how the
two functional metrics are defined, how differential tests are done, and the
figure style. Nothing here reads or writes outside the repository.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
OUT = REPO / "outputs"
OUT.mkdir(exist_ok=True)

# --- analysis constants -------------------------------------------------------
# A module counts as present in a sample at this stepwise completeness. Locked
# for the whole paper; changing it changes every number downstream.
PRESENCE = 0.75
# A MAG counts as encoding a module at this pathwise completeness (Fig 3E).
ENCODING = 0.75
# Modules whose pathwise completeness can pass ENCODING on enzymes shared with
# other pathways. The Calvin cycle reaches 9/11 steps in most heterotrophs
# through glycolysis and the pentose phosphate pathway; a genome encodes it only
# with RuBisCO and phosphoribulokinase. Nitrogen fixation needs nifH, nifD, nifK.
SIGNATURE_KOS = {"M00165": {"K01601", "K00855"},
                 "M00175": {"K02588", "K02586", "K02591"}}
# A MAG counts as present in a sample at this mean read coverage (Fig 3E, S4).
MAG_DETECT = 1.0
# Significance for the differential tests.
FDR = 0.05
MIN_L2FC = 0.5

SOURCE = ["W_01", "W_06", "W_07", "W_10", "W_13", "W_14", "W_15", "W_18"]
ENRICHMENT = ["M_63", "M_67", "M_71", "M_77", "M_132", "M_133", "M_134", "M_135"]
MEDIA_V = [f"V_P10_R{i}" for i in range(1, 6)]
MEDIA_NV = [f"NV_P10_R{i}" for i in (1, 2, 4, 5)]   # R3: only 2 bacterial populations (eukaryote-dominated)
MEDIA_NVNN = [f"NVNN_P10_R{i}" for i in range(1, 5)]

CYANO_PHYLUM = "Cyanobacteriota"

# Modules that require anoxic or microaerobic conditions, or that are inhibited
# by oxygen. Fixed in advance of testing, so the set-level test below is a
# hypothesis test and not a search. Sources: KEGG module definitions; the
# oxygen sensitivity of nitrogenase; cbb3 and bd oxidases as high-affinity,
# microaerobic terminal oxidases.
# Split into two sets, because the nitrogen modules would respond to a change in
# nitrogen supply whatever the oxygen regime, and testing them together with the
# oxygen-related ones would confound the two. OXYGEN_MODULES is the set that
# actually tests whether nitrate-free cultures develop more anoxic or
# microaerobic habitat; N_ANAEROBIC_MODULES is reported separately.
OXYGEN_MODULES = {
    "M00597": "Anoxygenic photosystem II",
    "M00612": "Anoxygenic photosynthesis, purple bacteria",
    "M00613": "Anoxygenic photosynthesis, green nonsulfur bacteria",
    "M00156": "Cytochrome c oxidase, cbb3-type (high affinity)",
    "M00153": "Cytochrome bd ubiquinol oxidase (high affinity)",
    "M00985": "Sulfide oxidation",
    "M00595": "Thiosulfate oxidation, SOX complex",
    "M00596": "Dissimilatory sulfate reduction",
    "M00579": "Phosphate acetyltransferase-acetate kinase (fermentative)",
    "M00357": "Methanogenesis, acetate to methane",
    "M00567": "Methanogenesis, CO2 to methane",
}
N_ANAEROBIC_MODULES = {
    "M00175": "Nitrogen fixation (nitrogenase is oxygen-sensitive)",
    "M00529": "Denitrification",
    "M00530": "Dissimilatory nitrate reduction",
}

ANOXIC_MODULES = {
    "M00597": "Anoxygenic photosystem II",
    "M00612": "Anoxygenic photosynthesis, purple bacteria",
    "M00613": "Anoxygenic photosynthesis, green sulfur bacteria",
    "M00156": "Cytochrome c oxidase, cbb3-type (high affinity)",
    "M00153": "Cytochrome bd ubiquinol oxidase (high affinity)",
    "M00985": "Sulfide oxidation",
    "M00595": "Thiosulfate oxidation, SOX complex",
    "M00596": "Dissimilatory sulfate reduction",
    "M00529": "Denitrification",
    "M00530": "Dissimilatory nitrate reduction",
    "M00175": "Nitrogen fixation",
    "M00579": "Phosphate acetyltransferase-acetate kinase (fermentative)",
    "M00357": "Methanogenesis, acetate to methane",
    "M00567": "Methanogenesis, CO2 to methane",
}

# --- figure style -------------------------------------------------------------
GROUP_COLOR = {"Source": "#3b4a6b", "Enrichment": "#6fbf9b",
               "V": "#2e7d32", "NV": "#f08c00", "NVNN": "#c62828"}
# Phylum colours of the EcoPhylo heatmap (Fig 2A, ColorBrewer Paired), reused in Fig 2B-C
PHYLUM_COLOR = {"Cyanobacteriota": "#33a02b", "Pseudomonadota": "#e21a1c",
                "Actinomycetota": "#a5cde2", "Bacteroidota": "#1f77b3",
                "Verrucomicrobiota": "#ff7e00", "Armatimonadota": "#fb9a99",
                "Planctomycetota": "#cab1d6", "Spirochaetota": "#4d7908",
                "Bacillota": "#b15928", "Gemmatimonadota": "#fdbf6e",
                "Chloroflexota": "#b1df89", "Other": "#bdbdbd", "Unknown": "#000000"}
CYANO_C = "#2e7d32"
HETERO_C = "#1565c0"
CLASS_COLOR = {"both metrics": "#37474f",
               "count-based only": "#b71c1c",
               "coverage-based only": "#8d6e63"}


def style() -> None:
    mpl.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9, "axes.titlesize": 10.5, "axes.titleweight": "bold",
        "axes.labelsize": 9.5, "axes.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "black", "axes.labelcolor": "black",
        "xtick.color": "black", "ytick.color": "black",
        "legend.frameon": False, "legend.fontsize": 8,
        "figure.facecolor": "white", "axes.facecolor": "white",
    })


def polish(fig, min_size: float = 6.0) -> None:
    """House style for every figure: grey text becomes black (coloured text is kept), and no
    text is smaller than `min_size` points."""
    from matplotlib.colors import to_rgb
    from matplotlib.text import Text
    for t in fig.findobj(Text):
        if not t.get_text():
            continue
        r, g, b = to_rgb(t.get_color())
        if max(r, g, b) - min(r, g, b) < 0.08 and 0.05 < (r + g + b) / 3 < 0.95:
            t.set_color("black")
        if t.get_fontsize() < min_size:
            t.set_fontsize(min_size)


def crop_pdf(path, margin: float = 4.0) -> None:
    """Crop a one-page PDF to its drawn content plus `margin` points."""
    import pymupdf
    doc = pymupdf.open(path)
    page = doc[0]
    area = page.rect.width * page.rect.height
    boxes = [pymupdf.Rect(b[1]) for b in page.get_bboxlog()
             if b[0] != "ignore-text" and pymupdf.Rect(b[1]).width * pymupdf.Rect(b[1]).height < 0.9 * area
             and not pymupdf.Rect(b[1]).is_empty]
    if boxes:
        r = boxes[0]
        for b in boxes[1:]:
            r |= b
        r = (r + (-margin, -margin, margin, margin)) & page.rect
        page.set_cropbox(r)
    doc.save(path, incremental=True, encryption=0)
    doc.close()


def save(fig, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.pdf")
    fig.savefig(OUT / f"{stem}.svg")
    print(f"wrote outputs/{stem}.pdf and .svg")


# --- data ---------------------------------------------------------------------
def load_modules(which: str) -> pd.DataFrame:
    """`which` is 'enrichments' or 'media'. Long format, one row per module x sample."""
    d = pd.read_csv(DATA / f"kegg_modules_{which}.csv")
    return d[d.stepwise_module_completeness >= PRESENCE].copy()


def population_counts(which: str, marker: str) -> pd.Series:
    """
    Denominator for the count-based metric: the number of coexisting populations
    in a sample, where a population is one distinct 97% ribosomal-protein
    sequence cluster. `marker` is one of:
        rpl19      RpL19 EcoPhylo clusters with non-zero coverage (enrichments
                   and reservoir only). The same count as population richness
                   in Fig 1B, so "population" means one thing in the paper.
        rps16      RpS16 EcoPhylo clusters, detection > 0.25 (both experiments;
                   the media experiment's denominator, the marker used for
                   composition there; RpL19 gives the same cyanobacteria once a
                   plastid-derived RpS16 cluster is set aside).
        mags       dereplicated MAGs detected (enrichments only).
        anvio_l19  the `bacterial_genomes_L19` column of the anvi'o metabolism
                   export. Used as the denominator until Sep 2026; an
                   assembly-level estimate whose exact rule is not recorded
                   (most likely the RpL19 HMM hit count behind anvi'o's
                   per-population copy number). Kept for the robustness table.
    """
    if marker == "rpl19":
        if which != "enrichments":
            raise ValueError("RpL19 clusters exist for the enrichment experiment only")
        cl = pd.read_csv(DATA / "rpl19_cluster_coverage.csv", index_col=0)
        cl.columns = [c.lstrip("\ufeff") for c in cl.columns]
        return (cl > 0).sum()
    d = pd.read_csv(DATA / f"kegg_modules_{which}.csv")
    if marker == "anvio_l19":
        return d.groupby("genome_name")["bacterial_genomes_L19"].first()
    if marker == "mags":
        s = d.groupby("genome_name")["N_MAGs"].first()
        if s.isna().all():
            raise ValueError("no MAG counts for this experiment")
        return s
    if marker == "rps16":
        p = pd.read_csv(DATA / "rps16_population_counts.csv")
        return p.set_index("sample")["n_rps16_populations"]
    raise ValueError(f"unknown marker {marker!r}")


def metric_matrices(which: str, marker: str = "rpl19"):
    """
    Return (coverage, count, names): two module x sample matrices and the module
    descriptions.

    coverage-based  per sample, anvi'o's module coverage (`avg_coverage`, the mean
                    read coverage of the genes annotated to the module), for
                    modules at stepwise completeness >= PRESENCE. Raw, not
                    depth-normalised: follows how abundant the carriers are, and
                    also how deeply the sample was sequenced.
    count-based     per sample, the module's stepwise copy number divided by the
                    number of coexisting populations (anvi'o's per-population
                    copy number, with the denominator from `population_counts`).
                    A proxy for how widely the function is spread, largely
                    independent of how abundant its carriers are.

    Both are computed per sample; groups are compared afterwards.
    """
    d = load_modules(which)
    npop = population_counts(which, marker)
    cov = d.pivot_table(index="module", columns="genome_name",
                        values="avg_coverage", aggfunc="max")
    cnt = d.pivot_table(index="module", columns="genome_name",
                        values="stepwise_copy_number", aggfunc="max")
    cnt = cnt.div(npop.reindex(cnt.columns).astype(float), axis=1)
    names = d.drop_duplicates("module").set_index("module")[["module_name", "module_category"]]
    return cov.fillna(0.0), cnt.fillna(0.0), names


def mag_tables():
    """
    The one MAG identity table used everywhere (Fig 3E, S4, MAG control).

    Returns (completeness, coverage, taxonomy):
        completeness  module x MAG pathwise completeness
        coverage      MAG x sample mean read coverage, linked to each MAG by
                      sequence (contig headers, see pipeline/05_link_mags_to_bins.py)
        taxonomy      MAG x rank GTDB taxonomy (mag_quality.tsv)
    restricted to the 88 dereplicated MAGs present in all three. A MAG is
    detected in a sample at coverage >= MAG_DETECT.
    """
    comp = pd.read_csv(DATA / "mag_module_pathwise_completeness.tsv", sep="\t").set_index("module")
    cov = pd.read_csv(DATA / "mag_sample_coverage.tsv", sep="\t").set_index("MAG_name")
    tax = pd.read_csv(DATA / "mag_quality.tsv", sep="\t").set_index("MAG_name")
    mags = [g for g in comp.columns if g in cov.index and g in tax.index]
    return comp[mags], cov.loc[mags], tax.loc[mags]


# --- statistics ---------------------------------------------------------------
def log2fc(a: np.ndarray, b: np.ndarray, all_values: np.ndarray) -> float:
    """
    Group-mean log2 fold change with a per-module half-minimum pseudocount.

    The pseudocount is half the smallest non-zero value of that module across all
    samples in the comparison. A single global pseudocount (the alternative) gives
    modules that are absent from one group an arbitrarily large fold change that
    depends on the smallest value anywhere in the dataset; a per-module one keeps
    the fold change on the module's own scale. This convention is used for every
    fold change in the paper.
    """
    pos = all_values[all_values > 0]
    eps = pos.min() / 2 if len(pos) else 1e-9
    return float(np.log2((a.mean() + eps) / (b.mean() + eps)))


def _exact_ranksum_p(a: np.ndarray, b: np.ndarray, n_mc: int = 200_000,
                     seed: int = 0) -> float:
    """
    Two-sided exact permutation p for the Wilcoxon rank-sum statistic.

    Mid-ranks are used for ties and the null distribution is enumerated over
    every assignment of the pooled values to the two groups (12,870 for 8 vs 8,
    715 for 4 vs 9), so ties are handled exactly rather than by the normal
    approximation SciPy falls back on when values tie. The statistic is the
    distance of group 1's rank sum from its expectation. Falls back to Monte
    Carlo only if the enumeration would exceed n_mc assignments.
    """
    from itertools import combinations
    from math import comb
    from scipy.stats import rankdata
    x = np.concatenate([a, b])
    r = rankdata(x)
    n, k = len(x), len(a)
    expected = k * (n + 1) / 2
    obs = abs(r[:k].sum() - expected)
    if comb(n, k) <= n_mc:
        idx = np.array(list(combinations(range(n), k)))
    else:
        rng = np.random.default_rng(seed)
        idx = np.array([rng.permutation(n)[:k] for _ in range(n_mc)])
    null = np.abs(r[idx].sum(axis=1) - expected)
    return float((null >= obs - 1e-9).mean())


def differential(mat: pd.DataFrame, g1: list[str], g2: list[str]) -> pd.DataFrame:
    """
    Two-sided exact rank-sum test per module, Benjamini-Hochberg across the
    modules tested on that metric.

    A module is tested when it is present (non-zero) in at least one sample of
    the comparison, so each metric has its own family size. Rank-based because
    group sizes are 4 to 9 and the metrics are not normal; exact because at
    these sizes, with ties, the normal approximation is poor. `separated`
    marks modules where every g1 value lies outside the range of every g2 value.
    """
    g1 = [s for s in g1 if s in mat.columns]
    g2 = [s for s in g2 if s in mat.columns]
    rows = []
    for m in mat.index:
        a = mat.loc[m, g1].to_numpy(float)
        b = mat.loc[m, g2].to_numpy(float)
        if not (a.any() or b.any()):
            continue
        rows.append(dict(module=m, l2fc=log2fc(a, b, mat.loc[m].to_numpy(float)),
                         p=_exact_ranksum_p(a, b),
                         p_asymptotic=mannwhitneyu(a, b, alternative="two-sided").pvalue,
                         mean_g1=a.mean(), mean_g2=b.mean(),
                         separated=bool(a.min() > b.max() or a.max() < b.min())))
    out = pd.DataFrame(rows).set_index("module")
    out["padj"] = multipletests(out.p, method="fdr_bh")[1]
    out["padj_asymptotic"] = multipletests(out.p_asymptotic, method="fdr_bh")[1]
    out["sig"] = (out.padj < FDR) & (out.l2fc.abs() > MIN_L2FC)
    return out


def niche_class(cov: pd.DataFrame, cnt: pd.DataFrame) -> pd.Series:
    """
    The two-metric classification, applied to a pair of differential tables.

    both metrics         increased in abundance and in how widely it is held:
                         the behaviour expected of a function the medium imposes
    count-based only     more widely held without becoming more abundant:
                         the behaviour expected of a function selected by
                         conditions the community generates for itself
    coverage-based only  more abundant but no more widely held
    neither              not enriched

    "Increased" means FDR < 0.05 and log2FC > +0.5 on that metric, the same
    threshold used everywhere else in the paper.
    """
    up_cov = (cov.padj < FDR) & (cov.l2fc > MIN_L2FC)
    up_cnt = (cnt.padj < FDR) & (cnt.l2fc > MIN_L2FC)
    idx = cov.index.union(cnt.index)
    up_cov = up_cov.reindex(idx, fill_value=False)
    up_cnt = up_cnt.reindex(idx, fill_value=False)
    return pd.Series(np.select(
        [up_cov & up_cnt, ~up_cov & up_cnt, up_cov & ~up_cnt],
        ["both metrics", "count-based only", "coverage-based only"],
        default="neither"), index=idx, name="niche_class")


def set_test(table: pd.DataFrame, members: set[str], label: str) -> dict:
    """
    Is a pre-specified set of modules shifted relative to all the others?

    Testing a curated set as a set, rather than each module separately, is the
    right test when per-module power is low: it asks one question instead of a
    few hundred, so nothing is spent on multiple-testing correction. Two-sided
    rank-sum on the log2 fold changes, plus the median of each group.
    """
    inside = table.loc[table.index.isin(members), "l2fc"].dropna()
    outside = table.loc[~table.index.isin(members), "l2fc"].dropna()
    row = dict(set=label, n_in=len(inside), n_out=len(outside),
               median_in=float(inside.median()) if len(inside) else np.nan,
               median_out=float(outside.median()) if len(outside) else np.nan,
               p=np.nan)
    if len(inside) >= 3:
        row["p"] = float(mannwhitneyu(inside, outside, alternative="two-sided").pvalue)
    return row


def compositional(mat: pd.DataFrame) -> pd.DataFrame:
    """
    Per-sample normalisation of a coverage matrix: every sample's module
    coverages are divided by that sample's total.

    Used for the ordinations (Fig 3B, Fig 4C) and nowhere else. An ordination of
    raw coverage is dominated by how deeply a sample was sequenced and how much
    biomass it had; normalising removes that and leaves the shape of the
    functional profile, which is what the panel is about. The differential tests
    and the two-metric classification deliberately keep raw coverage, because the
    point of pairing the two metrics is that one is abundance-sensitive and the
    other is not. Normalising the coverage metric collapses that distinction:
    on this dataset it moves 45 of the 60 modules enriched on both metrics into
    the count-based-only class, which is an artefact of the normalisation and not
    a result.
    """
    return mat.div(mat.sum(axis=0), axis=1)


def bray_pcoa(mat: pd.DataFrame, groups: dict, seed: int = 0) -> dict:
    """
    PCoA on Bray-Curtis dissimilarity between samples (columns of a relative-abundance
    matrix), with PERMANOVA (location) and PERMDISP (dispersion) on the same distances
    and each group's mean distance to its centroid. Used for the source vs enrichment
    ordinations of taxonomy (Fig 2C-D), function (Fig 3B) and the media (Fig 4C-E), so
    that the spread of a group in the plot reflects how different its samples are: a
    PCA on log abundance instead pulls low-richness samples together, whatever their
    composition. supplementary.ordination_robustness() compares the alternatives.
    """
    from skbio.diversity import beta_diversity
    samples = [s for g in groups.values() for s in g if s in mat.columns]
    d = beta_diversity("braycurtis", mat[samples].T.values.astype(float), ids=samples)
    return pcoa_stats(d, groups, seed)


def pcoa_stats(d, groups: dict, seed: int = 0) -> dict:
    """PCoA of a skbio DistanceMatrix, PERMANOVA and PERMDISP on it, and each sample's
    distance to its group centroid in the full PCoA space (Anderson 2006). A PCA is the
    special case of Euclidean distances."""
    from skbio.stats.ordination import pcoa
    from skbio.stats.distance import permanova as sk_permanova, permdisp
    samples = [s for g in groups.values() for s in g if s in d.ids]
    labels = [k for k, g in groups.items() for s in g if s in d.ids]
    d = d.filter(samples)
    o = pcoa(d, number_of_dimensions=2)
    pm = sk_permanova(d, labels, permutations=9999, seed=seed)
    pdisp = permdisp(d, labels, test="centroid", permutations=9999, seed=seed)
    X = pcoa(d).samples.values
    lab = np.array(labels)
    dist = pd.Series(0.0, index=samples)
    for k in groups:
        dist[lab == k] = np.linalg.norm(X[lab == k] - X[lab == k].mean(0), axis=1)
    cent = {k: dist[lab == k].mean() for k in groups}
    k, n = len(groups), len(samples)      # PERMANOVA R2 = SS_between / SS_total, from F
    r2 = pm["test statistic"] * (k - 1) / (pm["test statistic"] * (k - 1) + n - k)
    return dict(scores=pd.DataFrame(o.samples.values[:, :2], index=samples, columns=["PCo1", "PCo2"]),
                explained=100 * o.proportion_explained.values[:2], F=pm["test statistic"], R2=r2,
                p=pm["p-value"], disp_p=pdisp["p-value"], centroid_dist=cent, sample_dist=dist,
                pairwise=pairwise_permanova(d, groups) if len(groups) > 2 else {})


def pairwise_permanova(d, groups: dict) -> dict:
    """PERMANOVA p for every pair of groups on the same distances, Benjamini-Hochberg adjusted."""
    import itertools
    from skbio.stats.distance import permanova as sk_permanova
    out = {}
    for a, b in itertools.combinations(groups, 2):
        sa = [s for s in groups[a] if s in d.ids]; sb = [s for s in groups[b] if s in d.ids]
        out[(a, b)] = sk_permanova(d.filter(sa + sb), [a] * len(sa) + [b] * len(sb),
                                   permutations=9999, seed=0)["p-value"]
    if len(out) > 1:     # Benjamini-Hochberg across the pairs
        from statsmodels.stats.multitest import multipletests
        out = dict(zip(out, multipletests(list(out.values()), method="fdr_bh")[1]))
    return out


def permanova(mat: pd.DataFrame, groups: dict, n_perm: int = 999, seed: int = 0):
    """
    Permutational pseudo-F on Euclidean distances between two or more groups of
    samples, with the samples as rows of a z-scored module matrix. Returns
    (F, p). Used for every group-separation test in the paper.
    """
    from sklearn.preprocessing import StandardScaler
    rng = np.random.default_rng(seed)
    samples, labels = [], []
    for i, (_, members) in enumerate(groups.items()):
        for s_ in members:
            if s_ in mat.columns:
                samples.append(s_)
                labels.append(i)
    labels = np.array(labels)
    X = StandardScaler().fit_transform(mat[samples].T.values)
    D = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
    k = len(set(labels))

    def pseudo_f(lab):
        n = len(lab)
        total = (D ** 2).sum() / (2 * n)
        within = 0.0
        for g in set(lab):
            idx = np.where(lab == g)[0]
            within += (D[np.ix_(idx, idx)] ** 2).sum() / (2 * len(idx))
        return ((total - within) / (k - 1)) / (within / (n - k))

    obs = pseudo_f(labels)
    null = [pseudo_f(rng.permutation(labels)) for _ in range(n_perm)]
    return float(obs), (1 + sum(f >= obs for f in null)) / (n_perm + 1)


def tex_escape(s) -> str:
    return (str(s).replace("&", r"\&").replace("%", r"\%")
            .replace("_", r"\_").replace("=>", r"$\Rightarrow$"))


def mag_representation() -> pd.Series:
    """
    Per enrichment, cyanobacterial coverage in detected MAGs as a percentage of
    cyanobacterial RpL19 coverage (capped at 100): how much of the
    cyanobacterial community the MAG-based attribution can see.
    """
    _, cov, tax = mag_tables()
    clusters = pd.read_csv(DATA / "rpl19_cluster_coverage.csv", index_col=0)
    ctax = rpl19_taxonomy().reindex(clusters.index)
    marker = clusters.loc[ctax.t_phylum == CYANO_PHYLUM, ENRICHMENT].sum()
    in_mags = pd.Series({s: cov.loc[(cov[s] >= MAG_DETECT) & (tax.phylum == CYANO_PHYLUM), s].sum()
                         for s in ENRICHMENT})
    return (100 * in_mags / marker).clip(upper=100)


# --- EcoPhylo cluster taxonomy, with tree placement ---------------------------
def ecophylo_taxonomy(run: str) -> pd.DataFrame:
    """
    GTDB taxonomy of every cluster in an EcoPhylo run (data/ecophylo_clusters.csv).

    Clusters whose representative has no SCG taxonomy are placed on the run's
    marker tree (data/ecophylo_trees/<run>.nwk): they take the phylum and family
    of the nearest named tip (smallest patristic distance), never its genus, and
    are flagged `placed`. In the media experiment these unnamed clusters carry a
    large share of the coverage (up to ~40% in NV), including a dominant
    cyanobacterium, so dropping them would misstate composition.
    """
    k = pd.read_csv(DATA / "ecophylo_clusters.csv")
    k = k[k.run == run].drop_duplicates("cluster").set_index("cluster")
    out = k[["t_phylum", "t_class", "t_order", "t_family", "t_genus"]].copy()
    return place_on_tree(out, DATA / "ecophylo_trees" / f"{run}.nwk")


def rpl19_taxonomy() -> pd.DataFrame:
    """GTDB taxonomy of the clusters of the final RpL19 run (source vs enrichment,
    Fig 1-2), with unnamed clusters placed on that run's tree."""
    cov = pd.read_csv(DATA / "rpl19_cluster_coverage.csv", index_col=0)
    t = pd.read_csv(DATA / "rpl19_cluster_taxonomy.tsv", sep="\t", index_col="split_name")
    out = t.reindex(cov.index)[["t_phylum", "t_class", "t_order", "t_family", "t_genus"]]
    return place_on_tree(out, DATA / "rpl19_cluster_tree.nwk")


def place_on_tree(out: pd.DataFrame, tree_path: Path) -> pd.DataFrame:
    """Unnamed clusters take phylum to family (never genus) of the nearest named tip."""
    from skbio import TreeNode
    out = out.copy()
    out["placed"], out["patristic"] = False, 0.0
    tree = TreeNode.read([tree_path.read_text()])
    tips = {t.name.replace(" ", "_"): t for t in tree.tips()}  # newick reads unquoted _ as space
    named = [c for c in out.index if isinstance(out.loc[c, "t_phylum"], str) and c in tips]
    for c in out.index[out.t_phylum.isna()]:
        if c not in tips or not named:
            continue
        dist = {n: tips[c].distance(tips[n]) for n in named}
        best = min(dist, key=dist.get)
        out.loc[c, ["t_phylum", "t_class", "t_order", "t_family"]] = out.loc[best, ["t_phylum", "t_class", "t_order", "t_family"]].values
        out.loc[c, ["placed", "patristic"]] = [True, dist[best]]
    return out
