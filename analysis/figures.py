"""
Figure 3 B, D, E and Figure 4 C, D in the style of the published panels.

Run:  python figures.py      (from analysis/, after the modules before it in
      the README order)
"""
from __future__ import annotations

import re

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import common as C

# ============================================================================
# figure_panels.py
# ----------------------------------------------------------------------------
# Figure 3 B, D, E and Figure 4 C, D, drawn in the style of the published panels
# (boxed axes, centred bold titles, the Source/Enrichment and V/NV/NVNN palettes,
# split-triangle cells) from the current analysis.
#
# Each panel is written on its own, sized for the Affinity composites:
#     outputs/panels/fig3B.pdf   functional community structure (PCA)
#     outputs/panels/fig3D.pdf   the two metrics against each other, by class
#     outputs/panels/fig3D_volcano.pdf   the same data as a coverage volcano
#     outputs/panels/fig3E.pdf   who encodes each function, per replicate
#     outputs/panels/fig4C.pdf   functional community structure, media experiment
#     outputs/panels/fig4D.pdf   response to nitrate withdrawal (count-based)
# plus the CSV behind each. Labels never overlap: they are stacked in a column
# beside the points they name, with a thin leader line.
#
# Run after 02 and 03.
# ============================================================================

import function as enr  # noqa: E402


import media as med  # noqa: E402


OUT = C.OUT / "panels"


OUT.mkdir(exist_ok=True)


# the published palettes
SRC, ENR = "#55558f", "#8fd0b0"


MEDIA = {"V": "#6aae6a", "NV": "#f28468", "NVNN": "#e377b8"}


CLS = {"both metrics": "#2f3e46", "count-based only": "#c0392b",
       "coverage-based only": "#b8a48a", "neither": "#d5d8dc"}


GRID = "#e6e6e6"


# order of the published Fig 3E columns
E_ORDER = ["M_132", "M_133", "M_134", "M_135", "M_63", "M_67", "M_71", "M_77"]


# Fig 3E rows, by function, as in the published panel
E_GROUPS = [
    ("Photo-\nsynthesis", [("M00161", "Photosystem II"), ("M00163", "Photosystem I"),
                           ("M00597", "Anoxygenic photosystem II")]),
    ("Respiration", [("M00156", "Cytochrome c oxidase (cbb3)"),
                     ("M00417", "Cytochrome o oxidase")]),
    ("Nitrogen", [("M00615", "Nitrate assimilation"),
                  ("M00531", "Assimilatory nitrate reduction"),
                  ("M00529", "Denitrification")]),
    ("Vitamins and\nco-factors", [("M00122", "Cobalamin (aerobic)"),
                                  ("M00924", "Cobalamin (anaerobic)"),
                                  ("M00123", "Biotin"),
                                  ("M00899", "Thiamine salvage"),
                                  ("M00116", "Menaquinone")]),
    ("Amino acids", [("M00134", "Polyamine (ornithine)"),
                     ("M00133", "Polyamine (arginine)"),
                     ("M00970", "Proline degradation")]),
    ("Osmolytes and\ndefence", [("M00555", "Betaine"),
                                ("M00118", "Glutathione"),
                                ("M00642", "Multidrug efflux (MexJK-OprM)")]),
]


E_LABEL = {m: lab for _, rows in E_GROUPS for m, lab in rows}


# short names for labelled points; everything else uses the first KEGG clause
SHORT = dict(E_LABEL)


SHORT.update({
    "M00930": "Menaquinone via futalosine", "M00931": "Menaquinone via futalosine (modified)",
    "M00162": "Cytochrome b6f", "M00145": "NDH-1 complex",
    "M00932": "Phylloquinone", "M00112": "Tocopherol", "M00097": "beta-Carotene",
    "M00611": "Oxygenic photosynthesis", "M00165": "Calvin cycle",
    "M00175": "Nitrogen fixation", "M00580": "Pentose phosphate (archaeal)",
    "M00345": "Formaldehyde assimilation (RuMP)", "M00950": "Biotin (BioU)",
    "M00912": "NAD biosynthesis", "M00119": "Pantothenate",
    "M00064": "ADP-heptose", "M00612": "Anoxygenic photosynthesis (purple)",
    "M00034": "Methionine salvage", "M00133": "Polyamine (arginine)",
    "M00046": "Pyrimidine degradation",
    "M00642": "Multidrug efflux (MexJK-OprM)", "M00156": "Cytochrome cbb3 oxidase",
    "M00597": "Anoxygenic photosystem II",
    "M00530": "Dissimilatory nitrate reduction", "M00163": "Photosystem I",
})


def fmt_p(p, name="p"):
    """$p$ = 0.012, or $p$ < 0.001 (permutation p-values never reach 0)."""
    return f"${name}$ = {p:.3f}" if p >= 0.001 else f"${name}$ < 0.001"


def short(m, names, brackets=True):
    s = SHORT.get(m, str(names.get(m, m)).split(",")[0])
    return s if brackets else re.sub(r"\s*\(.*?\)", "", s)


def pub_style():
    C.style()
    mpl.rcParams.update({
        "axes.spines.top": True, "axes.spines.right": True,
        "axes.linewidth": 1.0, "axes.edgecolor": "#222222",
        "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlepad": 8,
        "axes.labelsize": 10, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    })


def ellipse(ax, xy, color, dashed=False, n_std=2.0):
    """Covariance ellipse at n_std, filled and edged, as in the published PCAs."""
    if len(xy) < 3:
        return
    vals, vecs = np.linalg.eigh(np.cov(xy.T))
    ang = np.degrees(np.arctan2(vecs[1, -1], vecs[0, -1]))
    w, h = 2 * n_std * np.sqrt(vals[::-1])
    for fc, ec, lw in ((color, "none", 0), ("none", color, 1.6 if dashed else 0.9)):
        ax.add_patch(mpatches.Ellipse(xy.mean(0), w, h, angle=ang, facecolor=fc,
                                      alpha=0.18 if fc != "none" else 1, edgecolor=ec, lw=lw,
                                      ls=(0, (5, 3)) if dashed and lw else "-", zorder=1 if lw == 0 else 2))


def label_column(ax, pts, x_col, ha, fontsize=7.2, min_gap=None, color="#222222",
                 y_range=None):
    """
    Stack labels in one column at x_col, each as close as possible to its
    point's height, never closer than min_gap to the next (in points), with a
    leader line to the point. `pts` is a list of (x, y, text).
    """
    if not pts:
        return
    if y_range is not None:
        # Labels in a block above or below their points. The label nearest the
        # points goes to the point nearest the column, so leader lines nest
        # instead of crossing.
        above = np.mean(y_range) > np.mean([p[1] for p in pts])
        pts = sorted(pts, key=lambda p: abs(p[0] - x_col))
        ys = np.linspace(*y_range, len(pts)) if len(pts) > 1 else [np.mean(y_range)]
        ys = ys if above else ys[::-1]
        for (x, y, text), ty in zip(pts, ys):
            ax.annotate(text, xy=(x, y), xytext=(x_col, ty), ha=ha, va="center",
                        fontsize=fontsize, color=color, zorder=6,
                        arrowprops=dict(arrowstyle="-", lw=0.5, color="#999999",
                                        shrinkA=1.5, shrinkB=3))
        return
    fig = ax.figure
    fig.canvas.draw()
    to_disp = ax.transData.transform
    to_data = ax.transData.inverted().transform
    gap = (min_gap or fontsize * 1.35) * fig.dpi / 72
    pts = sorted(pts, key=lambda p: p[1])
    ys = [to_disp((p[0], p[1]))[1] for p in pts]
    # forward then backward pass: keep order, enforce spacing, stay centred
    for i in range(1, len(ys)):
        ys[i] = max(ys[i], ys[i - 1] + gap)
    shift = (np.mean([to_disp((p[0], p[1]))[1] for p in pts]) - np.mean(ys))
    ys = [y + min(shift, 0) for y in ys]
    lo = to_disp((0, ax.get_ylim()[0]))[1] + gap / 2
    if ys[0] < lo:
        ys = [y + lo - ys[0] for y in ys]
    xc = to_disp((x_col, 0))[0]
    for (x, y, text), yd in zip(pts, ys):
        _, ty = to_data((xc, yd))
        ax.annotate(text, xy=(x, y), xytext=(x_col, ty), ha=ha, va="center",
                    fontsize=fontsize, color=color, zorder=6,
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="#888888",
                                    shrinkA=1.5, shrinkB=2.5))


def label_stair(ax, pts, extend, heights, fontsize=7.2, color="#222222"):
    """
    Label points from a block above or below them, one label per point, each
    directly over (or under) its point on a vertical leader line. `extend` is
    the side the text runs towards ('right' or 'left'); `heights` are the label
    rows from the one nearest the points outwards. Labels are ordered so that
    no leader line passes through another label.
    """
    if not pts:
        return
    pts = sorted(pts, key=lambda p: -p[0] if extend == "right" else p[0])
    heights = list(heights)[:len(pts)]
    above = heights[0] > max(p[1] for p in pts)
    used = []
    for (x, y, text), ly in zip(pts, heights):
        while any(abs(x - u) < 0.05 for u in used):     # shared x: nudge apart
            x += -0.07 if extend == "right" else 0.07
        used.append(x)
        ax.plot([x, x], [y + (0.06 if above else -0.06), ly], color="#9a9a9a", lw=0.5, zorder=2)
        ax.text(x, ly, " " + text if extend == "right" else text + " ",
                ha="left" if extend == "right" else "right",
                va="bottom" if above else "top", fontsize=fontsize, color=color, zorder=6)


def pca(mat, samples):
    comp = C.compositional(mat)
    X = StandardScaler().fit_transform(comp[samples].T.values)
    p = PCA(n_components=2).fit(X)
    return p.transform(X), p.explained_variance_ratio_ * 100, comp


# ------------------------------------------------------------------ Figure 3
def fig3B(cov_m):
    samples = [s for s in C.SOURCE + C.ENRICHMENT if s in cov_m.columns]
    sc, var, comp = pca(cov_m, samples)
    F, p = C.permanova(comp, {"Enrichment": C.ENRICHMENT, "Source": C.SOURCE})
    fig, ax = plt.subplots(figsize=(4.6, 3.7))
    ax.grid(color=GRID, lw=0.6, zorder=0)
    for grp, members, col in [("SOURCE", C.SOURCE, SRC), ("ENRICHMENT", C.ENRICHMENT, ENR)]:
        idx = [i for i, s in enumerate(samples) if s in members]
        ellipse(ax, sc[idx], col)
        ax.scatter(sc[idx, 0], sc[idx, 1], s=46, color=col, edgecolor="#333333",
                   lw=0.6, label=grp, zorder=3)
    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    ax.set_title("Functional community structure")
    ax.legend(loc="upper center", frameon=True, framealpha=0.9, edgecolor="#cccccc",
              handletextpad=0.2, borderpad=0.4)
    ax.text(0.98, 0.97, f"PERMANOVA $p$ = {p:.3f}", transform=ax.transAxes,
            ha="right", va="top", fontsize=7.5, color="#555555")
    fig.savefig(OUT / "fig3B.pdf"); fig.savefig(OUT / "fig3B.svg")
    pd.DataFrame(sc, index=samples, columns=["PC1", "PC2"]).to_csv(OUT / "fig3B.csv")
    plt.close(fig)
    return var, F, p


def fig3D(t):
    """Coverage-based against count-based log2 fold change, one point per module."""
    names = t.module_name.to_dict()
    d = t.dropna(subset=["cov_l2fc", "cnt_l2fc"])
    fig, ax = plt.subplots(figsize=(4.9, 4.2))
    ax.grid(color=GRID, lw=0.6, zorder=0)
    ax.axhline(0, color="#999999", lw=0.7, zorder=1)
    ax.axvline(0, color="#999999", lw=0.7, zorder=1)
    for cls in ["neither", "coverage-based only", "both metrics", "count-based only"]:
        s = d[d.niche_class == cls]
        big = cls in ("both metrics", "count-based only")
        ax.scatter(s.cov_l2fc, s.cnt_l2fc, s=34 if big else 14, color=CLS[cls],
                   edgecolor="white" if big else "none", lw=0.5, zorder=3 if big else 2,
                   label=cls)
    ax.set_xlim(d.cov_l2fc.min() - 1.6, d.cov_l2fc.max() + 0.3)
    ax.set_ylim(d.cnt_l2fc.min() - 0.3, d.cnt_l2fc.max() + 0.4)
    ax.set_xlabel("Coverage-based log$_2$FC (enrichment / source)")
    ax.set_ylabel("Count-based log$_2$FC")
    ax.set_title("Differentially abundant functions")

    co = d[d.niche_class == "count-based only"]
    label_stair(ax, [(r.cov_l2fc, r.cnt_l2fc, short(m, names)) for m, r in co.iterrows()],
                extend="left", heights=np.linspace(2.05, 4.3, len(co)),
                color=CLS["count-based only"])
    bo = d[d.niche_class == "both metrics"].sort_values("cov_l2fc", ascending=False).head(6)
    label_stair(ax, [(r.cov_l2fc, r.cnt_l2fc, short(m, names)) for m, r in bo.iterrows()],
                extend="left", heights=np.linspace(1.3, -1.75, len(bo)),
                color=CLS["both metrics"])
    ax.legend(loc="upper right", fontsize=7, frameon=True, framealpha=0.9,
              edgecolor="#cccccc", handletextpad=0.1, borderpad=0.4, labelspacing=0.3)
    fig.savefig(OUT / "fig3D.pdf"); fig.savefig(OUT / "fig3D.svg")
    d[["module_name", "cov_l2fc", "cov_padj", "cnt_l2fc", "cnt_padj", "niche_class"]].to_csv(OUT / "fig3D.csv")
    plt.close(fig)


def fig3D_volcano(t):
    names = t.module_name.to_dict()
    d = t.dropna(subset=["cov_l2fc"])
    y = -np.log10(d.cov_padj.clip(1e-5))
    fig, ax = plt.subplots(figsize=(4.9, 4.2))
    ax.grid(color=GRID, lw=0.6, zorder=0)
    for cls in ["neither", "coverage-based only", "both metrics", "count-based only"]:
        s = d.niche_class == cls
        big = cls in ("both metrics", "count-based only")
        ax.scatter(d.cov_l2fc[s], y[s], s=34 if big else 14, color=CLS[cls],
                   edgecolor="white" if big else "none", lw=0.5, zorder=3 if big else 2,
                   label=cls)
    ax.axhline(-np.log10(C.FDR), color="#e0605e", lw=1, ls="--", label="FDR 5%")
    for v in (-C.MIN_L2FC, C.MIN_L2FC):
        ax.axvline(v, color="#aaaaaa", lw=1, ls="--")
    ax.set_xlabel("Coverage-based log$_2$FC (enrichment / source)")
    ax.set_ylabel("$-$log$_{10}$ $q$")
    ax.set_title("Differentially abundant functions")
    co = d[d.niche_class == "count-based only"]
    label_column(ax, [(r.cov_l2fc, y[m], short(m, names)) for m, r in co.iterrows()],
                 x_col=ax.get_xlim()[0] + 0.1, ha="left", color=CLS["count-based only"])
    ax.legend(loc="upper left", fontsize=7, frameon=True, framealpha=0.9, edgecolor="#cccccc")
    fig.savefig(OUT / "fig3D_volcano.pdf"); fig.savefig(OUT / "fig3D_volcano.svg")
    plt.close(fig)


def fig3E(t):
    enr.PANEL_E_MODULES = list(E_LABEL)
    dol = enr.division_of_labour()
    dol.to_csv(OUT / "fig3E.csv", index=False)
    piv = dol.pivot_table(index="module", columns="sample", values="cyano_percent", dropna=False)
    piv = piv.reindex(index=list(E_LABEL), columns=E_ORDER)

    rows = []
    for gi, (g, members) in enumerate(E_GROUPS):
        for m, lab in members:
            rows.append((gi, g, m, lab))
    n = len(rows) + len(E_GROUPS) - 1
    fig = plt.figure(figsize=(7.6, 0.27 * n + 1.5))
    ax = fig.add_axes([0.36, 0.12, 0.47, 0.83])
    greens, blues = mpl.colormaps["Greens"], mpl.colormaps["Blues"]
    norm = Normalize(0, 100)
    y, prev, group_y, cell_rows = n - 1, None, {}, []
    for gi, g, m, lab in rows:
        if prev is not None and gi != prev:
            ax.plot([-0.5, len(E_ORDER) - 0.5], [y + 0.5, y + 0.5], color="black", lw=2.2,
                    clip_on=False)
            y -= 1
        prev = gi
        group_y.setdefault(g, []).append(y)
        cell_rows.append((g, m, y))
        for j, smp in enumerate(E_ORDER):
            v = piv.loc[m, smp]
            x0, y0 = j - 0.5, y - 0.45
            if np.isnan(v):
                ax.add_patch(mpatches.Rectangle((x0, y0), 1, 0.9, facecolor="#f2f2f2",
                                                edgecolor="#bbbbbb", lw=0.5))
                continue
            ax.add_patch(mpatches.Polygon([[x0, y0], [x0, y0 + 0.9], [x0 + 1, y0 + 0.9]],
                         facecolor=greens(0.08 + 0.85 * norm(v)), edgecolor="none"))
            ax.add_patch(mpatches.Polygon([[x0, y0], [x0 + 1, y0], [x0 + 1, y0 + 0.9]],
                         facecolor=blues(0.08 + 0.85 * norm(100 - v)), edgecolor="none"))
            ax.add_patch(mpatches.Rectangle((x0, y0), 1, 0.9, facecolor="none",
                                            edgecolor="#555555", lw=0.5))
            ax.text(x0 + 0.06, y0 + 0.86, f"{v:.0f}", ha="left", va="top", fontsize=6.3,
                    color="white" if v > 60 else "#222222")
            ax.text(x0 + 0.94, y0 + 0.04, f"{100 - v:.0f}", ha="right", va="bottom",
                    fontsize=6.3, color="white" if v < 40 else "#222222")
        ax.text(-0.95, y, lab, ha="right", va="center", fontsize=8)
        cls = t.niche_class.get(m, "neither")
        if cls in ("both metrics", "count-based only"):
            ax.scatter([-0.72], [y], s=22, color=CLS[cls], clip_on=False, zorder=5)
        elif cls == "coverage-based only":
            ax.scatter([-0.72], [y], s=22, facecolor="white", edgecolor=CLS["both metrics"],
                       lw=0.8, clip_on=False, zorder=5)
        y -= 1
    for g, ys in group_y.items():
        ax.text(-10.3, np.mean(ys), g, ha="center", va="center", fontsize=9,
                fontweight="bold")
    # how much of the cyanobacterial community the MAGs capture, per replicate
    rep = C.mag_representation().reindex(E_ORDER)
    ys = -1.0
    for j, smp in enumerate(E_ORDER):
        v = rep[smp]
        ax.add_patch(mpatches.Rectangle((j - 0.5, ys - 0.3), 1, 0.6,
                     facecolor=mpl.colormaps["Greys"](0.05 + 0.5 * v / 100),
                     edgecolor="#555555", lw=0.5))
        ax.text(j, ys, f"{v:.0f}", ha="center", va="center", fontsize=6.5)
        if v < 25:
            for _, _, yy in cell_rows:
                ax.add_patch(mpatches.Rectangle((j - 0.5, yy - 0.45), 1, 0.9, facecolor="white",
                             alpha=0.6, edgecolor="none", zorder=4))
        ax.text(j, -1.55, smp, ha="center", va="top", fontsize=8)
    ax.text(-0.95, ys, "Cyanobacteria in MAGs (%)", ha="right", va="center", fontsize=7.5,
            style="italic", color="#444444")
    ax.text(len(E_ORDER) / 2 - 0.5, -2.3, "ENRICHMENT Sample", ha="center", va="top", fontsize=10)
    ax.set_xlim(-0.5, len(E_ORDER) - 0.5)
    ax.set_ylim(-1.4, n - 0.4)
    ax.axis("off")
    ax.set_title("Functional partitioning within communities", pad=6)

    for k, (cmap, lab) in enumerate([(greens, "Cyanobacteria contribution (%)"),
                                     (blues, "Other bacteria contribution (%)")]):
        cax = fig.add_axes([0.875 + 0.05 * k, 0.25, 0.022, 0.55])
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=mpl.colors.ListedColormap(
            cmap(np.linspace(0.08, 0.93, 256)))), cax=cax)
        cb.set_label(lab, fontsize=8)
        cb.ax.tick_params(labelsize=7)
        if k == 0:
            cb.ax.yaxis.set_label_position("left")
            cb.ax.tick_params(labelright=False, right=False)
    handles = [plt.Line2D([], [], marker="o", ls="", color=CLS["both metrics"], label="both metrics"),
               plt.Line2D([], [], marker="o", ls="", color=CLS["count-based only"], label="count-based only"),
               plt.Line2D([], [], marker="o", ls="", markerfacecolor="white",
                          markeredgecolor=CLS["both metrics"], label="coverage-based only"),
               mpatches.Patch(facecolor="#f2f2f2", edgecolor="#bbbbbb", label="no encoding MAG"),
               mpatches.Patch(facecolor="#e9f3ec", edgecolor="#bbbbbb",
                              label="faded: <25% of cyanobacteria in MAGs")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.02, 0.0), ncol=5,
               fontsize=6.8, handletextpad=0.2, columnspacing=1.0)
    fig.savefig(OUT / "fig3E.pdf"); fig.savefig(OUT / "fig3E.svg")
    plt.close(fig)
    return piv


# ------------------------------------------------------------------ Figure 4
def fig4C(cov_m):
    groups = {g: [s for s in med.GROUPS[g] if s in cov_m.columns] for g in ("V", "NV", "NVNN")}
    samples = sum(groups.values(), [])
    sc, var, comp = pca(cov_m, samples)
    F, p = C.permanova(comp, groups)
    fig, ax = plt.subplots(figsize=(4.6, 4.0))
    ax.grid(color=GRID, lw=0.6, zorder=0)
    lab = {"V": "BG11 + vitamins + NaNO$_3$ (V)", "NV": "BG11 $-$ vitamins + NaNO$_3$ (NV)",
           "NVNN": "BG11 $-$ vitamins $-$ NaNO$_3$ (NVNN)"}
    mk = {"V": "o", "NV": "s", "NVNN": "^"}
    for g, members in groups.items():
        idx = [samples.index(s) for s in members]
        ellipse(ax, sc[idx], MEDIA[g], dashed=True)
        ax.scatter(sc[idx, 0], sc[idx, 1], s=48, marker=mk[g], color=MEDIA[g],
                   edgecolor="#333333", lw=0.5, label=lab[g], zorder=3)
    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    ax.set_title("Functional community structure")
    ax.legend(loc="upper right", fontsize=7, frameon=True, framealpha=0.9, edgecolor="#cccccc")
    ax.text(0.98, 0.74, f"PERMANOVA $p$ = {p:.3f}", transform=ax.transAxes, fontsize=7.5,
            color="#555555", ha="right")
    fig.savefig(OUT / "fig4C.pdf"); fig.savefig(OUT / "fig4C.svg")
    pd.DataFrame(sc, index=samples, columns=["PC1", "PC2"]).to_csv(OUT / "fig4C.csv")
    plt.close(fig)
    return var, F, p


def fig4D(res, names):
    """NVNN against the pooled nitrate-replete cultures, count-based metric."""
    d = res["nitrate"]["cnt"].join(names).dropna(subset=["l2fc"])
    nm = d.module_name.to_dict()
    y = -np.log10(d.p)
    sig = d.sig
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    ax.grid(color=GRID, lw=0.6, zorder=0)
    col = np.where(~sig, "#b9bcc0", np.where(d.l2fc > 0, MEDIA["NVNN"], MEDIA["NV"]))
    ax.scatter(d.l2fc, y, s=np.where(sig, 36, 13), c=col,
               edgecolor=np.where(sig, "#333333", "none"), lw=0.5, zorder=3)
    for v in (-C.MIN_L2FC, C.MIN_L2FC):
        ax.axvline(v, color="#aaaaaa", lw=1, ls="--")
    ax.set_xlim(min(d.l2fc.min() - 0.6, -6.4), max(d.l2fc.max() + 0.6, 6.4))
    top = y.max() + 1.6
    ax.set_ylim(-0.1, top)
    if sig.any():
        cut = -np.log10(d.p[d.padj < C.FDR].max())
        ax.axhline(cut, color="#e0605e", lw=1, ls="--")
        ax.text(ax.get_xlim()[1] - 0.1, cut - 0.05, "FDR 5%", fontsize=7, color="#e0605e",
                ha="right", va="top")
    s = d[sig]
    rows = np.linspace(y.max() + 0.12, top - 0.2, 7)
    label_stair(ax, [(r.l2fc, y[m], short(m, nm)) for m, r in s[s.l2fc > 0].iterrows()],
                extend="right", heights=rows)
    label_stair(ax, [(r.l2fc, y[m], short(m, nm)) for m, r in s[s.l2fc < 0].iterrows()],
                extend="left", heights=rows)
    ax.set_xlabel("Count-based log$_2$FC  [NVNN / (V + NV)]")
    ax.set_ylabel("$-$log$_{10}$ $p$")
    ax.set_title("Differentially abundant functions")
    ax.annotate("", xy=(ax.get_xlim()[0] + 0.2, -0.2), xytext=(-0.15, -0.2),
                xycoords=("data", "axes fraction"), textcoords=("data", "axes fraction"),
                arrowprops=dict(arrowstyle="-|>", color=MEDIA["NV"], lw=2.5))
    ax.annotate("", xy=(ax.get_xlim()[1] - 0.2, -0.2), xytext=(0.15, -0.2),
                xycoords=("data", "axes fraction"), textcoords=("data", "axes fraction"),
                arrowprops=dict(arrowstyle="-|>", color=MEDIA["NVNN"], lw=2.5))
    x0 = (0 - ax.get_xlim()[0]) / np.diff(ax.get_xlim())[0]
    ax.text(x0 / 2, -0.175, "higher with nitrate", transform=ax.transAxes, ha="center",
            va="bottom", fontsize=7.5)
    ax.text((1 + x0) / 2, -0.175, "higher without nitrate", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=7.5)
    fig.savefig(OUT / "fig4D.pdf"); fig.savefig(OUT / "fig4D.svg")
    d[["module_name", "l2fc", "p", "padj", "sig", "separated"]].to_csv(OUT / "fig4D.csv")
    plt.close(fig)


def run_panels():
    pub_style()
    t, cov_m, _ = enr.two_metric("rpl19")
    v3, F3, p3 = fig3B(cov_m)
    fig3D(t)
    fig3D_volcano(t)
    piv = fig3E(t)
    res, mcov, _, names = med.contrasts("rps16")
    v4, F4, p4 = fig4C(mcov)
    fig4D(res, names[["module_name"]])
    print(f"3B: PC1 {v3[0]:.1f}%, PC2 {v3[1]:.1f}%, PERMANOVA F={F3:.2f} p={p3:.3f}")
    print(f"4C: PC1 {v4[0]:.1f}%, PC2 {v4[1]:.1f}%, PERMANOVA F={F4:.2f} p={p4:.3f}")
    print("3E cyanobacterial %:")
    print(piv.round(0).to_string())
    print(f"wrote {OUT}")

# ============================================================================
# Composite Figures 3 and 4
# ----------------------------------------------------------------------------
# The full figures as they go into the paper. Data panels are drawn here;
# the schematics (3A, 3C, 4A) and the anvi'o EcoPhylo heatmap (4B) are vector
# crops of the Affinity composites, kept in data/schematics/ and placed with
# PyMuPDF so that nothing is rasterised.
#     outputs/FIG3.pdf, outputs/FIG4.pdf
# Layout is given in points from the top-left corner of the page.
# ============================================================================
from scipy.stats import chi2  # noqa: E402

SCHEM = C.DATA / "schematics"
CI95 = float(np.sqrt(chi2.ppf(0.95, 2)))     # 95% confidence ellipse, in SDs
PHOTO = {"M00161", "M00163", "M00597", "M00598", "M00162"}
NITRO = {"M00175", "M00531", "M00530", "M00529", "M00528", "M00804", "M00973", "M00615"}
GREEN, MAGENTA = "#3fa34d", "#d16bd6"

# Fig 3E rows, as in the caption: 12 modules in four groups. Menaquinone (M00116) is left
# out: it rose on the coverage-based metric only and fell on the count-based metric.
E3_GROUPS = [
    ("Photo-\nsynthesis", [("M00161", "Photosystem II"), ("M00163", "Photosystem I"),
                           ("M00597", "Anoxygenic photosystem II")]),
    ("Nitrate\nmetabolism", [("M00529", "Denitrification"),
                             ("M00531", "Assimilatory nitrate reduction"),
                             ("M00615", "Nitrate assimilation")]),
    ("Vitamins and\nco-factors", [("M00122", "Cobalamin biosynthesis"),
                                  ("M00123", "Biotin biosynthesis"),
                                  ("M00899", "Thiamine salvage pathway")]),
    ("Amino acid\nmetabolism", [("M00118", "Glutathione biosynthesis"),
                                ("M00133", "Polyamine biosynthesis"),
                                ("M00970", "Proline degradation")]),
]


class Page:
    """A matplotlib figure addressed in points from the top-left corner."""

    def __init__(self, w, h):
        self.w, self.h = w, h
        self.fig = plt.figure(figsize=(w / 72, h / 72))

    def ax(self, x0, y0, x1, y1, **kw):
        return self.fig.add_axes([x0 / self.w, 1 - y1 / self.h, (x1 - x0) / self.w,
                                  (y1 - y0) / self.h], **kw)

    def text(self, x, y, s, **kw):
        kw.setdefault("va", "top")
        self.fig.text(x / self.w, 1 - y / self.h, s, **kw)

    def letter(self, x, y, s):
        self.text(x, y, s, fontsize=10, fontweight="bold", ha="left")

    def title(self, x, y, s, size=8):
        self.text(x, y, s, fontsize=size, fontweight="bold", ha="center")

    def save(self, path, crops):
        """Save, then place each schematic crop {name: (x0, y0, x1, y1)} as vector."""
        import pymupdf
        tmp = path.with_suffix(".tmp.pdf")
        C.polish(self.fig)
        self.fig.savefig(tmp, bbox_inches=None)
        plt.close(self.fig)
        doc = pymupdf.open(tmp)
        for name, r in crops.items():
            r, clip = (r if isinstance(r[0], tuple) else (r, None))
            src = pymupdf.open(SCHEM / f"{name}.pdf")
            doc[0].show_pdf_page(pymupdf.Rect(*r), src, 0, keep_proportion=True,
                                 clip=pymupdf.Rect(*clip) if clip else None)
        doc.save(path, garbage=4, deflate=True)
        tmp.unlink()
        C.crop_pdf(path)
        doc = pymupdf.open(path)
        doc[0].get_pixmap(dpi=200).save(path.with_suffix(".png"))


def small_style():
    pub_style()
    mpl.rcParams.update({"font.size": 7.5, "axes.labelsize": 8, "xtick.labelsize": 7,
                         "ytick.labelsize": 7, "legend.fontsize": 6.8, "axes.linewidth": 0.7,
                         "xtick.major.size": 3, "ytick.major.size": 3,
                         "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                         "xtick.minor.size": 1.8, "ytick.minor.size": 1.8,
                         "xtick.minor.width": 0.45, "ytick.minor.width": 0.45,
                         "xtick.direction": "out", "ytick.direction": "out",
                         "xtick.major.pad": 2, "ytick.major.pad": 2,
                         "savefig.bbox": "standard"})


def score_ticks(ax, labelsize=5.6):
    """Ordination axes keep their scale: three or four ticks at round values on both axes."""
    from matplotlib.ticker import MaxNLocator, FormatStrFormatter
    for a in (ax.xaxis, ax.yaxis):
        a.set_major_locator(MaxNLocator(nbins=3, steps=[1, 2, 5, 10]))
        a.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v:g}".replace("-", "\u2212")))
    ax.tick_params(labelsize=labelsize, pad=1.5)


def similarity_ticks(ax, labelsize=6):
    """x axis of the pair strips: 0, 0.5, 1 labelled, quarters as minor ticks."""
    from matplotlib.ticker import MultipleLocator
    ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"])
    ax.xaxis.set_minor_locator(MultipleLocator(0.25))
    ax.tick_params(axis="x", labelsize=labelsize)


SIMILARITY_LABEL = "Similarity of a pair of samples (1 \u2212 Bray\u2013Curtis)"


def group_var(sc, idx):
    """Total variance of a group's scores in PC1 + PC2."""
    return float(sc[idx].var(axis=0, ddof=1).sum())


# ------------------------------------------------------------------ Figure 3
def f3_pca(ax, cov_m):
    """Bray-Curtis PCoA of the coverage-based module profiles (each sample divided by its
    total), as for taxonomy in Fig 2C-D. Group spread is reported as the mean distance to
    the group centroid in the full PCoA space, with PERMDISP: the 2-D spread of a PCA on
    z-scored modules (the earlier inset) depends on the scaling and is not a dispersion."""
    groups = {"SOURCE": C.SOURCE, "ENRICHMENT": C.ENRICHMENT}
    r = C.bray_pcoa(C.compositional(cov_m), groups)
    sc = r["scores"]
    ax.grid(color=GRID, lw=0.5, zorder=0)
    for grp, col in (("SOURCE", SRC), ("ENRICHMENT", ENR)):
        xy = sc.loc[[s for s in groups[grp] if s in sc.index]].values
        ellipse(ax, xy, col, n_std=CI95)
        ax.scatter(xy[:, 0], xy[:, 1], s=18, color=col, edgecolor="#333333", lw=0.4,
                   label=grp, zorder=3)
    ax.set_xlabel(f"PCo1 ({r['explained'][0]:.1f}%)")
    ax.set_ylabel(f"PCo2 ({r['explained'][1]:.1f}%)")
    ax.legend(loc="lower left", frameon=True, framealpha=0.9, edgecolor="#cccccc",
              handletextpad=0.1, borderpad=0.3)
    cd = r["centroid_dist"]
    ax.text(0.97, 0.03, f"PERMANOVA {fmt_p(r['p'])}\n"
            f"dispersion S {cd['SOURCE']:.2f}, E {cd['ENRICHMENT']:.2f}\n"
            f"PERMDISP {fmt_p(r['disp_p'])}", transform=ax.transAxes, ha="right",
            va="bottom", fontsize=6.0, color="#555555", linespacing=1.25)
    ax.margins(0.12)
    return dict(explained=r["explained"], F=r["F"], p=r["p"], disp_p=r["disp_p"], **cd)


def f3_volcano(ax, t):
    d = t.dropna(subset=["cov_l2fc"]).copy()
    y = -np.log10(d.cov_p)
    sig = (d.cov_padj < C.FDR) & (d.cov_l2fc.abs() > C.MIN_L2FC)
    cat = np.where(d.index.isin(PHOTO), "photo", np.where(d.index.isin(NITRO), "nitro", "other"))
    ax.grid(color=GRID, lw=0.5, zorder=0)
    ax.scatter(d.cov_l2fc[~sig], y[~sig], s=4, color="#9a9a9a", lw=0, zorder=2)
    col = {"photo": GREEN, "nitro": MAGENTA, "other": "#6b6b6b"}
    lab = {"photo": "Photosynthesis", "nitro": "Nitrogen metabolism", "other": "Other"}
    for c in ("other", "photo", "nitro"):
        s = sig & (cat == c)
        ax.scatter(d.cov_l2fc[s], y[s], s=13 if c == "other" else 20, color=col[c],
                   edgecolor="#333333", lw=0.35, zorder=3 if c == "other" else 4, label=lab[c])
    cut = -np.log10(d.cov_p[d.cov_padj < C.FDR].max())
    ax.axhline(cut, color="#e0605e", lw=0.8, ls="--", label="FDR 5%")
    for v in (-C.MIN_L2FC, C.MIN_L2FC):
        ax.axvline(v, color="#aaaaaa", lw=0.7, ls="--")
    ax.set_xlabel("Coverage-based log$_2$FC (enrichment / source)")
    ax.set_ylabel("$-$log$_{10}$ $p$")
    # the exact test cannot go below p = 2 / C(16, 8): every module whose 8 enrichments all lie
    # above (or all below) the 8 sources sits on this line
    floor = y.max()
    n_floor = int((y >= floor - 1e-9).sum())
    ax.axhline(floor, color="#bbbbbb", lw=0.6, ls=":", zorder=1)
    ax.text(0.02, floor - 0.1, f"complete separation,\n8 vs 8 ({n_floor} modules)", transform=ax.get_yaxis_transform(),
            ha="left", va="top", fontsize=5.4, color="#555555")
    ax.legend(loc="lower right", frameon=True, framealpha=0.9, edgecolor="#cccccc",
              handletextpad=0.1, borderpad=0.3, labelspacing=0.25)
    # key functions, labelled from rows above the points on vertical leaders
    key = d.index[sig & (d.index.isin(PHOTO | NITRO) | d.index.isin(list(E_LABEL)))]
    key = key[:12]
    heights = [floor + 0.35 + 0.27 * i for i in range(len(key))]
    pts = [(d.loc[m, "cov_l2fc"], y[m], re.sub(r"\s*\(.*?\)", "", SHORT.get(m, m))) for m in key]
    label_stair(ax, pts, extend="left", heights=heights, fontsize=5.4)
    col = {re.sub(r"\s*\(.*?\)", "", SHORT.get(m, m)): GREEN if m in PHOTO else MAGENTA if m in NITRO else "#333333"
           for m in key}
    for t_ in ax.texts:
        if t_.get_text().strip() in col:
            t_.set_color(col[t_.get_text().strip()])
    ax.set_ylim(-0.1, (heights[-1] if heights else floor) + 0.35)
    ax.set_yticks(range(0, int(floor) + 1))
    return dict(n_sig=int(sig.sum()), n_floor=n_floor)


F3_KEY_BOTH = ["M00161", "M00163", "M00162", "M00615", "M00531"]


def f3_two_metrics(ax, t):
    """Fig 3D: coverage-based against count-based log2FC (enrichment / source), one point per
    module, coloured by two-metric class; count-based-only modules and key shared ones labelled."""
    d = t.dropna(subset=["cov_l2fc", "cnt_l2fc"])
    names = t.module_name.to_dict()
    ax.axhline(0, color="#999999", lw=0.6); ax.axvline(0, color="#999999", lw=0.6)
    order = [("neither", 5), ("coverage-based only", 9), ("both metrics", 16), ("count-based only", 22)]
    for cls, size in order:
        sub = d[d.niche_class == cls]
        ax.scatter(sub.cov_l2fc, sub.cnt_l2fc, s=size, color=CLS[cls], lw=0.35,
                   edgecolor="#333333" if cls != "neither" else "none", zorder=2 + order.index((cls, size)),
                   label=f"{cls} ({int((t.niche_class == cls).sum())})")
    lo = d.cov_l2fc.min() - 3.2
    ax.set_xlim(lo, d.cov_l2fc.max() + 0.4)
    ax.set_ylim(d.cnt_l2fc.min() - 0.3, d.cnt_l2fc.max() + 0.5)
    lab = d.index[(d.niche_class == "count-based only") | (d.index.isin(F3_KEY_BOTH) & (d.niche_class == "both metrics"))]
    pts = [(d.loc[m, "cov_l2fc"], d.loc[m, "cnt_l2fc"], short(m, names, False)) for m in lab]
    label_column(ax, pts, x_col=lo + 0.1, ha="left", fontsize=5.4)
    col = {short(m, names, False): CLS[d.loc[m, "niche_class"]] for m in lab}
    for t_ in ax.texts:
        if t_.get_text() in col:
            t_.set_color(col[t_.get_text()])
    ax.set_xlabel("Coverage-based log$_2$FC")
    ax.set_ylabel("Count-based log$_2$FC")
    ax.text(0.98, 0.98, "enrichment / source", transform=ax.transAxes, ha="right", va="top", fontsize=5.8, color="#444444")
    ax.legend(loc="lower right", fontsize=5.4, frameon=True, framealpha=0.9, edgecolor="#cccccc",
              handletextpad=0.1, borderpad=0.3, labelspacing=0.25)
    return d.niche_class.value_counts().to_dict()


F3_CYANO, F3_HETERO, F3_MID = "#1b7837", "#2166ac", "#f2f2f2"


def f3_partition(page, x0, y0, x1, y1, t, cbar_x=None):
    """Fig 3E: share of encoding-MAG coverage in Cyanobacteriota, per replicate, on one
    diverging scale (green, cyanobacterial; blue, heterotrophic; the two sum to 100).
    Hatched cell: no detected MAG encodes the module."""
    from matplotlib.colors import LinearSegmentedColormap
    enr.PANEL_E_MODULES = [m for _, rows in E3_GROUPS for m, _ in rows]
    dol = enr.division_of_labour()
    piv = dol.pivot_table(index="module", columns="sample", values="cyano_percent", dropna=False)
    piv = piv.reindex(index=enr.PANEL_E_MODULES, columns=E_ORDER)
    rep = C.mag_representation().reindex(E_ORDER)
    cmap = LinearSegmentedColormap.from_list("cyano_het", [F3_HETERO, F3_MID, F3_CYANO])
    norm = Normalize(0, 100)
    gap = 0.45
    n_rows = len(enr.PANEL_E_MODULES) + gap * (len(E3_GROUPS) - 1)
    ax = page.ax(x0, y0, x1, y1)
    y = n_rows - 1
    for gi, (g, members) in enumerate(E3_GROUPS):
        if gi:
            y -= gap
        ys = []
        for m, lab in members:
            ys.append(y)
            for j, smp in enumerate(E_ORDER):
                v = piv.loc[m, smp]
                if np.isnan(v):
                    ax.add_patch(mpatches.Rectangle((j - 0.5, y - 0.5), 1, 1, fc="white", ec="#d0d0d0",
                                                    lw=0.4, hatch="////"))
                    continue
                ax.add_patch(mpatches.Rectangle((j - 0.5, y - 0.5), 1, 1, fc=cmap(norm(v)), ec="white", lw=0.6))
                ax.text(j, y, f"{v:.0f}", ha="center", va="center", fontsize=5.4,
                        color="white" if abs(v - 50) > 30 else "#222222")
            ax.text(-0.65, y, lab, ha="right", va="center", fontsize=6.4)
            y -= 1
        gx = (8 - x0) / ((x1 - x0) / len(E_ORDER)) - 0.5      # group labels at 8 pt from the page edge
        ax.text(gx, np.mean(ys), g, ha="left", va="center", fontsize=6.6, fontweight="bold")
    for j, smp in enumerate(E_ORDER):
        ax.text(j, -0.65, smp, ha="center", va="top", fontsize=6.2)
        ax.text(j, -1.55, f"{rep[smp]:.0f}", ha="center", va="top", fontsize=5.6, color="#666666")
    ax.text(-0.65, -1.55, "cyanobacteria captured by MAGs (%)", ha="right", va="top",
            fontsize=5.6, color="#666666", style="italic")
    ax.set_xlim(-0.5, len(E_ORDER) - 0.5)
    ax.set_ylim(-0.5, n_rows - 0.5)
    ax.axis("off")
    cx = cbar_x if cbar_x is not None else x1 + 14
    cax = page.ax(cx, y0 + 6, cx + 4, y1 - 30)
    cb = page.fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, ticks=[0, 50, 100])
    cb.ax.tick_params(labelsize=5.8, length=1.5, width=0.4, pad=1)
    cb.outline.set_linewidth(0)
    cb.ax.text(0.5, 1.03, "cyano-\nbacteria", transform=cb.ax.transAxes, ha="center", va="bottom",
               fontsize=5.8, color=F3_CYANO)
    cb.ax.text(0.5, -0.03, "hetero-\ntrophs", transform=cb.ax.transAxes, ha="center", va="top",
               fontsize=5.8, color=F3_HETERO)
    cb.set_label("% of encoding-MAG coverage", fontsize=6, labelpad=3)
    cb.ax.yaxis.set_label_position("left")
    kx, ky = cx - 10, y1 - 4                       # key for the hatched cells
    page.fig.patches.append(mpatches.Rectangle((kx / page.w, 1 - (ky + 7) / page.h), 10 / page.w, 7 / page.h,
                            transform=page.fig.transFigure, fc="white", ec="#b0b0b0", lw=0.4, hatch="////"))
    page.text(kx + 13, ky + 0.5, "no detected MAG\nencodes the module", fontsize=5.6, ha="left")
    return piv


F3_BETWEEN = "#9a9a9a"


def f3_convergence(ax):
    """Fig 3B bottom: coverage shared by pairs of samples (1 - Bray-Curtis) within the
    enrichments, within the source, and between an enrichment and a source sample, for
    populations and for function on both metrics. The between-group pairs are the reference
    that shows the enrichments left the source's profile in the same direction. Dots, pairs;
    black bar, median (value at right); shaded bar, interquartile range."""
    import taxonomy as tx
    from scipy.spatial.distance import pdist, cdist
    lt = tx.level_tables()
    rows = [("Populations\n(RpL19)", "population"), ("Function,\ncoverage", "function, coverage"),
            ("Function,\ncount", "function, count")]
    rng = np.random.default_rng(3)
    out = {}
    for i, (lab, key) in enumerate(rows):
        y0 = len(rows) - 1 - i
        m = lt[key]
        sets = [(0.27, "Enrichment", 1 - pdist(m[C.ENRICHMENT].T.values, "braycurtis")),
                (0.0, "Source", 1 - pdist(m[C.SOURCE].T.values, "braycurtis")),
                (-0.27, "Between", (1 - cdist(m[C.ENRICHMENT].T.values, m[C.SOURCE].T.values, "braycurtis")).ravel())]
        for dy, g, v in sets:
            y = y0 + dy
            col = F3_BETWEEN if g == "Between" else C.GROUP_COLOR[g]
            q1, q3 = np.percentile(v, [25, 75])
            med = np.median(v)
            ax.plot([q1, q3], [y, y], color=col, lw=3.4, alpha=0.3, solid_capstyle="butt", zorder=2)
            ax.scatter(v, y + rng.uniform(-0.06, 0.06, len(v)), s=3.5, lw=0.25, zorder=3,
                       marker="D" if g == "Between" else F2_MARKER[g], color=col,
                       edgecolor="#333333" if g == "Enrichment" else "white")
            ax.plot([med, med], [y - 0.1, y + 0.1], color="#222222", lw=1.0, zorder=4)
            ax.text(1.03, y, f"{med:.2f}", va="center", fontsize=5.6, color="#333333",
                    transform=ax.get_yaxis_transform())
            out[(key, g)] = med
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, len(rows) - 0.5)
    similarity_ticks(ax)
    ax.set_yticks(range(len(rows))[::-1], [r[0] for r in rows], fontsize=6)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel(SIMILARITY_LABEL, fontsize=6.5, labelpad=2)
    ax.grid(axis="x", color=GRID, lw=0.5, zorder=0)
    ax.spines[["top", "right", "left"]].set_visible(False)
    return out


def compose_fig3():
    """Fig 3: A workflow and C conditions (schematics); B function separated from the source
    (PCoA, PERMANOVA R2) yet shared between enrichments about as much as between source
    samples, unlike populations (pair strips, as Fig 2C); D two-metric classification;
    E who encodes selected modules in each enrichment."""
    import taxonomy as tx
    small_style()
    t, cov_m, _ = enr.two_metric("rpl19")
    W, H = 510, 470
    pg = Page(W, H)
    pg.title(134, 0, "Functional annotation workflow")
    pg.title(389, 0, "Enrichment conditions")
    pg.letter(4, 70, "B"); pg.title(138, 70, "Separation and convergence of function")
    r = f2_rank_pcoa(pg.ax(102, 98, 216, 150), tx.level_tables()["function, coverage"],
                     "Function (coverage-based)")
    conv = f3_convergence(pg.ax(102, 176, 216, 252))
    hs = [plt.Line2D([], [], ls="", marker=F2_MARKER[g], ms=4, color=C.GROUP_COLOR[g],
                     mec="#333333" if g == "Enrichment" else "white", mew=0.4, label=g)
          for g in ("Source", "Enrichment")]
    hs.append(plt.Line2D([], [], ls="", marker="D", ms=3.5, color=F3_BETWEEN, mec="white", mew=0.4,
                         label="Enrichment vs\nsource (pairs)"))
    pg.fig.legend(handles=hs, loc="upper left", ncol=1, fontsize=6, frameon=False,
                  bbox_to_anchor=(4 / W, 1 - 100 / H), handletextpad=0.2, labelspacing=0.5)
    pg.letter(266, 70, "D"); pg.title(390, 70, "Functions selected in the enrichments")
    axd = pg.ax(308, 88, 503, 250)
    v = f3_two_metrics(axd, t)
    axd.xaxis.label.set_size(6.5); axd.yaxis.label.set_size(6.5)
    axd.tick_params(labelsize=6)
    pg.letter(4, 276, "E"); pg.title(W / 2, 276, "Partitioning of selected functions between cyanobacteria and heterotrophs")
    f3_partition(pg, 172, 294, 420, 444, t, cbar_x=438)
    pg.save(C.OUT / "FIG3.pdf", {"fig3A_workflow": ((0, 12, 268, 62), (0, 9.5, 243, 58.5)),
                                 "fig3C_conditions": ((268, 12, 510, 62), (0, 9.5, 213, 58.5))})
    return dict(R2=r["R2"], p=r["p"], disp_p=r["disp_p"],
                **{f"{k[0]}|{k[1]}": round(x, 3) for k, x in conv.items()}), v


# ------------------------------------------------------------------ Figure 4
MEDIA_LAB = {"V": "V (+ vitamins, + nitrate)", "NV": "NV ($-$ vitamins, + nitrate)",
             "NVNN": "NVNN ($-$ vitamins, $-$ nitrate)"}
MEDIA_MK = {"V": "o", "NV": "s", "NVNN": "^"}

def media_groups(cols):
    return {g: [s for s in med.GROUPS[g] if s in cols] for g in ("V", "NV", "NVNN")}


def f4_ordination(ax, mat, title, show_legend=False):
    """PCA of a feature x sample matrix (z-scored), 95% ellipses, PERMANOVA."""
    groups = media_groups(mat.columns)
    samples = sum(groups.values(), [])
    X = StandardScaler().fit_transform(mat[samples].T.values)
    p = PCA(n_components=2).fit(X)
    sc, var = p.transform(X), p.explained_variance_ratio_ * 100
    F, pv = C.permanova(mat, groups)
    ax.grid(color=GRID, lw=0.5, zorder=0)
    for g, members in groups.items():
        idx = [samples.index(s) for s in members]
        ellipse(ax, sc[idx], MEDIA[g], dashed=True, n_std=CI95)
        ax.scatter(sc[idx, 0], sc[idx, 1], s=16, marker=MEDIA_MK[g], color=MEDIA[g],
                   edgecolor="#333333", lw=0.35, label=MEDIA_LAB[g], zorder=3)
    ax.set_xlabel(f"PC1 ({var[0]:.0f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.0f}%)")
    ax.set_title(title, fontsize=7, fontweight="bold", pad=3)
    ax.text(0.97, 0.03, f"PERMANOVA\n$F$ = {F:.2f}, $p$ = {pv:.3f}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=5.5, color="#444444", linespacing=1.2)
    ax.margins(0.15)
    score_ticks(ax)
    return dict(F=F, p=pv, var=var)


def rps16_matrix():
    e = pd.read_csv(C.DATA / "ecophylo_coverage.csv")
    e = e[e.run == "rps16_media_endpoints"]
    m = e.pivot_table(index="cluster", columns="sample", values="mean_coverage_Q2Q3").fillna(0)
    m = m.div(m.sum(axis=0), axis=1)
    return m[m.var(axis=1) > 0]


def f4_two_metrics(ax, res, names):
    """NVNN vs nitrate-replete, one point per module, on both metrics.

    Highlighted: modules whose NVNN replicates are completely separated from all
    nine nitrate-replete ones (the strongest result an exact test can give at
    n = 4 vs 9, p = 0.0014). FDR is not used to split the metrics here: with so
    few replicates the BH q of a perfectly separated module depends mostly on how
    many other modules are also perfectly separated.
    """
    cov, cnt = res["nitrate"]["cov"], res["nitrate"]["cnt"]
    d = cov[["l2fc", "separated"]].join(cnt[["l2fc", "separated"]], lsuffix="_cov", rsuffix="_cnt",
                                        how="inner").dropna(subset=["l2fc_cov", "l2fc_cnt"])
    sc, sv = d.separated_cnt.astype(bool), d.separated_cov.astype(bool)
    ax.axhline(0, color="#999999", lw=0.6); ax.axvline(0, color="#999999", lw=0.6)
    lim = max(d.l2fc_cov.abs().max(), d.l2fc_cnt.abs().max()) * 1.08
    ax.plot([-lim, lim], [-lim, lim], color="#cccccc", lw=0.6, ls=":", zorder=1)
    ax.scatter(d.l2fc_cov[~(sc | sv)], d.l2fc_cnt[~(sc | sv)], s=5, color="#b9bcc0", lw=0, zorder=2)
    kinds = [("both metrics", sc & sv, "#2f3e46"), ("count-based only", sc & ~sv, "#c0392b"),
             ("coverage-based only", sv & ~sc, "#b8a48a")]
    for lab, mask, col in kinds:
        ax.scatter(d.l2fc_cov[mask], d.l2fc_cnt[mask], s=24, color=col, edgecolor="#333333", lw=0.4,
                   zorder=4, label=f"{lab} ({int(mask.sum())})")
    nm = names.module_name.to_dict()
    hi = d[sc | sv]
    up = [(r.l2fc_cov, r.l2fc_cnt, short(m, nm, False)) for m, r in hi[hi.l2fc_cnt + hi.l2fc_cov > 0].iterrows()]
    dn = [(r.l2fc_cov, r.l2fc_cnt, short(m, nm, False)) for m, r in hi[hi.l2fc_cnt + hi.l2fc_cov <= 0].iterrows()]
    label_column(ax, up, x_col=-lim + 0.15, ha="left", fontsize=5.0, color="#444444")
    label_column(ax, dn, x_col=lim - 0.15, ha="right", fontsize=5.0, color="#444444")
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("Coverage-based log$_2$FC")
    ax.set_ylabel("Count-based log$_2$FC")
    ax.text(0.98, 0.98, "NVNN / (V + NV)", transform=ax.transAxes, ha="right", va="top", fontsize=5.8, color="#444444")
    ax.legend(loc="lower left", fontsize=5.2, frameon=True, framealpha=0.9, edgecolor="#cccccc",
              title="4/4 NVNN vs 0/9 separated", title_fontsize=5.2, handletextpad=0.1, borderpad=0.3)
    return dict(both=int((sc & sv).sum()), count_only=int((sc & ~sv).sum()), coverage_only=int((sv & ~sc).sum()))


def f4_convergence(ax, tax, mcov, mcnt):
    """Bray-Curtis distance between every pair of P10 cultures, for composition (RpS16)
    and function (both metrics, each sample divided by its total): pairs within V, within
    NV, within NVNN, and pairs from different media. R2 and p: PERMANOVA on the same
    distances (fraction of the variation explained by the medium)."""
    from scipy.spatial.distance import pdist, squareform
    sets = [("Composition (RpS16)", tax), ("Function, coverage-based", C.compositional(mcov)),
            ("Function, relative pathway copies", C.compositional(mcnt))]
    kinds = ["V", "NV", "NVNN", "between"]
    rows, stats = [], []
    rng = np.random.default_rng(3)
    for lab, m in sets:
        g = media_groups(list(m.columns))
        cols = sum(g.values(), [])
        D = pd.DataFrame(squareform(pdist(m[cols].T.values, "braycurtis")), index=cols, columns=cols)
        of = {s_: k for k, v in g.items() for s_ in v}
        for i, a in enumerate(cols):
            for b in cols[i + 1:]:
                rows.append(dict(set=lab, kind=of[a] if of[a] == of[b] else "between", d=D.loc[a, b]))
        r = C.bray_pcoa(m, g)
        n, k = len(cols), len(g)
        stats.append((r["F"] * (k - 1) / (r["F"] * (k - 1) + (n - k)), r["p"]))
    d = pd.DataFrame(rows)
    for i, (lab, _) in enumerate(sets):
        x_set = i * 5.2
        for j, kind in enumerate(kinds):
            x = x_set + j
            sub = d[(d.set == lab) & (d.kind == kind)].d
            col = MEDIA.get(kind, "#8c8c8c")
            ax.boxplot(sub, positions=[x], widths=0.6, showfliers=False, patch_artist=True, zorder=2,
                       boxprops=dict(facecolor="#f2f2f2" if kind == "between" else "white",
                                     edgecolor="#555555", lw=0.6),
                       medianprops=dict(color="#222222", lw=1.2), whiskerprops=dict(color="#555555", lw=0.6),
                       capprops=dict(color="#555555", lw=0.6))
            ax.scatter(x + rng.uniform(-0.17, 0.17, len(sub)), sub, s=7, color=col, lw=0, alpha=0.9, zorder=3)
        pv = stats[i][1]
        ptxt = f"$p$ = {pv:.3f}" if pv >= 0.001 else "$p$ < 0.001"
        ax.text(x_set + 1.5, 1.13, lab, ha="center", va="bottom", fontsize=6.8, fontweight="bold")
        ax.text(x_set + 1.5, 1.03, f"medium explains $R^2$ = {stats[i][0]:.2f}, {ptxt}",
                ha="center", va="bottom", fontsize=5.8, color="#444444")
    ticks = [i * 5.2 + j for i in range(len(sets)) for j in range(len(kinds))]
    ax.set_xticks(ticks, ["within\nV", "within\nNV", "within\nNVNN", "between\nmedia"] * len(sets),
                  fontsize=5.8)
    for t, kind in zip(ax.get_xticklabels(), kinds * len(sets)):
        t.set_color(MEDIA.get(kind, "#555555"))
    ax.set_xlim(-0.7, (len(sets) - 1) * 5.2 + 3.7); ax.set_ylim(0, 1.0)
    ax.set_ylabel("Bray-Curtis distance\nbetween two cultures")
    ax.spines[["top", "right"]].set_visible(False)
    return d.groupby(["set", "kind"]).d.median().unstack().round(3), stats


# Fig 4E rows: the functions of Fig 3E that the media share, plus every module
# whose carriers differ at p < 0.05 between NVNN and nitrate-replete cultures
E4_GROUPS = [
    ("Photosynthesis,\ncarbon fixation", [("M00161", "Photosystem II"), ("M00165", "Calvin cycle")]),
    ("Nitrogen", [("M00175", "Nitrogen fixation"), ("M00615", "Nitrate assimilation"),
                  ("M00531", "Assimilatory nitrate reduction")]),
    ("Vitamins and\nco-factors", [("M00122", "Cobalamin"), ("M00119", "Pantothenate"),
                                  ("M00123", "Biotin"), ("M00899", "Thiamine salvage"),
                                  ("M00116", "Menaquinone")]),
    ("Lipids and\namino acids", [("M00091", "Phosphatidylcholine"), ("M00098", "Acylglycerol degradation"),
                                 ("M00020", "Serine biosynthesis"), ("M00035", "Methionine degradation")]),
]
CYANO_C, HETERO_C = "Greens", "Blues"


def f4_carriers(page, x0, y0, x1, y1):
    """Fig 4E: cyanobacterial share of the coverage of genomes encoding each module,
    per P10 replicate (shared genome reference, pipeline/09), as Fig 3E. Right:
    exact Mann-Whitney p, NVNN (n = 4) against V + NV (n = 9)."""
    share, test, comm, pc = med.partition_test()
    _, sh = med.shared_partition()
    cols = C.MEDIA_V + C.MEDIA_NV + C.MEDIA_NVNN
    xs, xpos = [], 0
    for j, smp in enumerate(cols):
        if j and smp.split("_")[0] != cols[j - 1].split("_")[0]:
            xpos += 0.4
        xs.append(xpos); xpos += 1
    n_rows = sum(len(r) for _, r in E4_GROUPS) + len(E4_GROUPS)
    ax = page.ax(x0 + 150, y0, x1 - 52, y1)
    greens, blues = mpl.colormaps[CYANO_C], mpl.colormaps[HETERO_C]
    y = n_rows - 1
    for gi, (g, members) in enumerate(E4_GROUPS):
        if gi:
            ax.plot([-6.9, xpos + 2.0], [y + 0.5, y + 0.5], color="#333333", lw=0.7, clip_on=False)
        ax.text(-6.9, y - 0.05, g.replace("\n", " "), ha="left", va="center", fontsize=6.5, fontweight="bold")
        y -= 1
        ys = []
        for m, lab in members:
            ys.append(y)
            for xj, smp in zip(xs, cols):
                v = share.loc[m, smp] if m in share.index else np.nan
                v = np.nan if np.isnan(v) else float(np.clip(np.round(v), 0, 100))
                a, b = xj, y - 0.45
                if np.isnan(v):
                    ax.add_patch(mpatches.Rectangle((a, b), 1, 0.9, fc="#f3f3f3", ec="#bbbbbb", lw=0.3))
                    continue
                ax.add_patch(mpatches.Polygon([[a, b], [a, b + 0.9], [a + 1, b + 0.9]],
                             fc=greens(0.08 + 0.85 * v / 100), ec="none"))
                ax.add_patch(mpatches.Polygon([[a, b], [a + 1, b], [a + 1, b + 0.9]],
                             fc=blues(0.08 + 0.85 * (100 - v) / 100), ec="none"))
                ax.add_patch(mpatches.Rectangle((a, b), 1, 0.9, fc="none", ec="#555555", lw=0.3))
                ax.text(a + 0.06, b + 0.84, f"{v:.0f}", ha="left", va="top", fontsize=4.3,
                        color="white" if v > 60 else "#222222")
                ax.text(a + 0.94, b + 0.06, f"{100 - v:.0f}", ha="right", va="bottom", fontsize=4.3,
                        color="white" if v < 40 else "#222222")
            ax.text(-0.25, y, f"{lab}  {m}", ha="right", va="center", fontsize=6.2)
            if m in test.index:
                pv = test.loc[m, "p"]
                ax.text(xpos + 0.3, y, f"{pv:.3f}" if pv >= 0.001 else "<0.001", ha="left", va="center",
                        fontsize=6, fontweight="bold" if pv < 0.05 else "normal",
                        color="#222222" if pv < 0.05 else "#777777")
            else:
                ax.text(xpos + 0.3, y, "n.t.", ha="left", va="center", fontsize=6, color="#aaaaaa")
            y -= 1
    ax.text(xpos + 0.3, n_rows - 0.1, "$p$", ha="left", va="bottom", fontsize=6.5, style="italic")
    for xj, smp in zip(xs, cols):
        ax.text(xj + 0.5, -0.65, smp.split("_")[-1], ha="center", va="top", fontsize=5.6)
        ax.text(xj + 0.5, -1.45, f"{sh.loc[smp, 'cyano_of_mapped']:.0f}", ha="center", va="top", fontsize=5.4,
                color="#2e7d32")
    ax.text(-0.25, -1.45, "cyanobacteria (% of genome coverage)", ha="right", va="top", fontsize=5.6,
            color="#2e7d32")
    for g in ("V", "NV", "NVNN"):
        idx = [xj for xj, smp in zip(xs, cols) if smp.startswith(g + "_")]
        ax.text(np.mean(idx) + 0.5, n_rows - 0.1, g, ha="center", va="bottom", fontsize=7, fontweight="bold",
                color=MEDIA[g])
    ax.set_xlim(-0.1, xpos + 0.1); ax.set_ylim(-2.1, n_rows + 0.5)
    ax.axis("off")
    norm = Normalize(0, 100)
    page.text(x0 + 16, y1 + 10, "Cell: % of the coverage of genomes encoding the module that is cyanobacterial "
              "(upper left) or heterotrophic (lower right).\nGrey: no detected genome encodes it.  "
              "$p$: exact Mann-Whitney, NVNN (n = 4) vs V + NV (n = 9); none significant after BH correction.",
              fontsize=5.6, color="#444444", ha="left", linespacing=1.4)
    for k, (cmap, lab) in enumerate([(greens, "Cyanobacteria (%)"), (blues, "Heterotrophs (%)")]):
        cax = page.ax(x0 + 150 + 120 * k, y1 + 34, x0 + 240 + 120 * k, y1 + 39)
        cb = page.fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=mpl.colors.ListedColormap(
            cmap(np.linspace(0.08, 0.93, 256)))), cax=cax, orientation="horizontal")
        cb.set_label(lab, fontsize=5.8, labelpad=1)
        cb.ax.tick_params(labelsize=5, length=1.5)
        cb.outline.set_linewidth(0.4)
    share.to_csv(C.OUT / "fig4E_carriers.csv")
    test.to_csv(C.OUT / "fig4E_carriers_test.csv")
    return dict(test=test, community=comm, community_p=pc)


PHYLUM_COL = {"Cyanobacteriota": "#2e9e4f", "Pseudomonadota": "#d1495b", "Bacteroidota": "#3a7ca5",
              "Actinomycetota": "#8d6cab", "Acidobacteriota": "#c08a2e", "Planctomycetota": "#7f8c8d",
              "Verrucomicrobiota": "#5a9e9a", "Unassigned": "#bdbdbd"}
PLACE_MAX = 0.6     # a placed cluster borrows a family name only this close to a named tip


def lineage_prefix(phylum, klass, family):
    """Fig 4B label prefix: the GTDB class, or the family for cyanobacteria, which share one class."""
    out = klass.where(phylum != C.CYANO_PHYLUM, family)
    return out.where(phylum != "Unassigned")


def media_genera(min_share=0.01, max_rows=20):
    """Mean relative RpS16 coverage per genus, P10, with tree placement for unnamed clusters."""
    run = "rps16_media_endpoints"
    tax = C.ecophylo_taxonomy(run)
    e = pd.read_csv(C.DATA / "ecophylo_coverage.csv")
    e = e[e.run == run]
    m = e.pivot_table(index="cluster", columns="sample", values="mean_coverage_Q2Q3").fillna(0)
    rel = m / m.sum()

    def label(c):
        r = tax.loc[c]
        if not r.placed:
            return r.t_genus if isinstance(r.t_genus, str) else f"{r.t_family} (unclass.)"
        if r.patristic <= PLACE_MAX:
            return f"unnamed {r.t_family}"
        if r.t_phylum == C.CYANO_PHYLUM:
            # Too far from any named tip to trust, and inside the cyanobacteria: in the cultures
            # where it dominates, reads map to no cyanobacterial genome (pipeline/09), so it is
            # most likely the plastid copy of a photosynthetic eukaryote. Shown as unassigned.
            return "unassigned (near Cyanobiaceae)"
        return f"unresolved (nearest {r.t_family})"
    lab = pd.Series({c: label(c) for c in rel.index})
    g = rel.groupby(lab).sum()
    fam = tax.groupby(lab.reindex(tax.index)).t_family.first()
    phy = tax.groupby(lab.reindex(tax.index)).t_phylum.first()
    phy[[l for l in phy.index if l == "unassigned (near Cyanobiaceae)"]] = "Unassigned"
    klass = tax.groupby(lab.reindex(tax.index)).t_class.first()
    groups = media_groups(g.columns)
    mean = pd.DataFrame({k: g[v].mean(axis=1) for k, v in groups.items()})
    prev = pd.DataFrame({k: (g[v] > 0.001).sum(axis=1) for k, v in groups.items()})
    nrep = prev.V + prev.NV
    cls = np.where((prev.NVNN >= 3) & (nrep <= 1), "only without nitrate",
                   np.where((prev.NVNN <= 1) & (nrep >= 3), "only with nitrate", "shared"))
    t = mean.assign(cls=cls, family=fam, phylum=phy, peak=mean.max(axis=1),
                    rank_above=lineage_prefix(phy, klass, fam))
    t = t[(t.peak >= min_share) | ((t.cls == "only without nitrate") & (t.NVNN >= 0.004))]
    order = {"only without nitrate": 0, "shared": 1, "only with nitrate": 2}
    t = t.assign(o=t.cls.map(order)).sort_values(["o", "peak"], ascending=[True, False])
    t = t.groupby("cls", group_keys=False).head(8)
    return g.loc[t.index], t, groups


def rpl19_nvnn_only():
    """Genera found only without nitrate with the second marker (RpL19, P10)."""
    meta = pd.read_csv(C.DATA / "media_samples.csv").set_index("sample")
    p10 = meta[meta.passage == "P10"].sra_name
    run = "rpl19_media_allpassages"
    tax = C.ecophylo_taxonomy(run)
    e = pd.read_csv(C.DATA / "ecophylo_coverage.csv")
    e = e[(e.run == run) & e["sample"].isin(p10.index)]
    m = e.pivot_table(index="cluster", columns="sample", values="mean_coverage_Q2Q3").fillna(0).rename(columns=p10)
    rel = m / m.sum()
    lab = pd.Series(np.where(tax.placed, "unnamed " + tax.t_family.astype(str), tax.t_genus.astype(str)), index=tax.index)
    g = rel.groupby(lab.reindex(rel.index)).sum()
    nv = [c for c in g.columns if c.startswith("NVNN_")]
    rep = [c for c in g.columns if not c.startswith("NVNN_")]
    return set(g.index[((g[nv] > 0.001).sum(axis=1) >= 3) & ((g[rep] > 0.001).sum(axis=1) <= 1)])


LAB_SLASH, LAB_GENUS = 62, 92      # Fig 4B: width of the class column and of the genus column, pt


def f4_genera(page, x0, y0, x1, y1):
    g, t, groups = media_genera()
    confirmed = rpl19_nvnn_only()
    cols = sum(groups.values(), [])
    ax = page.ax(x0 + LAB_SLASH + LAB_GENUS, y0, x1 - 8, y1)
    xs, xpos = [], 0
    for j, smp in enumerate(cols):
        if j and smp.split("_")[0] != cols[j - 1].split("_")[0]:
            xpos += 0.8
        xs.append(xpos); xpos += 1
    n = len(t)
    prev_cls = None
    y = n - 1
    ypos = {}
    for gen, r in t.iterrows():
        if prev_cls is not None and r.cls != prev_cls:
            ax.axhline(y + 0.5, color="#bbbbbb", lw=0.5, ls=":")
        prev_cls = r.cls
        ypos[gen] = y
        col = PHYLUM_COL.get(r.phylum, "#999999")
        v = g.loc[gen, cols].values
        ax.scatter(xs, [y] * len(xs), s=np.sqrt(v) * 150, color=col, alpha=0.85, edgecolor="#333333", lw=0.3, zorder=3)
        # one right-aligned label: "Class/" (family for cyanobacteria) in grey, then the genus
        show = gen[len("unnamed "):] + "*" if gen.startswith("unnamed ") else gen
        named = show[:1].isupper() and " " not in show and not any(ch.isdigit() for ch in show)
        ty = y0 + (n + 0.6 - y) / (n + 1.3) * (y1 - y0)          # row y in page points (ylim below)
        uncertain = gen.startswith(("unresolved", "unassigned"))   # no class: placement too uncertain
        if uncertain:
            show = show.replace("nearest", "near")
        x_end = x0 + LAB_SLASH + LAB_GENUS - 4
        tg = page.fig.text(x_end / page.w, 1 - ty / page.h, show, ha="right", va="center", fontsize=6,
                           style="italic" if named and not show.endswith("*") else "normal")
        if isinstance(r.rank_above, str) and not uncertain:
            w = tg.get_window_extent(page.fig.canvas.get_renderer()).width * 72 / page.fig.dpi
            page.fig.text((x_end - w) / page.w, 1 - ty / page.h, f"{r.rank_above}/", ha="right",
                          va="center", fontsize=6, color="#555555")
        y -= 1
    for cls_, ys in t.assign(y=pd.Series(ypos)).groupby("cls").y:
        ax.text(xpos + 0.3, ys.mean(), cls_.replace("only ", "predominantly\n"), ha="left", va="center", fontsize=6,
                color="#444444", style="italic")
    for gname, members in groups.items():
        idx = [xs[cols.index(s_)] for s_ in members]
        ax.text(np.mean(idx), n - 0.2, gname, ha="center", va="bottom", fontsize=7.5, fontweight="bold",
                color=MEDIA[gname])
        ax.axvspan(min(idx) - 0.5, max(idx) + 0.5, color=MEDIA[gname], alpha=0.07, lw=0, zorder=0)
    ax.set_xlim(-0.6, xpos + 6.8); ax.set_ylim(-0.7, n + 0.6)
    ax.axis("off")
    hs = [plt.Line2D([], [], marker="o", ls="", color=c, markersize=4.5, label=p_) for p_, c in PHYLUM_COL.items()
          if p_ in set(t.phylum)]
    hs += [plt.scatter([], [], s=np.sqrt(v_) * 150, color="#cccccc", edgecolor="#333333", lw=0.3, label=f"{v_:.0%}")
           for v_ in (0.01, 0.1, 0.4)]
    ax.legend(handles=hs, loc="upper left", bbox_to_anchor=(-0.95, -0.02), ncol=6, fontsize=5.4,
              frameon=False, handletextpad=0.1, columnspacing=0.7)
    t.assign(rpl19_confirmed=t.index.isin(confirmed)).to_csv(C.OUT / "fig4B_genera.csv")
    return t


F4E_ROWS = [("M00175", "Nitrogen fixation"), ("M00950", "Biotin"), ("M00119", "Pantothenate"),
            ("M00912", "NAD biosynthesis"), ("M00930", "Menaquinone"), ("M00124", "Pyridoxal phosphate"),
            ("M00161", "Photosystem II"), ("M00165", "Calvin cycle"), ("M00615", "Nitrate assimilation"),
            ("M00034", "Methionine salvage")]


def f4_who_does_what(page, x0, y0, x1, y1, mcov, mcnt):
    """Fig 4E: per P10 replicate, coverage-based and count-based abundance of key
    nitrate-responsive modules, each scaled to its highest replicate."""
    groups = media_groups(mcov.columns)
    cols = sum(groups.values(), [])
    ax = page.ax(x0 + 120, y0, x1, y1)
    n = len(F4E_ROWS)
    cmap = mpl.colormaps["Greys"]
    xpos, starts = 0, []
    for bi, (ttl, mat) in enumerate([("Coverage-based", mcov), ("Count-based", mcnt)]):
        if bi:
            xpos += 1.6
        starts.append(xpos)
        for j, smp in enumerate(cols):
            if j and smp.split("_")[0] != cols[j - 1].split("_")[0]:
                xpos += 0.35
            for i, (m, _) in enumerate(F4E_ROWS):
                yy = n - 1 - i
                top = mat.loc[m, cols].max() if m in mat.index else 0
                v = mat.loc[m, smp] / top if top > 0 else 0
                ax.add_patch(mpatches.Rectangle((xpos, yy - 0.45), 0.92, 0.9, fc=cmap(0.04 + 0.9 * v), ec="white", lw=0.4))
            ax.text(xpos + 0.46, -0.6, smp.split("_")[-1], ha="center", va="top", fontsize=5)
            xpos += 1
        ax.text((starts[-1] + xpos) / 2, n - 0.1, ttl, ha="center", va="bottom", fontsize=7, fontweight="bold")
    for st in starts:
        x = st
        for gname, members in groups.items():
            ax.text(x + len(members) / 2, -1.25, gname, ha="center", va="top", fontsize=6.3, color=MEDIA[gname],
                    fontweight="bold")
            x += len(members) + 0.35
    for i, (m, lab) in enumerate(F4E_ROWS):
        ax.text(-0.3, n - 1 - i, f"{lab}  {m}", ha="right", va="center", fontsize=6.4)
    ax.set_xlim(-0.2, xpos + 0.2); ax.set_ylim(-2.0, n + 0.7)
    ax.axis("off")
    cax = page.ax(x1 - 150, y1 + 10, x1 - 30, y1 + 15)
    cb = page.fig.colorbar(mpl.cm.ScalarMappable(cmap=mpl.colors.ListedColormap(cmap(np.linspace(0.04, 0.94, 256)))),
                           cax=cax, orientation="horizontal")
    cb.set_ticks([0, 0.5, 1]); cb.ax.tick_params(labelsize=5.5, length=1.5)
    cb.set_label("Module abundance, fraction of the highest replicate", fontsize=5.8, labelpad=2)
    cb.outline.set_linewidth(0.4)
    return {m: (mcov.loc[m, cols].round(3).to_dict() if m in mcov.index else None) for m in ("M00161", "M00165")}


def compose_si_media_replicates():
    """SI figure: per-replicate coverage- and count-based values of key nitrate-responsive modules (former Fig 4E)."""
    small_style()
    res, mcov, mcnt, names = med.contrasts("rps16")
    W, H = 510, 190
    pg = Page(W, H)
    pg.title(W / 2, 0, "Key nitrate-responsive modules per replicate")
    f4_who_does_what(pg, 4, 22, 506, 150, mcov, mcnt)
    (C.OUT / "si").mkdir(exist_ok=True)
    pg.save(C.OUT / "si" / "figS_media_replicates.pdf", {})


def compose_si_media_partition():
    """SI figure: who carries each function in the media, shared genome reference (former Fig 4E draft)."""
    small_style()
    W, H = 510, 330
    pg = Page(W, H)
    pg.title(W / 2, 0, "Cyanobacterial share of each function, P10 replicates")
    out = f4_carriers(pg, 4, 30, 506, 262)
    (C.OUT / "si").mkdir(exist_ok=True)
    pg.save(C.OUT / "si" / "figS_media_partition.pdf", {})
    return out


F2_ORDER = ["M_132", "M_133", "M_134", "M_135", "M_63", "M_67", "M_71", "M_77"] + C.SOURCE
F2_MARKER = {"Enrichment": "o", "Source": "s"}


def f2_bars(ax, ph):
    """Stacked phylum bars, enrichments then reservoir, colours as Fig 2A."""
    named = [p for p in C.PHYLUM_COLOR if p not in ("Other", "Unknown")]
    m = ph.reindex(named).fillna(0.0)
    m.loc["Other"] = ph.loc[[p for p in ph.index if p not in named]].sum()
    m = m[F2_ORDER]
    x = np.arange(len(F2_ORDER)) + np.where(np.isin(F2_ORDER, C.SOURCE), 0.6, 0)
    bottom = np.zeros(len(F2_ORDER))
    for p in m.index:
        ax.bar(x, m.loc[p], bottom=bottom, width=0.86, color=C.PHYLUM_COLOR[p], lw=0)
        bottom += m.loc[p].values
    ax.set_xticks(x, F2_ORDER, rotation=90, fontsize=5.4)
    ax.set_xlim(x[0] - 0.6, x[-1] + 0.6)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Relative coverage", fontsize=6.5, labelpad=2)
    ax.tick_params(axis="x", length=0, pad=2)
    ax.tick_params(axis="y", labelsize=6)
    for g, members in (("Enrichment", C.ENRICHMENT), ("Source", C.SOURCE)):
        xs = x[np.isin(F2_ORDER, members)]
        ax.add_patch(mpatches.FancyBboxPatch((xs[0] - 0.42, 1.025), xs[-1] - xs[0] + 0.84, 0.06,
                     boxstyle="round,pad=0,rounding_size=0.02", transform=ax.transData,
                     clip_on=False, facecolor=C.GROUP_COLOR[g], alpha=0.35, lw=0))
        ax.text(xs.mean(), 1.055, g.upper(), ha="center", va="center", fontsize=5.8)
    ax.spines[["top", "right"]].set_visible(False)


def f2_pcoa(ax, mat, arrows, label, color=None, min_share=0.01):
    """Bray-Curtis PCoA (common.bray_pcoa) with 95% ellipses. Arrows: correlation of a
    feature's relative abundance with PCo1 and PCo2, for the `arrows` features with the
    strongest correlation among those at >= `min_share` mean share in either group; one
    arrow per label, and none closer than 25 degrees to an arrow already drawn."""
    groups = {"Source": C.SOURCE, "Enrichment": C.ENRICHMENT}
    r = C.bray_pcoa(mat, groups)
    sc = r["scores"]
    for g in ("Source", "Enrichment"):
        xy = sc.loc[[s for s in groups[g] if s in sc.index]].values
        ellipse(ax, xy, C.GROUP_COLOR[g], n_std=CI95)
        ax.scatter(xy[:, 0], xy[:, 1], s=22, marker=F2_MARKER[g], color=C.GROUP_COLOR[g],
                   edgecolor="#333333" if g == "Enrichment" else "white", lw=0.4, zorder=3)
    m = mat[sc.index]
    big = m.index[(m[C.SOURCE].mean(1) >= min_share) | (m[C.ENRICHMENT].mean(1) >= min_share)]
    cor = pd.DataFrame({ax_: m.loc[big].T.corrwith(sc[ax_]) for ax_ in ("PCo1", "PCo2")})
    ext = 0.85 * min(np.abs(sc["PCo1"]).max(), np.abs(sc["PCo2"]).max() * 1.6)
    drawn, seen = [], set()
    for f in (cor ** 2).sum(1).sort_values(ascending=False).index:
        ang = np.degrees(np.arctan2(cor.loc[f, "PCo2"], cor.loc[f, "PCo1"]))
        if label(f) in seen or any(abs((ang - a + 180) % 360 - 180) < 25 for a in drawn):
            continue
        drawn.append(ang); seen.add(label(f))
        dx, dy = cor.loc[f] * ext
        col = color(f) if color else "#555555"
        ax.annotate("", xy=(dx, dy), xytext=(0, 0), zorder=2,
                    arrowprops=dict(arrowstyle="-|>", lw=0.8, color=col, mutation_scale=6))
        ax.annotate(label(f), xy=(dx, dy), xytext=(3 * np.cos(np.radians(ang)) * 1.4,
                                                     3 * np.sin(np.radians(ang)) * 1.4),
                    textcoords="offset points", fontsize=5.8, color=col, zorder=4,
                    ha="left" if abs(ang) < 70 else "right" if abs(ang) > 110 else "center",
                    va="bottom" if ang > 20 else "top" if ang < -20 else "center",
                    fontstyle="normal" if color else "italic", annotation_clip=False,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.75))
        if len(drawn) == arrows:
            break
    ax.set_xlabel(f"PCo1 ({r['explained'][0]:.1f}%)")
    ax.set_ylabel(f"PCo2 ({r['explained'][1]:.1f}%)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.margins(0.18)
    x0, x1 = ax.get_xlim()
    ax.set_xlim(x0 - 0.22 * (x1 - x0), x1)     # room for labels of left-pointing arrows
    ax.text(0.03, 0.03, f"PERMANOVA {fmt_p(r['p'])}\nPERMDISP {fmt_p(r['disp_p'])}",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=5.8, color="#555555",
            linespacing=1.25)
    return dict(PCo1=r["explained"][0], PCo2=r["explained"][1], F=r["F"], p=r["p"],
                disp_p=r["disp_p"], **{f"centroid_{k}": v for k, v in r["centroid_dist"].items()})


def f2_rank_pcoa(ax, mat, title):
    """Bray-Curtis PCoA at one rank with 95% ellipses; PERMANOVA R2 is the share of the
    variation between samples explained by reservoir vs enrichment (the separation)."""
    groups = {"Source": C.SOURCE, "Enrichment": C.ENRICHMENT}
    r = C.bray_pcoa(mat, groups)
    for g, members in groups.items():
        xy = r["scores"].loc[members].values
        ellipse(ax, xy, C.GROUP_COLOR[g], n_std=CI95)
        ax.scatter(xy[:, 0], xy[:, 1], s=14, marker=F2_MARKER[g], color=C.GROUP_COLOR[g],
                   edgecolor="#333333" if g == "Enrichment" else "white", lw=0.4, zorder=3)
    ax.set_title(title, fontsize=7, fontweight="bold", pad=9)
    ax.text(0.5, 1.0, f"$R^2$ = {r['R2']:.2f}, {fmt_p(r['p'])}", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=6, color="#444444")
    ax.set_xlabel(f"PCo1 ({r['explained'][0]:.0f}%)", fontsize=6.5, labelpad=1)
    ax.set_ylabel(f"PCo2 ({r['explained'][1]:.0f}%)", fontsize=6.5, labelpad=1)
    ax.margins(0.15)
    score_ticks(ax)
    ax.spines[["top", "right"]].set_visible(False)
    return r


def f2_pair_strip(ax, mat, ylabels=True):
    """Convergence on the distances behind the PCoA above: each dot is one pair of samples
    from the same group (28 enrichment pairs, 28 reservoir pairs), placed at the share of
    relative coverage the two have in common (1 - Bray-Curtis); tick: median; shaded bar:
    interquartile range."""
    from scipy.spatial.distance import pdist
    rng = np.random.default_rng(2)
    med = {}
    for y, g, members in ((1, "Enrichment", C.ENRICHMENT), (0, "Source", C.SOURCE)):
        v = (1 - pdist(mat[members].T.values, "braycurtis"))
        ax.scatter(v, y + rng.uniform(-0.17, 0.17, len(v)), s=5, marker=F2_MARKER[g], lw=0.25,
                   color=C.GROUP_COLOR[g], edgecolor="#333333" if g == "Enrichment" else "white",
                   alpha=0.85, zorder=3)
        med[g] = np.median(v)
        q1, q3 = np.percentile(v, [25, 75])
        ax.plot([q1, q3], [y, y], color=C.GROUP_COLOR[g], lw=4.5, alpha=0.3,
                solid_capstyle="butt", zorder=2)
        ax.plot([med[g]] * 2, [y - 0.3, y + 0.3], color="#222222", lw=1.1, zorder=4)
        up = g == "Enrichment"
        ax.text(med[g], y + (0.34 if up else -0.34), f"{med[g]:.2f}", ha="center",
                va="bottom" if up else "top", fontsize=6, color="#222222")
    ax.set_xlim(0, 1); ax.set_ylim(-0.85, 1.85)
    similarity_ticks(ax)
    ax.set_yticks([1, 0], ["Enrichment\npairs", "Source\npairs"] if ylabels else ["", ""], fontsize=6)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("")              # one shared label under both strips, set in compose_fig2
    ax.grid(axis="x", color=GRID, lw=0.5, zorder=0)
    ax.spines[["top", "right", "left"]].set_visible(False)
    return med


def f2_dominant(ax):
    """Heatmap of the populations that dominate an enrichment (> 20% in at least one):
    relative RpL19 coverage in each enrichment and each reservoir sample on a log colour
    scale (greys); white, not detected. Rows ordered by the
    first enrichment each population dominates. Returns the rows and the colour mappable."""
    import taxonomy as tx
    from matplotlib.colors import LogNorm
    rel = tx.load_clusters(relative=True)
    _, dom = tx.population_sharing()
    cols = ["M_63", "M_67", "M_132", "M_71", "M_133", "M_134", "M_135", "M_77"]
    dom["first"] = [min(cols.index(c) for c in d.split(",")) for d in dom.dominant_in]
    dom = dom.sort_values(["first", "max_enrichment"], ascending=[True, False]).reset_index(drop=True)
    samples = cols + C.SOURCE
    xs = np.arange(len(samples)) + np.where(np.isin(samples, C.SOURCE), 0.7, 0)
    from matplotlib.colors import ListedColormap      # grey, so as not to compete with the group and
    cmap = ListedColormap(mpl.colormaps["Greys"](np.linspace(0.22, 1, 256)))   # phylum hues; no near-white
    norm = LogNorm(vmin=0.01, vmax=100)
    for i, r in dom.iterrows():
        y = len(dom) - 1 - i
        for x, c in zip(xs, samples):
            v = 100 * rel.loc[r.population, c]
            ax.add_patch(mpatches.Rectangle((x - 0.5, y - 0.5), 1, 1, lw=0.4, ec="#e3e3e3",
                                            fc=cmap(norm(max(v, 0.01))) if v > 0 else "white"))
    ax.set_xlim(xs[0] - 0.5, xs[-1] + 0.5); ax.set_ylim(-0.5, len(dom) - 0.5)
    ax.set_xticks(xs, samples, rotation=90, fontsize=5)
    name = lambda r: r.genus if isinstance(r.genus, str) else "Leptolyngbyaceae*"
    ax.set_yticks(range(len(dom))[::-1], [name(r) for r in dom.itertuples()], fontsize=6)
    for t in ax.get_yticklabels():     # italics for named genera only, not GTDB placeholders
        t.set_fontstyle("italic" if re.fullmatch(r"[A-Z][a-z]+", t.get_text()) else "normal")
    ax.tick_params(length=0, pad=1.5)
    for sp in ax.spines.values():
        sp.set_visible(False)
    from matplotlib.transforms import blended_transform_factory
    tr = blended_transform_factory(ax.transData, ax.transAxes)     # strips placed exactly as in panel B
    for g, members in (("Enrichment", cols), ("Source", C.SOURCE)):
        xg = xs[np.isin(samples, members)]
        ax.add_patch(mpatches.FancyBboxPatch((xg[0] - 0.42, 1.025), xg[-1] - xg[0] + 0.84, 0.06,
                     boxstyle="round,pad=0,rounding_size=0.02", transform=tr, clip_on=False,
                     facecolor=C.GROUP_COLOR[g], alpha=0.35, lw=0))
        ax.text(xg.mean(), 1.055, g.upper(), ha="center", va="center", fontsize=5.8, transform=tr)
    return dom, mpl.cm.ScalarMappable(norm=norm, cmap=cmap)


F2_ECOPHYLO_LEGEND = [      # colours as drawn in data/schematics/fig2A_ecophylo.pdf
    ("Sample group", [("Enrichment (8)", "#75b79e"), ("Source (8)", "#31326f")]),
    ("RpL19 sequence from", [("contigs (508)", "#808080"), ("MAGs (151)", "#ff99ff")]),
    ("Phylum", [(p, c) for p, c in C.PHYLUM_COLOR.items() if p != "Other"]),
]


def f2_ecophylo_legend(pg, x, y, dy=11.4, gap=9):
    """Legends of the EcoPhylo panel at the figure's legend size, spread over the panel height."""
    for k, (head, items) in enumerate(F2_ECOPHYLO_LEGEND):
        pg.text(x, y, head, fontsize=6.6, fontweight="bold", ha="left")
        y += dy + 1.5
        for lab, col in items:
            pg.fig.patches.append(mpatches.Rectangle((x / pg.w, 1 - (y + 5.6) / pg.h), 8 / pg.w, 5.6 / pg.h,
                                                     transform=pg.fig.transFigure, color=col, lw=0))
            pg.text(x + 11, y - 0.8, lab, fontsize=6.2, ha="left")
            y += dy
        y += gap


def compose_fig2():
    """Fig 2: A EcoPhylo (anvi'o, vector crop); B phylum composition; C Bray-Curtis PCoA at
    class and at order (separation from the reservoir, PERMANOVA R2) with, under each, the
    coverage shared by every pair of enrichments and of reservoir samples on the same
    distances (convergence; at class the low enrichment pairs are the seven with M_135);
    D the dominant populations, rare or undetected in the reservoir. Panels run from coarse
    to fine. Second-row axes share one top (Y0) and bottom (Y1) edge."""
    import taxonomy as tx
    small_style()
    ph = tx.phylum_table()
    lt = tx.level_tables()
    W, H = 510, 448
    Y0, Y1 = 304, 414
    pg = Page(W, H)
    pg.letter(4, 272, "B"); pg.title(84, 272, "Phylum composition")
    f2_bars(pg.ax(30, Y0, 140, Y1), ph)
    pg.letter(146, 272, "C"); pg.title(266, 272, "Separation and convergence")
    rows = {}
    for lev, t, x0 in (("class", "Class", 192), ("order", "Order", 278)):
        r = f2_rank_pcoa(pg.ax(x0, Y0, x0 + 56, Y0 + 46), lt[lev], t)
        med = f2_pair_strip(pg.ax(x0, Y1 - 36, x0 + 56, Y1), lt[lev], ylabels=(lev == "class"))
        rows[lev] = dict(R2=r["R2"], p=r["p"], disp_p=r["disp_p"], shared_E=med["Enrichment"],
                         shared_S=med["Source"])
    pg.text(268, Y1 + 13, SIMILARITY_LABEL, fontsize=6.5, ha="center")
    pg.letter(346, 272, "D"); pg.title(430, 272, "Dominant populations")
    winners, sm = f2_dominant(pg.ax(392, Y0, 486, Y1))
    cax = pg.ax(491, Y0 + 6, 494, Y1 - 6)
    cb = pg.fig.colorbar(sm, cax=cax, ticks=[0.01, 1, 100])
    cb.ax.set_yticklabels(["0.01", "1", "100"])
    cb.ax.tick_params(labelsize=5.6, length=1.5, width=0.4, pad=1)
    cb.ax.minorticks_off()
    cb.outline.set_linewidth(0)
    pg.text(440, 285, "% of RpL19 coverage, log scale", fontsize=6, ha="center", color="#555555")
    # panel A: the anvi'o EcoPhylo heatmap, cropped from the published vector figure without its
    # legends, which are redrawn by f2_ecophylo_legend at the type sizes of the other panels
    pg.letter(4, 4, "A")
    f2_ecophylo_legend(pg, 10, 26)
    for t in pg.fig.findobj(mpl.text.Text):   # one step up the type scale, panel letters excepted
        if t.get_text() in "ABCD" or t.get_fontsize() >= 9:
            continue
        if t.get_fontsize() == 8 and t.get_fontweight() == "bold":
            continue                          # panel titles stay at 8 pt, as the EcoPhylo title
        t.set_fontsize(t.get_fontsize() * 1.08)
    pg.save(C.OUT / "FIG2.pdf", {"fig2A_ecophylo": ((112, 0, 496, 266), (118, 0, 586, 324))})
    tops = ph.loc[["Cyanobacteriota", "Pseudomonadota", "Actinomycetota", "Bacteroidota"]]
    out = pd.DataFrame({g: tops[m].apply(lambda r: f"{100 * r.min():.0f}-{100 * r.max():.0f}%", axis=1)
                        for g, m in (("Source", C.SOURCE), ("Enrichment", C.ENRICHMENT))})
    return out, pd.DataFrame(rows).T.round(3), winners


def f4_pcoa(ax, mat):
    """Bray-Curtis PCoA of the media end-points (as Fig 2C-D, 3B), 95% ellipses; PERMANOVA
    across the three media and PERMDISP on the same distances."""
    groups = media_groups(list(mat.columns))
    r = C.bray_pcoa(mat.fillna(0), groups)
    sc = r["scores"]
    for g, members in groups.items():
        xy = sc.loc[members].values
        ellipse(ax, xy, MEDIA[g], n_std=CI95)
        ax.scatter(xy[:, 0], xy[:, 1], s=16, color=MEDIA[g], edgecolor="#333333", lw=0.4,
                   zorder=3, label=g)
    ax.set_xlabel(f"PCo1 ({r['explained'][0]:.1f}%)")
    ax.set_ylabel(f"PCo2 ({r['explained'][1]:.1f}%)")
    ax.margins(0.15)
    pw = ", ".join(f"{a}\u2013{b} {v:.3f}" for (a, b), v in r["pairwise"].items())
    ax.text(0.0, 1.02, f"PERMANOVA {fmt_p(r['p'])}\npairwise $p$ (BH): {pw}", transform=ax.transAxes,
            linespacing=1.3,
            ha="left", va="bottom", fontsize=5.6, color="#555555")
    ax.spines[["top", "right"]].set_visible(False)
    return r


def compose_fig4():
    small_style()
    res, mcov, mcnt, names = med.contrasts("rps16")
    tax = rps16_matrix()
    W, H = 510, 440
    pg = Page(W, H)
    pg.letter(4, 0, "A")
    pg.title(104, 0, "Parallel enrichments in three media")
    pg.letter(196, 0, "B")
    pg.title(350, 0, "Nitrate removal selects a distinct set of genera")
    genera = f4_genera(pg, 196, 26, 510, 196)
    pg.letter(4, 234, "C"); pg.title(140, 234, "Function (relative pathway copies)")
    rc = f4_pcoa(pg.ax(44, 272, 236, 412), C.compositional(mcnt))
    rt = C.bray_pcoa(tax, media_groups(list(tax.columns)))                       # SI only
    rf = C.bray_pcoa(C.compositional(mcov), media_groups(list(mcov.columns)))   # SI only
    hs = [plt.Line2D([], [], ls="", marker="o", ms=4.5, color=MEDIA[g], mec="#333333", mew=0.4,
                     label=g) for g in ("V", "NV", "NVNN")]
    pg.fig.legend(handles=hs, loc="upper right", bbox_to_anchor=(236 / W, 1 - 274 / H), ncol=1,
                  fontsize=6.5, handletextpad=0.1, frameon=False)
    pg.letter(272, 234, "D")
    pg.title(392, 234, "Modules separated by nitrate removal")
    two = f4_two_metrics(pg.ax(318, 254, 500, 412), res, names[["module_name"]])
    # pairwise dissimilarities (former panel E) are Supplementary Figure S19
    # schematic without its own title line (drawn above at the common size)
    pg.save(C.OUT / "FIG4.pdf", {"fig4A_media_design": ((0, 16, 190, 200), (0, 14, 211, 219))})
    pcoa = pd.DataFrame({k: dict(PCo1=r["explained"][0], PCo2=r["explained"][1], F=r["F"], p=r["p"],
                                 disp_p=r["disp_p"], **r["centroid_dist"])
                         for k, r in (("composition", rt), ("function coverage", rf),
                                      ("function count", rc))}).round(3)
    return dict(two=two, genera=genera, pcoa=pcoa)


def compose_fig1():
    """Fig 1: sampling schematic (vector crop); B-C richness and Shannon diversity as
    source | lineages through passages | end-points; D rank abundance by stage. The
    cyanobacterial share is in fig1_samples.csv (composition is Fig 2)."""
    import taxonomy as tx
    small_style()
    cov, m, t = tx.run_lineage_trajectories()
    W, H = 510, 378
    pg = Page(W, H)
    pg.letter(272, 0, "B"); pg.title(406, 0, "Population richness")
    tx.panel_course(pg.ax(312, 18, 500, 166), m, "richness", "RpL19 populations")
    pg.letter(4, 204, "C"); pg.title(142, 204, "Shannon diversity")
    tx.panel_course(pg.ax(48, 222, 236, 350), m, "shannon", "Shannon diversity (H)")
    pg.letter(272, 204, "D"); pg.title(406, 204, "Rank abundance by stage")
    tx.panel_stage_rank(pg.ax(312, 222, 500, 350), cov, m)
    pg.save(C.OUT / "FIG1.pdf", {"fig1A_sampling": (0, 0, 262, 196)})
    return t


def main() -> None:
    run_panels()
    print('FIG1', compose_fig1().round(2).to_string())
    print('FIG2', *compose_fig2(), sep='\n')
    print('FIG3', compose_fig3())
    print('FIG4', compose_fig4())


if __name__ == "__main__":
    main()
