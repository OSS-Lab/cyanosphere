"""
Render the images shown in WORKFLOW.md from the figures that analysis/ writes to outputs/.

    cd analysis && python figures.py && python supplementary.py   # writes outputs/
    python docs/make_images.py                                   # writes docs/img/

Each image is a crop of a composite figure, given in points from its top-left corner.
"""
from pathlib import Path

import pymupdf

ROOT = Path(__file__).resolve().parents[1]
OUT, IMG = ROOT / "outputs", ROOT / "docs" / "img"
DPI = 200

CROPS = {
    "fig1_richness.png":       ("FIG1.pdf", (264, 0, 510, 200)),
    "fig1_diversity_rank.png": ("FIG1.pdf", (0, 200, 510, 378)),
    "fig2_ecophylo.png":       ("FIG2.pdf", (0, 0, 510, 262)),
    "fig2_convergence.png":    ("FIG2.pdf", (0, 263, 510, 444)),
    "fig3_function_pcoa.png":  ("FIG3.pdf", (0, 66, 266, 274)),
    "fig3_two_metrics.png":    ("FIG3.pdf", (264, 66, 510, 268)),
    "fig3_partition.png":      ("FIG3.pdf", (0, 273, 510, 470)),
    "fig4_genera.png":         ("FIG4.pdf", (188, 0, 506, 216)),
    "fig4_media.png":          ("FIG4.pdf", (0, 222, 506, 440)),
}
WHOLE = {
    "si_pcoa_by_rank.png":     "si/S_pcoa_by_rank.pdf",
    "si_count_vs_mags.png":    "fig_si_count_vs_mags.pdf",
}


def main() -> None:
    IMG.mkdir(parents=True, exist_ok=True)
    for name, (pdf, box) in CROPS.items():
        page = pymupdf.open(OUT / pdf)[0]
        page.get_pixmap(dpi=DPI, clip=pymupdf.Rect(*box)).save(IMG / name)
    for name, pdf in WHOLE.items():
        pymupdf.open(OUT / pdf)[0].get_pixmap(dpi=DPI).save(IMG / name)
    print(f"{len(CROPS) + len(WHOLE)} images in {IMG.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
