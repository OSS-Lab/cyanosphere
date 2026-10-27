# Phylogenomic tree of the 88 dereplicated MAGs (SI Fig S8)

Bacteria_71 single-copy genes found in >= 80% of the MAGs, best hit per gene, concatenated
and aligned with muscle (anvi'o `anvi-get-sequences-for-hmm-hits`), trimmed with trimal
(-gt 0.5), maximum likelihood with IQ-TREE (ModelFinder, 1000 ultrafast bootstraps).
GTDB names of the MAGs come from `anvi-estimate-scg-taxonomy`; this job adds the tree.

    sbatch --output=$HOME/mag_phylogeny-%j.out mag_phylogeny.sbatch    # on blythe

Run from `/springbrook/.../runs/mag_phylogeny/`. On blythe (Oct 2026), jobs whose Slurm log is
written under /springbrook/share die at start with signal 53 and no log, so the log goes to $HOME.
Output: `mags_bac71.treefile` (+ `.iqtree`, `.contree`).
