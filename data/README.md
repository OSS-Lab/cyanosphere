# Inputs

Derived tables only (about 20 MB): reads are at NCBI, and the anvi'o databases these
were exported from are archived separately. `pipeline/` has the jobs that made
them.

| File | Rows | What it is | Produced by |
|---|---|---|---|
| `kegg_modules_enrichments.csv` | 5,734 | one row per KEGG module per sample, 8 enrichments and 8 reservoir samples: completeness, `avg_coverage`, `stepwise_copy_number`, and two legacy denominators (`bacterial_genomes_L19`, `N_MAGs`) | `anvi-estimate-metabolism --add-coverage --add-copy-number` on the per-sample contigs databases, columns subset |
| `kegg_modules_media.csv` | 4,320 | the same, for the 13 media samples that passed sequencing depth | as above |
| `rpl19_cluster_coverage.csv` | 659 | coverage of each 97% RpL19 sequence cluster in each of the 16 source and enrichment samples | anvi'o EcoPhylo (mean coverage, confirmed by the pipeline/10 control) |
| `fig1_lineages/coverage_mean-COVs.txt` | 659 | the same clusters in the 11 earlier-passage samples (M_62, M_70, M_76, M_79-M_86) plus the M_132 control | `pipeline/10_fig1_lineages` |
| `rpl19_cluster_taxonomy.tsv` | | GTDB taxonomy of each RpL19 cluster representative | EcoPhylo `06_MISC_DATA/Ribosomal_L19_scg_taxonomy_data.tsv` |
| `rpl19_cluster_tree.nwk` | | FastTree phylogeny of the RpL19 cluster representatives | EcoPhylo |
| `rpl19_phylum_coverage.tsv` | 25 | superseded, not read by the code: per-sample `anvi-estimate-scg-taxonomy` phylum matrix (drops unclassified RpL19 genes); kept for provenance. Phylum composition is now `taxonomy.phylum_table()` on the clusters | `anvi-estimate-scg-taxonomy --matrix-format` (no job script survives) |
| `rps16_population_counts.csv` | 30 | RpS16 clusters detected per sample, both experiments | EcoPhylo RpS16 runs, detection > 0.25 |
| `mag_module_pathwise_completeness.tsv` | 461 | pathwise completeness of every module in every dereplicated MAG | `anvi-estimate-metabolism` on the MAG collection |
| `mag_quality.tsv` | 88 | completeness, redundancy, MIMAG class and GTDB taxonomy of the dereplicated MAGs | `anvi-estimate-genome-completeness`, GTDB r214 |
| `mag_sample_coverage.tsv` | 88 | mean read coverage of each MAG in each sample | read recruitment to the dereplicated MAGs |
| `mag_bin_map.tsv` | 88 | dereplicated MAG to its source bin, by >90% bidirectional contig-header overlap | `pipeline/05_link_mags_to_bins.py` |
| `mag_biosample_reference.tsv` | 88 | genome size, contig count and collection metadata as submitted to NCBI | NCBI BioSample submission |
| `mag_biosample_accessions.tsv` | 88 | MAG name to BioSample accession | NCBI |
| `sra_accessions.csv` | 22 | sample name to BioSample and SRA run accession | NCBI, fetched 2026-08-17 |
| `symclatron_results.tsv` | 88 | symclatron lifestyle scores per MAG | `pipeline/04_symclatron/run.sbatch` |
| `ecophylo_coverage.csv` | 22,373 | Q2Q3 mean coverage and detection of every ribosomal-protein cluster in every sample, for four EcoPhylo runs (`run` column): RpS16 source vs enrichment (Supplementary Fig. S5), RpS16 media end-points (Fig. 4B), RpS16 and RpL19 media across all passages (Supplementary Figs. S12-S13) | `pipeline/07_export_raw_tables.py`, from the merged EcoPhylo `PROFILE.db` |
| `ecophylo_clusters.csv` | 1,150 | cluster size and GTDB taxonomy of each cluster representative, same four runs | as above |
| `assembly_stats.csv` | 16 | per-sample metaSPAdes assembly statistics (length, N50, contigs by size, genes, SCG hits) | as above, from `anvi-display-contigs-stats` exports |
| `media_samples.csv` | 42 | media-experiment sample IDs (M_ numbers) with condition, passage, replicate, lab description and, for the P10 end-points, the SRA sample name | as above, from the sequencing log |
| `raw_kegg/modules_unbinned_enrichments.csv.gz` | 5,734 | full `anvi-estimate-metabolism` output on the unbinned contigs of the 16 source and enrichment samples: enzyme hits, gene caller IDs, per-gene coverage and detection, copy numbers. `kegg_modules_enrichments.csv` is a column subset | `pipeline/07_export_raw_tables.py` |
| `raw_kegg/modules_unbinned_media.csv.gz` | 4,570 | the same for the 14 media samples (including NV_P10_R3, which the analysis excludes: only two bacterial RpS16 populations beside a eukaryotic phototroph, too few for a per-population metric) | as above, from the per-sample module files |
| `strong_bins_all.csv` | 511 | every STRONG bin before dereplication: coverage per sample, completeness, redundancy, length, GTDB taxonomy (SI Fig S2) | STRONG + anvi'o, per co-assembly |
| `media_mag_quality.csv`, `media_mag_coverage.csv`, `media_mag_modules.csv` | | MAGs of the three media co-assemblies (completion >= 50%, redundancy <= 10%): quality and taxonomy, coverage in every passage, per-MAG module completeness (Fig 4F) | `pipeline/08_media_mag_tables.py` |
| `schematics/` | 4 | vector crops of the hand-drawn panels (Fig 3A, 3C, 4A) and the anvi'o EcoPhylo heatmap (Fig 4B), placed into the composites by `figures.py` | Affinity Designer, anvi-interactive |

## "Population"

A population is one distinct 97% ribosomal-protein sequence cluster. The
count-based metric divides by the number of RpL19 clusters with non-zero
coverage in the sample (45 to 94 in the enrichments, 159 to 224 in the
reservoir), which is also the richness in Fig 1C. For the media experiment it
divides by RpS16 clusters, the marker used for composition in that experiment.
RpL19 gives the same cyanobacterial composition once one RpS16 cluster, a
photosynthetic-eukaryote plastid placed inside the cyanobacteria, is set aside
(the earlier belief that the picocyanobacteria lack RpL19 was this artefact).

## Legacy columns

`kegg_modules_*.csv` keep `bacterial_genomes_L19` and `N_MAGs`, the
denominators used until September 2026. `bacterial_genomes_L19` is an
assembly-level estimate whose exact rule was not recorded (it does not equal the
EcoPhylo cluster count); it is kept only for the robustness table
(`common.population_counts(..., "anvio_l19")`). The strain-level STRONG tables
and the species-collapsed RpL19 table used by earlier versions of the analysis
are not needed any more and are archived with the raw project data.
