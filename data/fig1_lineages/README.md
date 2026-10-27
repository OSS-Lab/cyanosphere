# Fig 1 lineages on the final RpL19 populations

Why: Fig 1B-D (end points vs reservoir) use the final RpL19 EcoPhylo run
(`03_ecophylo/RpL19_source_vs_enrichment_FINAL`), but the lineage time course (Fig 1E)
used an older RpL19 run (ECOPHYLO2) with its own clusters, so "populations" were not the
same unit. Here the 11 earlier-passage samples are profiled against the final run's RpL19
cluster references, with the final run's settings (references mode; bowtie2 --no-unal;
samtools view -F 4; anvi-profile --profile-SCVs --min-contig-length 0), so every Fig 1
number comes from one set of populations.

Samples: M_79-M_82 (DW S4-S7, P0), M_83-M_86 (P3), M_62 (lineage 2, P3), M_70 (P0),
M_76 (P0), copied from the lab NAS (sftp from blythe) into 01.RawData/<sample>/.
Control: M_132 (an end point already in the final run) from SRA (SRR37528228); its new
coverage must reproduce the existing column of data/rpl19_cluster_coverage.csv.

    sbatch fig1_lineages.sbatch      # on blythe; waits for the NAS copies to complete

Output: coverage_mean.txt and coverage_Q2Q3.txt (cluster split x sample), control.txt.
