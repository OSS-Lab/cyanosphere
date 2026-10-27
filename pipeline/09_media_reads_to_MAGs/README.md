# Media reads against a shared MAG reference

Why: the media co-assemblies did not bin the dominant picocyanobacteria of V and NV
(MW101C3, *Cyanobium*): media MAGs hold 0% of the cyanobacterial coverage in V and
0-35% in NV, against 30-79% by RpS16. MAG-based attribution of function per medium
(Fig 4E) and functional trajectories need a reference that contains them. The 88
dereplicated enrichment MAGs do (`Cyanobium_sp016882315`, `MW101C3_sp002252635`,
`WH_5701`, ...), but not *Calothrix*, which is only in the NVNN co-assembly. So the
reference is the 88 enrichment MAGs plus the 181 media MAGs (completion >= 50%,
redundancy <= 10%), dereplicated at 95% ANI, and all 35 media samples (P0 parent,
P4, P5 parent line, P10) are mapped to it competitively.

## 1. Copy the genomes to blythe (from the laptop)

    bash pipeline/09_media_reads_to_MAGs/stage_genomes.sh

## 2. Run the job (on blythe)

    cd /springbrook/share/lifesci/lfsnhf/metagenomics/runs/media_reads_to_MAGs
    sbatch map_media.sbatch

It needs a conda env with dRep, bowtie2, samtools and CoverM; the sbatch creates
`media_mapping` from bioconda on first run. Output: `coverm_mean.tsv` and
`coverm_relative_abundance.tsv` (genome x sample), `drep/data_tables/Cdb.csv`
(cluster membership), `drep/data_tables/Wdb.csv` (representatives).

## 3. Bring the tables back

    scp blythe:/springbrook/share/lifesci/lfsnhf/metagenomics/runs/media_reads_to_MAGs/{coverm_mean.tsv,coverm_relative_abundance.tsv,drep/data_tables/Wdb.csv,drep/data_tables/Cdb.csv} data/media_mapping/     # from the repository root

`analysis/media.py` and `figures.py` then use them for Fig 4E (cyanobacterial share per
medium) and the functional trajectories.
