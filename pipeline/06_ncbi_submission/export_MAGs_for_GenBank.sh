#!/bin/bash
# Export the MAGs from their anvi'o contigs databases for GenBank: FASTA,
# completeness, SCG taxonomy and assembly statistics. Run in the folder of the .db files.

set -e

# Create output directories
mkdir -p MAG_fasta
mkdir -p MAG_metadata

# Create external genomes file for batch operations
echo "Creating external genomes file..."
echo -e "name\tcontigs_db_path" > external_genomes.txt

for db in *.db; do
    # Skip the main MAGs-GENOMES.db file
    if [[ "$db" == "MAGs-GENOMES.db" ]]; then
        continue
    fi
    name="${db%.db}"
    echo -e "${name}\t$(pwd)/${db}" >> external_genomes.txt
done

echo "Found $(( $(wc -l < external_genomes.txt) - 1 )) MAG databases"

# Export FASTA for each MAG
echo ""
echo "Exporting FASTA files..."

for db in *.db; do
    if [[ "$db" == "MAGs-GENOMES.db" ]]; then
        continue
    fi
    name="${db%.db}"
    echo "  Exporting: $name"
    anvi-export-contigs -c "$db" \
        -o "MAG_fasta/${name}.fasta" \
        --no-wrap
done

echo "Exported $(ls MAG_fasta/*.fasta | wc -l) FASTA files"

# Get completeness and redundancy for all MAGs
echo ""
echo "Estimating genome completeness..."

anvi-estimate-genome-completeness \
    -e external_genomes.txt \
    -o MAG_metadata/completeness_report.txt

echo "Completeness report saved to MAG_metadata/completeness_report.txt"

# Get GTDB taxonomy for all MAGs
echo ""
echo "Estimating SCG taxonomy..."

anvi-estimate-scg-taxonomy \
    -e external_genomes.txt \
    -o MAG_metadata/taxonomy_report.txt

echo "Taxonomy report saved to MAG_metadata/taxonomy_report.txt"

# Summary statistics
echo ""
echo "Generating assembly stats..."

# Create a summary file with basic stats per MAG
echo -e "MAG_name\tnum_contigs\ttotal_length_bp\tN50" > MAG_metadata/assembly_stats.txt

for fasta in MAG_fasta/*.fasta; do
    name=$(basename "$fasta" .fasta)
    # Count contigs and total length
    num_contigs=$(grep -c "^>" "$fasta")
    total_length=$(grep -v "^>" "$fasta" | tr -d '\n' | wc -c)
    
    # Calculate N50
    lengths=$(grep -v "^>" "$fasta" | awk '{print length}' | sort -rn)
    half=$(echo "$total_length / 2" | bc)
    cumsum=0
    n50=0
    while IFS= read -r len; do
        cumsum=$((cumsum + len))
        if [ "$cumsum" -ge "$half" ]; then
            n50=$len
            break
        fi
    done <<< "$lengths"
    
    echo -e "${name}\t${num_contigs}\t${total_length}\t${n50}" >> MAG_metadata/assembly_stats.txt
done

echo "Done: MAG_fasta/, MAG_metadata/, external_genomes.txt"
