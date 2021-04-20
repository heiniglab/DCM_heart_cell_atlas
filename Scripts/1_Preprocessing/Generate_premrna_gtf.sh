# Remove annotation per exon
awk 'BEGIN{FS="\t"; OFS="\t"} $3 == "transcript"{ $3="exon"; print}' ./genes.gtf > ./genes_premrna.gtf

# create new STAR Index:
cellranger mkref --genome=refdata-cellranger-GRCh38-3.0.0_premrna  --fasta=/cellranger-3.0.2/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa --genes=./GRCh38-3.0.0_genes_premrna.gtf
