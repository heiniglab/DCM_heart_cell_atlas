#!/bin/sh
#Script for Chromium data processing fill in the parameters:
sample=$1
fastq=Sample_${sample}

export PATH=cellranger-3.0.2:$PATH
cellranger count \
--id=${sample}_premrna \
--fastqs=${fastq} \
--transcriptome=refdata-cellranger-GRCh38-3.0.0_premrna \
--sample=${sample} \
--expect-cells=5000 \
--localcores=12 \
--localmem=90 \
--uiport=3600
