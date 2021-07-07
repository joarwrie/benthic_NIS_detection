#!/bin/bash

# Values of e for 18S: e=0.16, for COI: e=0.12, for 16S: e=0.14
# Requires cutadapt-2.8

for i in $(ls *R1_001.fastq.gz);do
	Nom=$(echo $i | cut -f8 -d"/" | cut -f1 -d"_")
	reverse=$(basename $i "R1_001.fastq.gz")"R2_001.fastq.gz"
	cutadapt -g file:Demult_barcodes/barcode_COI_forward -G file:Demult_barcodes/barcode_COI_reverse -e 0.12 -o ${Nom}_{name1}_{name2}_lib2_R1_cut.fastq -p ${Nom}_{name1}_{name2}_lib2_R2_cut.fastq $i $reverse
done

