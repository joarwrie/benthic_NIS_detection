########################
###  DEMULTIPLEXING  ###
########################

#!/bin/bash

# Requires cutadapt-2.8
# Values of e for 18S: e=0.16, for COI: e=0.12, for 16S: e=0.14

for i in $(ls *R1_001.fastq.gz);do
	Nom=$(echo $i | cut -f8 -d"/" | cut -f1 -d"_")
	reverse=$(basename $i "R1_001.fastq.gz")"R2_001.fastq.gz"
	cutadapt -g file:Demult_barcodes/barcode_COI_forward -G file:Demult_barcodes/barcode_COI_reverse -e 0.12 -o ${Nom}_{name1}_{name2}_lib2_R1_cut.fastq -p ${Nom}_{name1}_{name2}_lib2_R2_cut.fastq $i $reverse
done

# Select only files with proper tag combinations for further processing
mkdir inputDADA2
for i in $(ls *cut.fastq);do
	if [ $(echo $i | grep -c "_F1_R1_") -eq 1 ];then
		newname=$(echo $i | sed 's/_F1_R1_/_PCR1_/g')
		mv $i inputDADA2/$newname
	elif [ $(echo $i | grep -c "_F2_R2_") -eq 1 ];then
		newname=$(echo $i | sed 's/_F2_R2_/_PCR2_/g')
		mv $i inputDADA2/$newname
	elif [ $(echo $i | grep -c "_F3_R3_") -eq 1 ];then
		newname=$(echo $i | sed 's/_F3_R3_/_PCR3_/g')
		mv $i inputDADA2/$newname
	elif [ $(echo $i | grep -c "_F4_R4_") -eq 1 ];then
		newname=$(echo $i | sed 's/_F4_R4_/_PCR4_/g')
		mv $i inputDADA2/$newname
	elif [ $(echo $i | grep -c "_F5_R5_") -eq 1 ];then
		newname=$(echo $i | sed 's/_F5_R5_/_PCR5_/g')
		mv $i inputDADA2/$newname
	fi
done

#########################
###  READ PROCESSING  ###
#########################

# Requires R-3.5.1
Rscript DADA2.R

##############################
###  TAXONOMIC ASSIGNMENT  ###
##############################

# Requires blast-2.9.0
blastn -query List_ASVs_18S.fasta -db References/18S/BDD_18SFonseca_allMetazoa -out Blast_res_18S.txt -outfmt "6 qseqid qlen sseqid slen length gaps qcovs pident"

#############################
###  DOWNSTREAM ANALYSES  ###
#############################

Rscript Downstream_analyses.R

##################
###  GRAPHICS  ###
##################

Rscript Graphics.R

# END OF SCRIPT