# benthic_NIS_detection
Pipeline used to evaluate the efficacy of water eDNA metabarcoding for the detection and monitoring of benthic NIS in marinas

This pipeline was used in the article “Erroneous detection and undetected species still hinder the routine monitoring of non-indigenous species by eDNA metabarcoding” by Marjorie Couton, Laurent Lévêque, Claire Daguin-Thiébaut, Thierry Comtet & Frédérique Viard, submitted for publication in .

The pipeline can be performed by running the file named "file_to_run.sh". All softwares or R packages and their version number are indicated as comments. The parameters in this script are adapted to process 18S reads but the values for the COI reads and 16S reads are indicated as comments.

The references used to perform the taxonomic assignments for each marker can be found in the "References" folder, either in fasta format or in blast DB format.

Raw sequencing data on which this pipeline was used are stored as SRA in ncbi under the BioProject ID: PRJNA744352.
