# benthic_NIS_detection
Pipeline used to evaluate the efficacy of water eDNA metabarcoding for the detection and monitoring of benthic NIS in marinas

This pipeline was used in the article “Erroneous detection and undetected species still hinder the routine monitoring of non-indigenous species by eDNA metabarcoding” by Marjorie Couton, Laurent Lévêque, Claire Daguin-Thiébaut, Thierry Comtet & Frédérique Viard, submitted for publication in Scientific reports.

The pipeline starts by demultiplexing the raw data with the file named "Demult_cutadapt.sh". This step requires the files present in the Demult_barcodes folder. Then the data can be processed using the DADA2.R script.
Raw sequencing data on which this pipeline was used are stored as SRA in ncbi under the BioProject ID: PRJNA744352.
