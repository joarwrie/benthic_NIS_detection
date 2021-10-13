##################
# PIPELINE DADA2
##################

# Requires R-3.5.1
# Requires dada2-1.13.1

library(dada2)
library(data.table)

# fastq files reading
path <- "./input_DADA2"

fnFs <- sort(list.files(path, pattern="_R1_cut.fastq", full.names=T))
fnRs <- sort(list.files(path, pattern="_R2_cut.fastq", full.names=T))
sample.names <- sapply(strsplit(basename(fnFs), "_R1_cut"), `[`, 1)

# Quality visualisation
#plotQualityProfile(fnFs[20:21])

# Filtering and cutting reads according to their quality
filtFs <- file.path(".", "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(".", "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
# truncLen values for each marker : 18S=260,180; COI=267,230; 16S=90,90
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,180), maxN=0, truncQ=0, multithread=T, rm.phix=F)

# Error modelling
errF <- learnErrors(filtFs, multithread=T)
errR <- learnErrors(filtRs, multithread=T)

# Data dereplication
derepF <- derepFastq(filtFs)
derepR <- derepFastq(filtRs)

# Denoising
dadaF <- dada(derepF, err=errF, multithread=T, pool="pseudo")
dadaR <- dada(derepR, err=errR, multithread=T, pool="pseudo")

# Merging read pairs (minOverlap values for each marker : 18S=60; COI=170; 16S=60)
mergers <- mergePairs(dadaF, derepF, dadaR, derepR, minOverlap=20, maxMismatch=0)

# Length selection (18=335:500; COI=303:323; 16S=105:125)
seqtab <- makeSequenceTable(mergers)
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 335:500]

# Chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=T, verbose=T)

# Exporting data files
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaF, getN), sapply(dadaR, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.table(file="tab_track_18S_dada2.csv", track, col.names=T, row.names=T, sep="\t", dec=".", quote=F)

ASVs=data.frame(Code=paste(rep("ASV",length(colnames(seqtab.nochim))), 1:length(colnames(seqtab.nochim)), sep="_"), Sequence=getSequences(seqtab.nochim))
colnames(seqtab.nochim)=paste(rep("ASV",length(colnames(seqtab.nochim))), 1:length(colnames(seqtab.nochim)), sep="_")
write.table(file="List_ASVs_18S.fasta", ASVs, col.names=F, row.names=F, sep="\n", dec=".", quote=F)
write.table(file="tab_distri_ASVs_18S.csv", seqtab.nochim, col.names=T, row.names=T, sep="\t", dec=".", quote=F)

