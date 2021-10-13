##############################################
# Analyses of the data after read processing #
##############################################

# Contingency table modification
entree=as.data.frame(seqtab.nochim)
entree$Sample=row.names(entree)
entree=as.data.table(entree)
tableau=melt(entree, variable.name="ASV", value.name="Compte", id.vars="Sample")
tableau=tableau[tableau$Compte!=0,]
tableau$Echantillon=do.call(rbind, strsplit(as.character(tableau$Sample), "_"))[,1]
tableau$Replicat=do.call(rbind, strsplit(as.character(tableau$Sample), "_"))[,2]
tableau=tableau[tableau$Echantillon!="Undetermined",]
tableau$Port=tableau$Echantillon
tableau$Port=gsub("A0.","",tableau$Port)
tableau$Port=gsub("P0.rep.","",tableau$Port)
tableau$Saison="Automne"
tableau[grep("rep1", tableau$Echantillon),]$Saison="Printemps_1"
tableau[grep("rep2", tableau$Echantillon),]$Saison="Printemps_2"
tableau$Ponton=tableau$Echantillon
tableau$Ponton=gsub(".*A0", "", tableau$Ponton)
tableau$Ponton=gsub(".*P0", "", tableau$Ponton)
tableau$Ponton=gsub("rep.", "", tableau$Ponton)
tableau$Ponton=gsub("bis", "", tableau$Ponton)
tableau[tableau$Port=="TRbis",]$Port="TR"
tableau[tableau$Echantillon=="TRA01bis",]$Echantillon="TRA01"

# Calculation of the number of reads for each ASV or for each sample or both
tableau[,"TotASV":=sum(Compte), by=ASV]
tableau[,"TotEch":=sum(Compte), by=Echantillon]
tableau[,"TotASVEch":=sum(Compte), by=.(Echantillon,ASV)]
tableau[,"TotPonton":=sum(Compte), by=.(Port, Saison, Ponton)]

# Calculation of the percentage of each ASV in each sample
tableau$PourcentASV=tableau$TotASVEch/tableau$TotASV

# Identification of the maximum percentage of reads in a index control sample
neg=tableau[grep("^neg[0-9]", tableau$Echantillon),]
maxIndex=max(neg[neg$PourcentASV<1,]$PourcentASV)
# Correction for index jump (values used as maxIndex for our results : 18S=0.0182, COI=0.0109, 16S=0.0023)
tab_tagJump=tableau[tableau$PourcentASV>maxIndex,]
# Keep only ASVs that are present in at least 2 PCR replicates for each sample
tab_tagJump[,"RepCheck":=.N, by=.(Echantillon, ASV)]
tab_final=tab_tagJump[tab_tagJump$RepCheck>1,]

# Adding blast taxonomy
	# Preparing blast file
blast=read.csv("Blast_res_18S.txt", header=F, sep="\t", dec=".")
colnames(blast)=c("ASV", "qlen", "Ref", "Slen", "ali_length", "gaps", "qcov", "pident")
	# Keeping only alignments with at least 99% query cover or subject cover
blast$scov=blast$ali_length/blast$Slen*100
blastCov=blast[blast$qcov>99 | blast$scov>99,]
	# Finding the assignment with the highest identity percentage for each ASV
blastCov=as.data.table(blastCov)
blastCov[,"MaxIdent":=max(pident), by=ASV]
blastIdent=blastCov[blastCov$pident==blastCov$MaxIdent,]
	# Keeping only assignments where the identity percentage is above the selected threshold (18S=99%, COI=92%, 16S=97%)
blastIdent=blastIdent[blastIdent$MaxIdent>=99,]
	# Adding taxonomic info for each accession number
taxo=read.csv("References/18S/Taxo_BDD_18SFonseca_allMetazoa.txt", header=T, sep="\t", dec=".")
tab_taxo=merge(blastIdent, taxo, by.x="Ref", by.y="Accession.number", all.x=T, all.y=F)
colnames(tab_taxo)=c("Ref", "ASV", "Qlen", "Slen", "ali_length", "gaps", "qcov", "pident", "scov", "MaxIdent", "Species", "SPAcceptee", "Class", "Order", "Family")
tab_taxo$Genus=do.call(rbind, strsplit(as.character(tab_taxo$SPAcceptee), " "))[,1]
	# Collapsing taxonomic info if the ASV has several assignments possible
tab_taxo[,Assignation:=paste(unique(.SD[,SPAcceptee]),collapse="£"), by=ASV]
tab_taxo[,Assignation_genus:=paste(unique(.SD[,Genus]),collapse="£"), by=ASV]
tab_taxo[,Assignation_family:=paste(unique(.SD[,Family]),collapse="£"), by=ASV]
tab_taxo=unique(tab_taxo, by="ASV")
	# Going back to an assignment to the genus level or the family level if necessary. If one ASV was assigned to sequences from different families, the ASV is classified as unassigned.
tab_taxo[grep("£", tab_taxo$Assignation),]$Assignation=tab_taxo[grep("£", tab_taxo$Assignation),]$Assignation_genus
tab_taxo[grep("£", tab_taxo$Assignation),]$Assignation=tab_taxo[grep("£", tab_taxo$Assignation),]$Assignation_family
tab_taxo[grep("£", tab_taxo$Assignation),]$Assignation="Unassigned"
	# Adding ASV information
resultats=merge(tab_final, tab_taxo, by="ASV", all.x=T, all.y=F)

# Exporting final contingency tables
	# ASVs contingency table after processing
tab_sortie=dcast(tab_final, ASV~Sample, value.var="Compte", fun.aggregate=sum, fill=0)
write.table(file="Contingency_table_ASVs.csv", tab_sortie, col.name=T, row.name=F, sep="\t", dec=".", quote=F)
	# Taxa contingency table
tab_sortie2=dcast(resultats, Assignation~Sample, value.var="Compte", fun.aggregate=sum, fill=0)
write.table(file="Contingency_table_Taxa.csv", tab_sortie2, col.name=T, row.name=F, sep="\t", dec=".", quote=F)

# END OF SCRIPT
