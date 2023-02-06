#############
# Libraries
#############

library(data.table)
library(ggplot2)
library(vegan)
library(RColorBrewer)
library(stringr)
library(ade4)

#################################
# Contingency table modification
#################################
# Import file
entree=read.csv("tab_distri_ASVs.csv", header = T, sep="\t")
entree$Sample=row.names(entree)
entree=as.data.table(entree)
# Transform into table
tableau=melt(entree, variable.name="ASV", value.name="Compte", id.vars="Sample")
tableau=tableau[tableau$Compte!=0,]
# Seperating sample info
tableau$Site=do.call(rbind, strsplit(as.character(tableau$Sample), "-"))[,1]
tableau$Replicat=do.call(rbind, strsplit(as.character(tableau$Sample), "_"))[,2]
tableau$Filter=do.call(rbind, strsplit(as.character(tableau$Sample), "_"))[,1]
tableau$Filter=do.call(rbind, strsplit(as.character(tableau$Filter), "-"))[,2]
tableau=tableau[tableau$Filter!="Undetermined",]

#################
# ASVs filtering
#################
# Calculation of the number of reads for each ASV or for each sample or both
tableau[,"TotASV":=sum(Compte), by=ASV]
tableau[,"TotEch":=sum(Compte), by=Site]
tableau[,"TotASVEch":=sum(Compte), by=.(Site,ASV)]
tableau[,"TotASVFilter":=sum(Compte), by=.(ASV, Site, Filter)]

# Calculation of the percentage of each ASV in each sample (filter unit)
tableau$PourcentASV=tableau$TotASVFilter/tableau$TotASV

# Identification of the maximum percentage of reads in a index control sample
neg=tableau[tableau$Site=="NA",]
maxIndex=max(neg$PourcentASV)

# Correction for index jump (maxIndex=0.0004)
tab_tagJump=tableau[tableau$PourcentASV>maxIndex,]

# Keep only ASVs that are present in at least 2 PCR replicates for each sample
tab_tagJump[,"RepCheck":=.N, by=.(Site, ASV)]
tab_final=tab_tagJump[tab_tagJump$RepCheck>1,]

# Production of Table S2
tab_tmp=as.data.table(table(unique(tab_final[,c(2,4)])$Site))
colnames(tab_tmp)=c("Site", "ASVs")
tab_final[,"TotReads":=sum(Compte), by=Site]
TableS2=merge(unique(tab_final[order(TotReads),c(4,13)]), tab_tmp, by="Site", all=T)
write.table(file="TableS2.txt", TableS2, col.names=T, row.names=F, sep="\t")

##########################################
# Checking for efficiency of the protocol
##########################################
# Accumulation curves
  # Per site including both filter and PCR replicates (20 replicates per site)
mat_asvs=reshape2::dcast(tab_final, Sample~ASV, fill=0, value.var="Compte", fun.aggregate = sum)
List_ech=mat_asvs$Sample
mat_asvs=as.matrix(mat_asvs[,-1])
row.names(mat_asvs)=List_ech
List_ech=do.call(rbind, strsplit(as.character(List_ech), "-"))[,1]
graph=vector("list", 20)
names(graph)=unique(List_ech)
compt=0
tab_graph=data.frame(Site="", Replicates="", Richness=0, sd=0)
for (i in unique(List_ech)){
  compt=compt+1
  graph[[compt]]=specaccum(mat_asvs[grep(i, row.names(mat_asvs)),])
  tab_tmp=cbind(rep(i, max(graph[[compt]]$sites)), graph[[compt]]$sites, graph[[compt]]$richness, graph[[compt]]$sd)
  colnames(tab_tmp)=c("Site", "Replicates", "Richness", "sd")
  tab_graph=rbind(tab_graph, tab_tmp)
}
tab_graph=tab_graph[-1,]
tab_graph$Replicates=as.integer(tab_graph$Replicates)
tab_graph$Richness=as.numeric(tab_graph$Richness)
tab_graph=as.data.table(tab_graph)
tab_graph[,"MaxReads":=max(Richness), by=Site]

p1=ggplot(tab_graph, aes(x=Replicates, y=Richness, group=Site)) +
  geom_line() +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25)) +
  xlab("\nNumber of replicates") + ylab("Number of ASVs\n") +
  geom_text(data=unique(tab_graph[,c(1,5)]), aes(x=21, y=MaxReads, label=Site))

pdf("Figure_2b.pdf", width=8, height=12)
p1
dev.off()

  # Per site displaying only the filter replicates (all 5 PCR replicates of one filter pooled)  
mat_asvs=tab_final
mat_asvs$filt_num=do.call(rbind, strsplit(as.character(mat_asvs$Sample), "_"))[,1]
mat_asvs=reshape2::dcast(mat_asvs, filt_num~ASV, fill=0, value.var="Compte", fun.aggregate = sum)
List_ech=mat_asvs$filt_num
mat_asvs=as.matrix(mat_asvs[,-1])
row.names(mat_asvs)=List_ech
List_ech=do.call(rbind, strsplit(as.character(List_ech), "-"))[,1]
graph=vector("list", 20)
names(graph)=unique(List_ech)
compt=0
tab_graph=data.frame(Site="", Replicates="", Richness=0, sd=0)
for (i in unique(List_ech)){
  compt=compt+1
  graph[[compt]]=specaccum(mat_asvs[grep(i, row.names(mat_asvs)),])
  tab_tmp=cbind(rep(i, max(graph[[compt]]$sites)), graph[[compt]]$sites, graph[[compt]]$richness, graph[[compt]]$sd)
  colnames(tab_tmp)=c("Site", "Replicates", "Richness", "sd")
  tab_graph=rbind(tab_graph, tab_tmp)
}
tab_graph=tab_graph[-1,]
tab_graph$Replicates=as.integer(tab_graph$Replicates)
tab_graph$Richness=as.numeric(tab_graph$Richness)
tab_graph=as.data.table(tab_graph)
tab_graph[,"MaxReads":=max(Richness), by=Site]

p1bis=ggplot(tab_graph, aes(x=Replicates, y=Richness, group=Site)) +
  geom_line() +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25)) +
  xlab("\nNumber of Filters") + ylab("Number of ASVs\n") +
  geom_text(data=unique(tab_graph[,c(1,5)]), aes(x=4.2, y=MaxReads, label=Site))

pdf("Figure_2c.pdf", width=8, height=12)
p1bis
dev.off()

  # Per Filter displaying the effect of the 5 PCR replicates
mat_asvs=reshape2::dcast(tab_final, Sample~ASV, fill=0, value.var="Compte", fun.aggregate = sum)
List_ech=mat_asvs$Sample
mat_asvs=as.matrix(mat_asvs[,-1])
row.names(mat_asvs)=List_ech
List_ech=do.call(rbind, strsplit(as.character(List_ech), "_"))[,1]
graph=vector("list", 80)
names(graph)=unique(List_ech)
compt=0
tab_graph=data.frame(Site="", Replicates="", Richness=0, sd=0)
for (i in unique(List_ech)){
  compt=compt+1
  graph[[compt]]=specaccum(mat_asvs[grep(i, row.names(mat_asvs)),])
  tab_tmp=cbind(rep(i, max(graph[[compt]]$sites)), graph[[compt]]$sites, graph[[compt]]$richness, graph[[compt]]$sd)
  colnames(tab_tmp)=c("Site", "Replicates", "Richness", "sd")
  tab_graph=rbind(tab_graph, tab_tmp)
}
tab_graph=tab_graph[-1,]
tab_graph$Replicates=as.integer(tab_graph$Replicates)
tab_graph$Richness=as.numeric(tab_graph$Richness)
tab_graph=as.data.table(tab_graph)
tab_graph[,"MaxReads":=max(Richness), by=Site]

p1ter=ggplot(tab_graph, aes(x=Replicates, y=Richness, group=Site)) +
  geom_line() +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  theme(axis.text=element_text(size=20), axis.title=element_text(size=25)) +
  xlab("\nNumber of PCR replicates") + ylab("Number of ASVs\n") +
  geom_text(data=unique(tab_graph[,c(1,5)]), aes(x=5.2, y=MaxReads, label=Site))

pdf("Figure_2d.pdf", width=8, height=12)
p1ter
dev.off()

# Rarefaction curves of reads for each site
mat_asvs=reshape2::dcast(tab_final, Site~ASV, fill=0, value.var="Compte", fun.aggregate = sum)
List_ech=mat_asvs$Site
mat_asvs=as.matrix(mat_asvs[,-1])
row.names(mat_asvs)=List_ech

pdf("Figure_2a.pdf", width=8, height=12)
  rarecurve(mat_asvs, step=100, ylab="ASVs", label=T)
dev.off()

####################################
# Taxonomical assignment with Blast
####################################
# Importing blast results
blast=read.csv("Blast_res.txt", header = F, sep = "\t", dec = ".")
colnames(blast)=c("ASV", "qlen", "GIAccNum", "slen", "Ali_length", "Gaps", "Qcovs", "Pident", "Taxid")
blast$AccNum=do.call(rbind, strsplit(as.character(blast$GIAccNum), "\\|"))[,4]
  # Keeping only alignments with at least 99% of query cover
blast=blast[blast$Qcovs>=99,]

# Importing taxonomy
# This file is a table with the full taxonomy for each accession number that matched with one of our sequences
taxo=rbind(read.csv("Tab_taxo_part1.txt", header = F, sep = "\t"), read.csv("Tab_taxo_part2.txt", header = F, sep = "\t"), read.csv("Tab_taxo_part3.txt", header = F, sep = "\t"))
colnames(taxo)=c("AccNum", "Family", "Species", "Phylum", "Taxid", "Genus", "Class", "Kingdom", "Order")
taxo=unique(taxo)

# Combining information
tab_taxo=merge(blast, taxo, by=c("AccNum", "Taxid"), all.x = T, all.y = F)
tab_taxo[grep("sp\\.", tab_taxo$Species), ]$Species="None"
# Elimination of unspecific information: sequences from environmental samples, from uncultured organisms or assigned to a higher rank than the family
tab_taxo=tab_taxo[grep("nvironmental", tab_taxo$Species, invert=T),]
tab_taxo=tab_taxo[grep("ncultured", tab_taxo$Species, invert=T),]
tab_taxo=tab_taxo[!(tab_taxo$Family=="None" & tab_taxo$Genus=="None" & tab_taxo$Species=="None"),]

# Choosing assignments for each ASV
tab_ass=as.data.table(tab_taxo)
# Calculating max identity for each ASV
tab_ass[,"MaxIdent":=max(Pident), by=ASV]
# Keeping all assignments at less than 1% below the max identity
tab_ass=tab_ass[tab_ass$Pident>=(tab_ass$MaxIdent-1),]
# Combining all assignments and finding the last common taxon
tab_ass[,"Assignment":=paste(unique(.SD[,Species]),collapse="£"), by=ASV]
tab_ass$level="Species"
tab_ass[Assignment%like%"£" | Assignment=="None", "level":="Genus"]
tab_ass[Assignment%like%"£" | Assignment=="None", "Assignment":=paste(unique(.SD[,Genus]),collapse="£"), by=ASV]
tab_ass[Assignment%like%"£" | Assignment=="None", "level":="Family"]
tab_ass[Assignment%like%"£" | Assignment=="None", "Assignment":=paste(unique(.SD[,Family]),collapse="£"), by=ASV]
tab_ass[Assignment%like%"£" | Assignment=="None", "level":="Order"]
tab_ass[Assignment%like%"£" | Assignment=="None", "Assignment":=paste(unique(.SD[,Order]),collapse="£"), by=ASV]
tab_ass[Assignment%like%"£" | Assignment=="None", "level":="Class"]
tab_ass[Assignment%like%"£" | Assignment=="None", "Assignment":=paste(unique(.SD[,Class]),collapse="£"), by=ASV]
tab_ass[Assignment%like%"£" | Assignment=="None", "level":="Phylum"]
tab_ass[Assignment%like%"£" | Assignment=="None", "Assignment":=paste(unique(.SD[,Phylum]),collapse="£"), by=ASV]
tab_ass[Assignment%like%"£" | Assignment=="None", "level":="Kingdom"]
tab_ass[Assignment%like%"£" | Assignment=="None", "Assignment":=paste(unique(.SD[,Kingdom]),collapse="£"), by=ASV]
Finale=unique(tab_ass, by="ASV")
Finale[grep("£", Finale$Assignment),]$Assignment="Unassigned"
Finale[Finale$Assignment=="Unassigned",]$level="Unassigned"
Finale[Finale$Assignment=="None",]$Assignment="Unassigned"

# Choosing thresholds
Finale=Finale[Finale$MaxIdent>=80,]
Finale[Finale$MaxIdent<90 & Finale$level=="Species",]$Assignment=do.call(rbind, strsplit(as.character(Finale[Finale$MaxIdent<90 & Finale$level=="Species",]$Assignment), " "))[,1]
Finale[Finale$MaxIdent<90 & Finale$level=="Species",]$level="Genus"

# Manual modifications for ASVs assigned to a higher rank than the species with 100% identity
Finale[Finale$ASV=="ASV_1277",]$Assignment="Niphargus fontanus"
Finale[Finale$ASV=="ASV_1277",]$level="Species"
Finale[Finale$ASV=="ASV_3113",]$Assignment="Niphargus fontanus"
Finale[Finale$ASV=="ASV_3113",]$level="Species"
Finale[Finale$ASV=="ASV_7070",]$Assignment="Tipula paludosa"
Finale[Finale$ASV=="ASV_7070",]$level="Species"
Finale[Finale$ASV=="ASV_1984",]$Assignment="Aporrectodea"
Finale[Finale$ASV=="ASV_1984",]$level="Genus"
Finale[Finale$ASV=="ASV_1331",]$Assignment="Phylloneta impressa"
Finale[Finale$ASV=="ASV_1331",]$level="Species"

# Linking ASVs information
tab_export=merge(Finale, tab_final, by="ASV", all.x=F, all.y=T)
tab_export[is.na(tab_export$Assignment),]$Assignment="Unassigned"
tab_export[tab_export$Assignment=="Unassigned",]$level="Unassigned"
  # Exporting the blast assignment results in a table
write.table(tab_export, file="Tab_assignment_blast_allGB_80.txt", col.names = T, row.names = F, sep="\t", dec=".", quote=F)

##############################
# Plotting assignment results
##############################
# Barplot assignment levels
tab_export[,"Reads":=sum(Compte), by=level]
tab_ASVs=unique(tab_export[,c(1,20,26,33)])
tab_ASVs$level=factor(tab_ASVs$level, levels=c("Unassigned", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
tab_ASVs[,"NbASVs":=.N, by=level]
Nbreads=unique(tab_ASVs[,c(2,4,5)])
Nbreads$X=Nbreads$level
tmp=Nbreads[Nbreads$X!="Unassigned" & Nbreads$X!="Kingdom",]
tmp$X="Kingdom"
Nbreads=rbind(Nbreads, tmp)
tmp=Nbreads[Nbreads$X!="Unassigned" & Nbreads$X!="Kingdom" & Nbreads$X!="Phylum",]
tmp$X="Phylum"
Nbreads=rbind(Nbreads, tmp)
tmp=Nbreads[Nbreads$X!="Unassigned" & Nbreads$X!="Kingdom" & Nbreads$X!="Phylum" & Nbreads$X!="Class",]
tmp$X="Class"
Nbreads=rbind(Nbreads, tmp)
tmp=Nbreads[Nbreads$X=="Species" | Nbreads$X=="Genus" | Nbreads$X=="Family",]
tmp$X="Order"
Nbreads=rbind(Nbreads, tmp)
tmp=Nbreads[Nbreads$X=="Species" | Nbreads$X=="Genus",]
tmp$X="Family"
Nbreads=rbind(Nbreads, tmp)
tmp=Nbreads[Nbreads$X=="Species",]
tmp$X="Genus"
Nbreads=rbind(Nbreads, tmp)
Nbreads[,"Somme":=sum(Reads), by=X]
Nbreads[,"TotASV":=sum(NbASVs), by=X]


p2=ggplot(Nbreads, aes(x=X, y=NbASVs, fill=level)) +
    geom_bar(stat="identity", position="stack", color="black") +
    theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
    theme(axis.text.x=element_text(size=14, angle=45, vjust=0.5), axis.text.y=element_text(size=16), axis.title=element_text(size=18, face="bold")) +
    xlab("") + ylab("Number of ASVs\n") +
    theme(legend.text = element_text(size=16), legend.title = element_blank()) +
    scale_fill_manual(values=c("gray60", brewer.pal(9, "Spectral")[c(1:3,5,7:9)])) +
    geom_label(data=unique(Nbreads, by="X"), aes(x=X, y=TotASV, label=Somme), fill="white", nudge_y = 350)

pdf("Figure3.pdf", width=9, height=6)
  p2
dev.off()

# Table S3 Contingency table: distribution of ASVs with their assignment
tabS3=dcast(tab_export, ASV~Site, fun.aggregate = sum, value.var = "Compte", fill=0)
tabS3=merge(tabS3, unique(tab_export[,c(1,11:20)]), by="ASV", all=T)
write.table(tabS3, file="Tab_distri_ASV_assigned_blast.txt", col.names = T, row.names = F, quote = F, sep = "\t", dec = ".")

# Table 1 ASVs assigned to metazoan species
tab1=tab_export[tab_export$Kingdom=="Metazoa",c(1,11:23)]
tab1[,"Reads":=sum(Compte), by=ASV]
tab1=unique(tab1, by=c("ASV", "Site"))
tab1[,"NbSites":=.N, by=ASV]
tab1=unique(tab1, by="ASV")
tab1=tab1[tab1$level=="Species" | tab1$level=="Genus" | tab1$level=="Family",]
write.table(tab1, file="Tab_assignment_metazoa.txt", col.names = T, row.names = F, quote = F, sep = "\t", dec = ".")

#####################
# Diversity analyses
#####################
# PCA based on hellinger transformation
conting=dcast(tab_final, ASV~Site, fun.aggregate = sum, fill=0, value.var = "Compte")
Noms=conting$ASV
conting=as.matrix(conting[,2:21])
row.names(conting)=Noms
tab_hellinger=decostand(t(conting), method="hellinger")
tab_pca=rda(tab_hellinger, scale=F)

tab_sc1=as.data.frame(scores(tab_pca, display = "sites", scaling=1))
tab_sc2=as.data.frame(scores(tab_pca, display = "species", scaling=1))
tab_sc2$length=sqrt(tab_sc2$PC1^2+tab_sc2$PC2^2)

tab_graph=rbind(tab_sc1, tab_sc2[tab_sc2$length>0.32465,1:2])
tab_graph$Type=c(rep("Site", 20), rep("ASV", 10))
tab_graph$aquifer=c("Unconsolidated", rep("Fissured", 3), rep("Unconsolidated", 3), rep("Fissured", 3), "Unconsolidated", rep("Fissured", 4), "Unconsolidated", rep("Fissured", 4), rep("None", 10))
tab_graph$Land=c(rep(c("Forest", "Field"), 2), rep("Forest", 4), "Field", rep("Forest", 3), "Field", "Field", "Forest", "Forest", "Field", "Forest", "Field", "Forest", rep("None", 10))

p3=ggplot(tab_graph[tab_graph$Type=="Site",], aes(x=PC1, y=PC2, color=aquifer, fill=aquifer, shape=Land)) +
  geom_vline(xintercept=0, color="black") +
  geom_hline(yintercept=0, color="black") + 
  geom_point(size=4) +
  scale_color_manual(values=c("Darkorange", "gray40", "Purple4")) +
  scale_fill_manual(values=c("Darkorange", "gray40", "Purple4")) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  theme(legend.position = "top", legend.title = element_text(size=20, face="bold")) +
  scale_shape_manual(values = c(24,21,22)) +
  xlab("\nPC1 (14.1%)") + ylab("PC2 (8.6%)\n") +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=20, face="bold"), legend.text = element_text(size=18)) +
  geom_segment(data=tab_graph[tab_graph$Type=="ASV",], aes(x=0, y=0, xend=PC1, yend=PC2), size=1, lineend = "round", linejoin = "round", arrow=arrow(length=unit(0.5, "cm"))) +
  geom_text(data=tab_graph, aes(x=PC1, y=PC2), label=row.names(tab_graph), nudge_x = -0.03, nudge_y = 0.02, color="black", size=5)

pdf("Figure6.pdf", width=8, height=7)
  p3
dev.off()

# Permanova
tab_env=data.frame(Sites=row.names(tab_hellinger), Aquifer=c("Unconsolidated", rep("Fissured", 3), rep("Unconsolidated", 3), rep("Fissured", 3), "Unconsolidated", rep("Fissured", 4), "Unconsolidated", rep("Fissured", 4)), Land=c(rep(c("Forest", "Field"), 2), rep("Forest", 4), "Field", rep("Forest", 3), "Field", "Field", "Forest", "Forest", "Field", "Forest", "Field", "Forest"))
adonis(tab_hellinger~Aquifer*Land, data=tab_env, permutations=9999, method="euclidean")

# Heatmap for all PCR replicates
conting=dcast(tab_final, ASV~Sample, fun.aggregate = sum, fill=0, value.var = "Compte")
Noms=conting$ASV
conting=as.matrix(conting[,-1])
row.names(conting)=Noms
tab_hellinger=decostand(t(conting), method="hellinger")
tab_graph=vegdist(tab_hellinger, method = "euclidean")
tab_graph=reshape2::melt(as.matrix(tab_graph), value.name = "Dist")
tab_graph$Site1=do.call(rbind, strsplit(as.character(tab_graph$Var1), "-"))[,1]
tab_graph$Site1=factor(tab_graph$Site1, levels=c("Ste1", "ReiL", "Sen", "Ror1", "Mul", "Ber", "Bru6", "Hina", "Unt", "Pfu", "ObeD", "Sam3", "BraN", "Sib3", "Pfa", "Gan", "Sem", "Bru8", "Abs", "GeiL"))
tab_graph=tab_graph[order(tab_graph$Site1),]
tab_graph$Var1=factor(tab_graph$Var1, levels=unique(tab_graph$Var1))
tab_graph$Var2=factor(tab_graph$Var2, levels=unique(tab_graph$Var1))

p4=ggplot(tab_graph, aes(x=Var1, y=Var2)) +
  geom_tile(aes(fill=Dist)) +
  scale_fill_continuous(high="khaki", low="black") +
  theme(axis.text.x=element_text(size=6, angle=45), axis.text.y=element_text(size=6), axis.title=element_text(size=18)) +
  xlab("") + ylab("") +
  theme(legend.text = element_text(size=16), legend.title = element_text(size=18))

pdf("Figure5", width=10, height=10)
  p4
dev.off()

########################
# Plotting Darn results
########################
# Importing Darn table (results of the tool DARN that performs phylogenetic placement)
DARN=read.csv("DARN_res_per_query.tsv", header=T, sep="\t", dec=".")
DARN=as.data.table(unique(DARN[,1:8]))
colnames(DARN)=c("ASV", "LWR", "fract", "aLWR", "afract", "Kingdom", "Phylum", "Class")

# Calculating max weight ratio for each ASV and keeping assignment with best LWR
DARN[,"MaxLWR":=max(LWR), by=ASV]
DARN=DARN[DARN$MaxLWR==DARN$LWR,]

# Dealing with double assignments:
DARN[,"cbKingdom":=paste(unique(.SD[,Kingdom]),collapse="£"), by=ASV]
DARN[,"cbPhylum":=paste(unique(.SD[,Phylum]),collapse="£"), by=ASV]
DARN[,"cbClass":=paste(unique(.SD[,Class]),collapse="£"), by=ASV]
DARN[grep("£", DARN$cbKingdom),]$cbKingdom="Unassigned"
DARN[grep("£", DARN$cbPhylum),]$cbPhylum="Unassigned"
DARN[grep("£", DARN$cbClass),]$cbClass="Unassigned"
DARN=unique(DARN[,-(3:9)])

# Applying a selection threshold of 75% for LWR
DARN[DARN$LWR<0.5,]$cbKingdom="Unassigned"
DARN[DARN$cbKingdom=="DISTANT",]$cbKingdom="Unassigned"
DARN[DARN$LWR<0.5,]$cbPhylum="Unassigned"
DARN[DARN$cbPhylum=="",]$cbPhylum="Unassigned"
DARN[DARN$LWR<0.5,]$cbClass="Unassigned"
DARN[DARN$cbClass=="",]$cbClass="Unassigned"

# Adding Kingdom info based on phylum (Since there are some conflict of taxonomy at the kingdom level in genbank, we upload here a correction file based on WoRMS taxonomy)
Taxo=read.csv("Taxo_correction.txt", header=T, sep="\t")
DARN=merge(DARN, Taxo, by="cbPhylum", all=T)
DARN[is.na(DARN$Kingdom),]$Kingdom=DARN[is.na(DARN$Kingdom),]$cbKingdom
DARN[is.na(DARN$Phylum),]$Phylum=DARN[is.na(DARN$Phylum),]$cbPhylum

# Keeping only ASVs without assignments
# The DARN tool was used on the whole list of ASVs but we look at the results only for ASVs that were not assigned with the blast method
List_ASVs=unique(tab_export[tab_export$level!="Species" & tab_export$level!="Genus" & tab_export$level!="Family",]$ASV)
DARN=DARN[DARN$ASV%in%List_ASVs,]
Manquants=length(List_ASVs[!(List_ASVs%in%DARN$ASV)])
Discarded=data.frame(cbPhylum=rep("Unassigned",Manquants), ASV=List_ASVs[!(List_ASVs%in%DARN$ASV)], LWR=rep(0,Manquants), cbKingdom=rep("Unassigned",Manquants), cbClass=rep("Unassigned",Manquants), Kingdom=rep("Unassigned",Manquants), Phylum=rep("Unassigned",Manquants))
DARN=rbind(DARN, Discarded)

# Merging both types of assignments and adding read numbers
DARN=DARN[,c(2,3,5:7)]
colnames(DARN)=c("ASV", "Pident", "Class", "Kingdom", "Phylum")
DARN$Type="DARN"
Blast=unique(tab_export[!(tab_export$ASV%in%List_ASVs),c(1,10,15,16,13)])
Blast$Type="Blast"
tab_graph=rbind(DARN, Blast)
tab_export[,"TotASVbis":=sum(Compte), by=ASV]
tab_graph=merge(tab_graph, unique(tab_export[,c(1,34)]), by="ASV", all.x=T, all.y=F)

# Taxonomic harmonisation
tab_graph[tab_graph$Phylum=="Discosea",]$Kingdom="Protozoa"
tab_graph[tab_graph$Phylum=="Discosea",]$Phylum="Amoebozoa"
tab_graph[tab_graph$Phylum=="Oomycota",]$Kingdom="Chromista"
tab_graph[tab_graph$Phylum=="Imbricatea",]$Kingdom="Chromista"
tab_graph[tab_graph$Phylum=="Imbricatea",]$Phylum="Cercozoa"
tab_graph[tab_graph$Phylum=="Rhodophyta",]$Kingdom="Plantae"
tab_graph[tab_graph$Phylum=="Endomyxa",]$Kingdom="Chromista"
tab_graph[tab_graph$Phylum=="Endomyxa",]$Phylum="Cercozoa"
tab_graph[tab_graph$Phylum=="Cercozoa",]$Kingdom="Chromista"
tab_graph[tab_graph$Phylum=="Tubulinea",]$Kingdom="Protozoa"
tab_graph[tab_graph$Phylum=="Tubulinea",]$Phylum="Amoebozoa"
tab_graph[tab_graph$Phylum=="Bacillariophyta",]$Kingdom="Chromista"
tab_graph[tab_graph$Phylum=="Proteobacteria",]$Kingdom="Bacteria"
tab_graph[tab_graph$Class=="Raphidophyceae" | tab_graph$Class=="Chrysophyceae" | tab_graph$Class=="Pelagophyceae" | tab_graph$Class=="Phaeophyceae" | tab_graph$Class=="Dictyochophyceae",]$Kingdom="Chromista"
tab_graph[tab_graph$Class=="Raphidophyceae" | tab_graph$Class=="Chrysophyceae" | tab_graph$Class=="Pelagophyceae" | tab_graph$Class=="Phaeophyceae" | tab_graph$Class=="Dictyochophyceae",]$Phylum="Ochrophyta"
tab_graph[tab_graph$Class=="Choanoflagellata",]$Kingdom="Protozoa"
tab_graph[tab_graph$Class=="Choanoflagellata",]$Phylum="Choanozoa"
tab_graph[tab_graph$Kingdom=="Eukaryota",]$Kingdom="Unassigned eukaryote"
tab_graph[tab_graph$Kingdom=="Viridiplantae",]$Kingdom="Plantae"

# Barplot for assignment distributions
tab_kingdom=tab_graph
tab_kingdom[,"SommeK":=sum(TotASVbis), by=.(Kingdom, Type)]
tab_kingdom[,"NbASVs":=.N, by=.(Kingdom, Type)]
tab_kingdom=unique(tab_kingdom[,c(4,6,8,9)])
tab_kingdom$Level="Kingdom"
tab_kingdom=rbind(tab_kingdom, list(Kingdom="Unassigned", Type="Blast", SommeK=sum(tab_kingdom[tab_kingdom$Type=="DARN",]$SommeK), NbASVs=4323, Level="Kingdom"))
tab_Phylum=tab_graph[tab_graph$Kingdom=="Metazoa",]
tab_Phylum[,"SommeK":=sum(TotASVbis), by=.(Phylum, Type)]
tab_Phylum[,"NbASVs":=.N, by=.(Phylum, Type)]
tab_Phylum=unique(tab_Phylum[,c(5,6,8,9)])
tab_Phylum$Level="Phylum"
colnames(tab_Phylum)[1]="Kingdom"
tab_graph=rbind(tab_kingdom, tab_Phylum)
tab_graph$Kingdom=factor(tab_graph$Kingdom, levels=c("Archaea","Bacteria","Unassigned eukaryote", "Protozoa", "Chromista", "Plantae", "Fungi", "Metazoa", "Other", "Unassigned", "Phylum", "Chordata", "Echinodermata", "Arthropoda", "Tardigrada", "Placozoa", "Mollusca", "Annelida", "Porifera", "Cnidaria", "Nemertea", "Nematoda", "Platyhelminthes", "Acanthocephala", "Gnathostomulida", "Gastrotricha", "Rotifera"))

p5=ggplot(tab_graph, aes(x=Level, y=SommeK, fill=Kingdom)) +
  geom_bar(color="black", stat="identity", position="fill") +
  facet_wrap(~Type, nrow=1) +
  theme(panel.background=element_rect(fill="white", colour="black"), axis.line=element_line(colour="black")) +
  theme(axis.text.x=element_text(size=18), axis.text.y=element_text(size=18), axis.title=element_text(size=20, face="bold")) +
  xlab("") + ylab("Proportion of reads (%)\n") +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), labels = c(0,25,50,75,100)) +
  theme(legend.text = element_text(size=18), legend.title = element_text(size=20, face="bold")) +
  theme(strip.background = element_rect(colour = "black", fill="white"), strip.text = element_text(size=20, face="bold")) +
  scale_fill_manual(values=c("lightsalmon", "indianred2", brewer.pal(8, "YlGnBu")[2:8], "gray50", "white", brewer.pal(11, "Set3")[1], "dodgerblue3", brewer.pal(11, "Set3")[6], "lightsalmon4", brewer.pal(11, "Set3")[3], "hotpink4", brewer.pal(11, "Set3")[7], "forestgreen", brewer.pal(11, "Set3")[2], "gold", brewer.pal(11, "Set3")[4], "firebrick", brewer.pal(11, "Set3")[c(8,10,11,5)]), drop=F)

pdf("Figure4.pdf", width=12, height=8)
  p5
dev.off()

#END
