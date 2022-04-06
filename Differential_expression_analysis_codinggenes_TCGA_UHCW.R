#Differential expression analysis of combined RNA-seq databases from UHCW and TCGA

###############
rm(list=ls())
getwd()
setwd("/Volumes/Jing.Y4/177")
setwd("/Users/u1566273/Projects/TCGA_datasets")

####################################################################
# load libraries
####################################################################
# load the required functions
source(file='~/bins/dataFunctions.r')
library(DESeq2)
library(baySeq)
library(limma)
library(ggplot2)
library(gplots)
library(reshape2)
library(pheatmap)
library(genefilter)
library(ggrepel)


#you must alter these paths to make them useful!
htseqDir = "counting"
sampleDataFilename = "sampletable.txt"

####################################################################
# Do DESeq2 -> create dds object
####################################################################
##  Read in the file identifying the different samples
sampleTable = read.table(sampleDataFilename,header=TRUE)
sampleTable = sampleTable[-c(59,30,115,112,68,104,133,20,60,71,142,96,90,136,57,130,48,107,42,162,101,63,8,13,156),]
head(sampleTable)

##  Read in the results from the LibiNorm analysis (the counts files)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory = htseqDir,design = ~ Exp + Tumour)
##  And perform the analysis (details in the manual)

##  And perform the analysis (details in the manual)
dds<-DESeq(ddsHTSeq)

####################################################################
save.image("first_save.RData")
save.image("~/Projects/TCGA_datasets/12:11:2018:heatmap.RData")
load("first_save.RData")
load("loaded.RData")



####################################################################
# run a PCA
####################################################################
rld <- vst(dds)
rld2 = rld

plotPCAEx(rld2,1,2,cond = c("Exp", "Tumour"),ntop = 1000)
plotPCAEx(rld2,2,3,cond = c("Exp", "Tumour"),ntop = 1000)
plotPCAEx(rld2,1,3,cond = c("Exp", "Tumour"),ntop = 1000)
plotPCAEx(rld2,1,1,cond = c("Exp", "Tumour"),ntop = 1000)

# calculate the variance for each gene
rv <- rowVars(assay(rld2))

# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[1:700]


# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(rld2)[select,]))

# analysis of PC1
hist(pca$rotation[,"PC1"])

#tail(sub_sampletable,10)
#x = c(rep("old", 143), rep("new",9))

#make a suitable dataframe
PCx = 1
PCy = 10

pca.df =  c()
for (i in PCx:PCy){
  pca.coord = unclass(pca$x)[, c(i)]
  pca.number = rep(i,length = length(pca.coord),mode = numeric)
  names = row.names(unclass(pca$x))
  dff = data.frame(names,pca.number,pca.coord)
  pca.df = rbind(pca.df,dff)
  
}

pca.df$gender = sampleTable$Gender
pca.df$race = sampleTable$Race
pca.df$yob = sampleTable$Yearofbirth
pca.df$tumour = sampleTable$Tumour 
pca.df$exp = sampleTable$Exp
#pca.df$exp = x

pc1 = pca.df$pca.number = as.factor(pca.df$pca.number)
ggplot(pca.df, aes(x=pca.number, y=pca.coord,colour = gender,shape = tumour)) +
  geom_jitter(position=position_jitter(0.05)) +
  coord_flip() +
  ylab(label = "PCA Coordinate Value") +
  xlab(label = "Principle Component") +
  ggtitle("")+
  theme_bw()

ggplot(pca.df, aes(x=pca.number, y=pca.coord,colour = race,shape = tumour)) +
  geom_jitter(position=position_jitter(0.05)) +
  coord_flip() +
  ylab(label = "PCA Coordinate Value") +
  xlab(label = "Principle Component") +
  ggtitle("")+
  theme_bw()

ggplot(pca.df, aes(x=pca.number, y=pca.coord,colour = yob,shape = tumour)) +
  geom_jitter(position=position_jitter(0.05)) +
  coord_flip() +
  ylab(label = "PCA Coordinate Value") +
  xlab(label = "Principle Component") +
  ggtitle("")+
  theme_bw()

#
write.table(names(pca$rotation[,1][pca$rotation[,1] > 0.07]),quote=F, row.names = F)
write.table(names(pca$rotation[,3][pca$rotation[,3] < -0.05]),quote=F, row.names = F)
#

pca.df2 = c()
for (i in PCx:PCy){
  pca.coord = unclass(pca$x)[, c(i)]
  pca.number = rep(i,length = length(pca.coord),mode = numeric)
  names = row.names(unclass(pca$x))
  dff = data.frame(names,pca.number,pca.coord)
  
  
  pca.df2 = cbind(pca.df2,pca.coord)
}
colnames(pca.df2) = paste("PC",seq(PCx:PCy),sep="")
pca.df2 = as.data.frame(pca.df2)

pca.df2$gender = sampleTable$Gender
pca.df2$samples = sampleTable$PatientID
pca.df2$race = sampleTable$Race
pca.df2$yob = sampleTable$Yearofbirth
pca.df2$tumour = sampleTable$Tumour
pca.df2$exp = sampleTable$Exp
#pca.df$exp = x

PCx = 3
PCy = 4


ggplot(pca.df2, aes_string("PC1", "PC3",colour = "yob", shape = "tumour")) + geom_point(size=5) +
  coord_fixed()+geom_point() +
  theme_bw()#+
#geom_text_repel(aes(label = samples),
#                data = pca.df2,
#                color = "gray20",
#                force = 10)


ggplot(pca.df2, aes_string("PC1", "PC3",colour = "race",shape = "tumour")) + geom_point(size=5) +
  coord_fixed()+geom_point() +
  theme_bw()#+
#geom_text_repel(aes(label = samples),
#                data = pca.df2,
#                color = "gray20",
#                force = 10)

ggplot(pca.df2, aes_string("PC1", "PC3",colour = "gender",shape ="tumour")) + geom_point(size=5) +
  coord_fixed()+geom_point() +
  theme_bw()#+
#geom_text_repel(aes(label = samples),
#                data = pca.df2,
#                color = "gray20",
#                force = 10)



####################################################################
# do the differential expression analysis 
####################################################################
res <- results(dds,contrast = c("Tumour", "Y", "N"),  alpha=0.001)
summary(res)

resSig = res[which(!is.na(res$padj)),]
res2 = resSig[resSig$padj < 0.001, ]

#order the results table
#get the top 30 names
#get the counts table for all genes
#extract the correct rows
row.names(res0rdered)[seq(1,30)]
counts = counts(dds, normalized = T)
topGeneCounts = counts[row.names(counts) %in% row.names(res0rdered)[seq(1,30)], ]
write.table(topGeneCounts, file = "fpkm.txt",quote = F, row.names = T, col.names = T)

#all gene counts#
row.names(res0rdered)
counts = counts(dds, normalized = T)
allGeneCounts = counts[row.names(counts) %in% row.names(res0rdered), ]
write.table(allGeneCounts, file = "countAllgenes.txt",quote = F, row.names = T, col.names = T)


#################################################################################
#generate a heatmap
#get the differentially epressed gene list
#make a matrix with gene_IDs and samples
#give the row names with samples names
#give the coloumn names with gene_IDs
#give the read counts of genes in the database
#give the read counts of each differentially expressed genes in the gene list
#################################################################################
all = row.names(res2)
mat = matrix(0,nrow = length(all), ncol = length(colnames(assay(dds))) )
row.names(mat) = all
colnames(mat) = colnames(assay(dds))
counts = counts(dds, normalized = T)
counts = counts[row.names(counts) %in% all ,]
pheatmap(log(counts+1),scale = "row",border_color = NA)
pheatmap(log(counts+1),border_color = NA)

##################################################################


#################################################################
#cytoscape section#
#################################################################
setwd("/Users/u1566273/Desktop/gene_lists_new/heatmap3")

#get geneIDs cell cycle#
cell_cycle_list = read.table("/Volumes/Jing.Y4/gene_lists_new/Cell_cycle_list.txt")

#gene expression values (normalised)#
row.names(res0rdered)
counts = counts(dds, normalized = T)

#extract revelant rows from previous step#
counts_cellcycle = counts[sub("\\.\\d*","",row.names(counts)) %in% cell_cycle_list$V1,]

#generate heat map#
pdf("cellcycle_scale.pdf",height = 100, width = 35)
pheatmap(log(counts_cellcycle+1),scale = "row",border_color = NA)
dev.off()

pdf("cellcycle.pdf",height = 100, width = 35)
pheatmap(log(counts_cellcycle+1),border_color = NA)
dev.off()


####cytoskeleton####
cytoskeleton_list = read.table("/Volumes/Jing.Y4/gene_lists_new/Cytoskeleton_list.txt")
counts_cytoskeleton = counts[sub("\\.\\d*","",row.names(counts)) %in% cytoskeleton_list$V1,]

pdf("cytoskeleton_scale.pdf",height = 200, width = 35)
pheatmap(log(counts_cytoskeleton+1),scale = "row",border_color = NA)
dev.off()

pdf("cytoskeleton.pdf",height = 220, width = 35)
pheatmap(log(counts_cytoskeleton+1),border_color = NA)
dev.off()



####epigenetics####
epigenetics_list = read.table("/Users/u1566273/Desktop/gene_lists_new/Epigenetics_list.txt")
counts_epigenetics = counts[sub("\\.\\d*","",row.names(counts)) %in% epigenetics_list$V1,]

pdf("epigenetics_scale.pdf",height = 230, width = 60)
pheatmap(log(counts_epigenetics+1),scale = "row",border_color = NA)
dev.off()

pdf("epigenetics.pdf",height = 240, width = 65)
pheatmap(log(counts_epigenetics+1),border_color = NA)
dev.off()



####cell binding####
cell_binding_list = read.table("/Users/u1566273/Desktop/gene_lists_new/cell_binding_list.txt")
counts_cellbinding = counts[sub("\\.\\d*","",row.names(counts)) %in% cell_binding_list$V1,]

pdf("cellbinding_scale.pdf",height = 350, width = 60)
pheatmap(log(counts_cellbinding+1),scale = "row",border_color = NA)
dev.off()

pdf("cellbinding.pdf",height = 350, width = 65)
pheatmap(log(counts_cellbinding+1),border_color = NA)
dev.off()



####cancer pathway####
cancer_pathway_list = read.table("/Users/u1566273/Desktop/gene_lists_new/Cancer_pathway_list.txt")
counts_cancerpathway = counts[sub("\\.\\d*","",row.names(counts)) %in% cancer_pathway_list$V1,]

pdf("cancerpawthway_scale.pdf",height = 50, width = 35)
pheatmap(log(counts_cancerpathway+1),scale = "row",border_color = NA)
dev.off()

pdf("cancerpathway.pdf",height = 60, width = 50)
pheatmap(log(counts_cancerpathway+1),border_color = NA)
dev.off()



####DNA repair####
DNA_repair_list = read.table("/Users/u1566273/Desktop/gene_lists_new/DNA_repair_list.txt")
counts_DNArepair = counts[sub("\\.\\d*","",row.names(counts)) %in% DNA_repair_list$V1,]

pdf("DNArepair_scale.pdf",height = 50, width = 45)
pheatmap(log(counts_DNArepair+1),scale = "row",border_color = NA)
dev.off()

pdf("DNArepair.pdf",height = 55, width = 65)
pheatmap(log(counts_DNArepair+1),border_color = NA)
dev.off()



####Homebox####
homebox_list = read.table("/Users/u1566273/Desktop/gene_lists_new/Homebox_list.txt")
counts_homebox = counts[sub("\\.\\d*","",row.names(counts)) %in% homebox_list$V1,]

pdf("homebox_scale.pdf",height = 30, width = 80)
pheatmap(log(counts_homebox+1),scale = "row",border_color = NA)
dev.off()

pdf("homebox.pdf",height = 30, width = 70)
pheatmap(log(counts_homebox+1),border_color = NA)
dev.off()


####other####
other_list = read.table("/Users/u1566273/Desktop/gene_lists_new/other_list.txt")
counts_other = counts[sub("\\.\\d*","",row.names(counts)) %in% other_list$V1,]

pdf("other_scale.pdf",height = 55, width = 70)
pheatmap(log(counts_other+1),scale = "row",border_color = NA)
dev.off()

pdf("other.pdf",height = 55, width = 70)
pheatmap(log(counts_other+1),border_color = NA)
dev.off()


####plot normalized counts for a single gene####
plot_gene = function(ID){
  
  count = counts[row.names(counts) %in% ID,]
  names = c(rep("Y", 143), rep("N",5), rep("Y", 4))
  data = cbind(count, names)
  ggplot(as.data.frame(data)) +
    geom_boxplot(aes(x = names, y = as.numeric(as.character(count)))) + 
    theme_bw() 
  
  
}

#############################################################################

plot_gene("ENSG00000229425.1")
plot_gene("ENSG00000156970.11")
plot_gene("ENSG00000138778.10")
plot_gene("ENSG00000117724.11")
plot_gene("ENSG00000169679.13")
plot_gene("ENSG00000080986.11")
plot_gene("ENSG00000112742.8")
plot_gene("ENSG00000013810.17")
plot_gene("ENSG00000126787.11")
plot_gene("ENSG00000175063.15")
plot_gene("ENSG00000186185.12")
plot_gene("ENSG00000112742.8")
plot_gene("ENSG00000135476.10")


plot_gene("ENSG00000080986.11")
plot_gene("ENSG00000126787.11")
plot_gene("ENSG00000169679.13")
plot_gene("ENSG00000135476.10")
plot_gene("ENSG00000121621.6")
plot_gene("ENSG00000178999.11")


############Select candidate genes by fdr values##############

###cell_cycle###
cell_cycle_list [cell_cycle_list %in% row.names(res0rdered)]
row.names(res0rdered[sub("\\.\\d*", "", row.names(res0rdered)) %in% cell_cycle_list$V1,])
(res0rdered[sub("\\.\\d*", "", row.names(res0rdered)) %in% c("ENSG00000171241"),])


###other###
other_list [other_list %in% row.names(res0rdered)]
row.names(res0rdered[sub("\\.\\d*", "", row.names(res0rdered)) %in% other_list$V1,])
(res0rdered[sub("\\.\\d*", "", row.names(res0rdered)) %in% c("ENSG00000137825"),])


###homeobox###
homebox_list [homebox_list %in% row.names(res0rdered)]
row.names(res0rdered[sub("\\.\\d*", "", row.names(res0rdered)) %in% homebox_list$V1,])
(res0rdered[sub("\\.\\d*", "", row.names(res0rdered)) %in% c("ENSG00000137825"),])


###DNA repair###
DNA_repair_list [DNA_repair_list %in% row.names(res0rdered)]
row.names(res0rdered[sub("\\.\\d*", "", row.names(res0rdered)) %in% DNA_repair_list$V1,])
(res0rdered[sub("\\.\\d*", "", row.names(res0rdered)) %in% c("ENSG00000137825"),])


###cancer pathway###
cancer_pathway_list [cancer_pathway_list %in% row.names(res0rdered)]
row.names(res0rdered[sub("\\.\\d*", "", row.names(res0rdered)) %in% cancer_pathway_list$V1,])
(res0rdered[sub("\\.\\d*", "", row.names(res0rdered)) %in% c("ENSG00000137825"),])


###cytoskeleton###
cytoskeleton_list [cytoskeleton_list %in% row.names(res0rdered)]
row.names(res0rdered[sub("\\.\\d*", "", row.names(res0rdered)) %in% cytoskeleton_list$V1,])
(res0rdered[sub("\\.\\d*", "", row.names(res0rdered)) %in% c("ENSG00000137825"),])


###epigenetics###
epigenetics_list [epigenetics_list %in% row.names(res0rdered)]
row.names(res0rdered[sub("\\.\\d*", "", row.names(res0rdered)) %in% epigenetics_list$V1,])
(res0rdered[sub("\\.\\d*", "", row.names(res0rdered)) %in% c("ENSG00000123119"),])


###cell binding###
cell_binding_list [cell_binding_list %in% row.names(res0rdered)]
row.names(res0rdered[sub("\\.\\d*", "", row.names(res0rdered)) %in% cell_binding_list$V1,])
(res0rdered[sub("\\.\\d*", "", row.names(res0rdered)) %in% c("ENSG00000078902"),])

