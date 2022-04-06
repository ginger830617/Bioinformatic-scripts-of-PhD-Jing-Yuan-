###Differential expression of non-coding RNAs in UHCW###

rm(list=ls())
getwd()
setwd("/Volumes/LncRNA/DEseq_lncRNA/Analysis2")
setwd("/Users/u1566273/Projects/LncRNA-seq/DeSeq2/Analysis2")


####################################################################
# load libraries
####################################################################
# load the required functions
#library("abind", lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")
#install.packages("clipr","csvy","feather","fst","readODS","reshape2","rmatio","stringr","tidyr","xml2",lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")

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
library(gridExtra)
library(gtable)
library(cowplot)
library(grid)
library(dplyr)


#round count tables made by rsem
H2 = read.table("counting/WTCHG_109791_02.genes.rsem.counts",header = F, sep = "\t", row.names = 1)
head(H2)
H2 = round(H2)
write.table(H2, file ="WTCHG_109791_02.genes.rsem.round.counts",sep = "\t", quote = F, row.names = T, col.names = T)


T2 = read.table("counting/WTCHG_109791_04.genes.rsem.counts",header = F, sep = "\t", row.names = 1)
head(T2)
T2 = round(T2)
write.table(T2, file ="WTCHG_109791_04.genes.rsem.round.counts", sep = "\t", quote = F, row.names = T, col.names = T)


H3 = read.table("counting/WTCHG_110734_16.genes.rsem.counts",header = F, sep = "\t", row.names = 1)
H3 = round(H3)
write.table(H3, file = "WTCHG_110734_16.genes.rsem.round.counts", sep = "\t", quote = F, row.names = T, col.names = T)


T3 = read.table("counting/WTCHG_109791_05.genes.rsem.counts",header = F, sep = "\t", row.names = 1)
T3 = round(T3)
write.table(T3, file = "WTCHG_109791_05.genes.rsem.round.counts", sep = "\t", quote = F, row.names = T, col.names = T )


H4 = read.table("counting/WTCHG_110084_07.genes.rsem.counts", header = F, sep = "\t", row.names = 1 )
H4 = round(H4)
write.table(H4, file = "WTCHG_110084_07.genes.rsem.round.counts", sep = "\t", quote = F, row.names = T, col.names = T)

T4 = read.table("counting/WTCHG_109791_06.genes.rsem.counts", header = F, sep = "\t", row.names = 1)
T4 = round(T4)
write.table(T4, file = "WTCHG_109791_06.genes.rsem.round.counts", sep = "\t", quote = F, row.names = T, col.names = T)


H6 = read.table("counting/WTCHG_110734_18.genes.rsem.counts", header = F, sep = "\t", row.names = 1)
H6 = round(H6)
write.table(H6, file = "WTCHG_110734_18.genes.rsem.round.counts", sep = "\t", quote = F, row.names = T, col.names = T)


T6 = read.table("counting/WTCHG_110084_12.genes.rsem.counts",header = F, sep = "\t", row.names = 1)
T6 = round(T6)
write.table(T6, "WTCHG_110084_12.genes.rsem.round.counts", sep = "\t", quote = F, row.names = T, col.names = T)


H8 = read.table("counting/WTCHG_110084_13.genes.rsem.counts", header = F, sep = "\t", row.names = 1)
H8 = round(H8)
write.table(H8, "WTCHG_110084_13.genes.rsem.round.counts", sep = "\t", quote = F, row.names = T, col.names = T)


T8 = read.table("counting/WTCHG_110734_19.genes.rsem.counts", header = F, sep = "\t", row.names = 1)
T8 = round(T8)
write.table(T8, "WTCHG_110734_19.genes.rsem.round.counts", sep = "\t", quote = F, row.names = T, col.names = T)


H9 = read.table("counting/WTCHG_110084_14.genes.rsem.counts", header = F, sep = "\t", row.names = 1)
H9 = round(H9)
write.table(H9, "WTCHG_110084_14.genes.rsem.round.counts", sep = "\t", quote = F, row.names = T, col.names = T)


T9 = read.table("counting/WTCHG_110734_15.genes.rsem.counts", header = F, sep = "\t", row.names = 1)
T9 = round(T9)
write.table(T9, "WTCHG_110734_15.genes.rsem.round.counts", sep = "\t", quote = F, row.names = T, col.names = T)



#you must alter these paths to make them useful!
htseqDir = "counting"
sampleDataFilename = "sampletable.txt"

####################################################################
# Do DESeq2 -> create dds object
####################################################################
##  Read in the file identifying the different samples
sampleTable = read.table(sampleDataFilename,header=TRUE)
sampleTable



##  Read in the results from the LibiNorm analysis (the counts files)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory = htseqDir,design = ~ patientID + Tumour)
##  And perform the analysis (details in the manual)

##  And perform the analysis (details in the manual)
dds<-DESeq(ddsHTSeq)



####################################################################
# run a PCA
####################################################################
rld <- rlog(dds)

plotPCAEx(rld,1,2,cond = c("Tumour"), ntop=1000)
plotPCAEx(rld,2,3,cond = c("Tumour"), ntop=1000)
plotPCAEx(rld,1,3,cond = c("Tumour"), ntop=1000)
plotPCAEx(rld,1,1,cond = c("Tumour"), ntop=1000)

# calculate the variance for each gene
rv <- rowVars(assay(rld))

# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[1:700]


# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(rld)[select,]))

# analysis of PC1
hist(pca$rotation[,"PC1"])



####################################################################
# do the differential expression analysis 
####################################################################
res <- results(dds,contrast = c("Tumour", "Y", "N"),  alpha=0.001)
summary(res)
resSig = res[which(!is.na(res$padj)),]
res2 = resSig[resSig$padj < 0.001, ]
row.names(res2)
write.table(row.names(res2), file = "Differentially_expressed_gene_lncRNA.txt",quote = F, row.names = F, col.names = F)
res0rdered <- res2[order(res2$padj),]
row.names(head(res0rdered,20))
write.table(row.names(res0rdered), file = "Differentially_expressed_gene_lncRNA1.txt",quote = F, row.names = F, col.names = F)


write.table(res2_v, file = "res_DE_lncRNA.txt",quote = F, row.names = T, col.names = T)


#order the results table
#get the top 30 names
#get the counts table for all genes
#extract the correct rows
row.names(res0rdered)[seq(1,30)]
counts = counts(dds, normalized = T)
topGeneCounts = counts[row.names(counts) %in% row.names(res0rdered)[seq(1,30)], ]
write.table(topGeneCounts, file = "fpkm_lncRNA.txt",quote = F, row.names = T, col.names = T)

#all gene counts#
row.names(res0rdered)
counts = counts(dds, normalized = T)
allGeneCounts = counts[row.names(counts) %in% row.names(res0rdered), ]
write.table(allGeneCounts, file = "countAllgenes_lncRNA.txt",quote = F, row.names = T, col.names = T, sep = "\t")


#Remove extensions in gene IDs
lncRNA_9 = read.table("Differentially_expressed_gene_lncRNA.txt")
lncRNA_9 = lncRNA_9[,1]
lncRNA_9_ne = sub("\\.\\d*","",lncRNA_9)
write.table(lncRNA_9_ne, file = "Differentially_expressed_gene_lncRNA_ne.txt",quote = F, row.names = F, col.names = F)



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
pheatmap(log(counts+1),scale = "row",border_color = NA, cluster_rows = F)
pheatmap(log(counts+1),border_color = NA)
##################################################################

write.table(counts, file = "counts_p0001.txt",quote = F, row.names = T, col.names = T,sep = "\t")
##################################################################################################




###### Plots across the chromosome ###############################################################
merge_table <- read.table("merge_table.txt", header = T, sep = "\t")
head(merge_table)

##order the samples
order = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
          "chr10", "chr11", "chr12", "chr14", "chr15", "chr16", 
          "chr17", "chr18", "chr19", "chr20", "chr21", "chr22","chr23")
merge_table$seq_Name = factor(merge_table$seq_Name, levels=order)
rownames(merge_table) = paste0(merge_table$seq_Name,":",merge_table$start,"-",merge_table$end,"-",merge_table$gene_ID)
merge_table2 = merge_table[order(merge_table$seq_Name, merge_table$start),]

rownames(merge_table) = paste0(merge_table$seq_Name,":",merge_table$gene_ID)

head(merge_table2)
write.table(merge_table2, file = "merge_table2.txt",quote = F, row.names = T, col.names = T, sep = "\t")

######density plot#####################################
ggplot(merge_table2) + 
  geom_freqpoly(aes(start),binwidth = 50000000) +
  scale_x_continuous(breaks = seq(0,3e+8,1e+8)) +
  facet_grid(.~seq_Name,scales="free_x", space="free_x") + 
  ggtitle("DE_lncRNA") +
  xlab("Position across the genome") +  ylab("counts") +
  theme_bw()

ggsave("DE_lncRNA.pdf",limitsize = FALSE) # width = 20, height = 50)
########################################################


#Heatmap in each chromosome###############
q = rownames(merge_table2) = paste0(merge_table2$seq_Name,":",merge_table2$start,"-",merge_table2$end,"-",merge_table2$gene_ID)
pheatmap(log(merge_table2[,5:ncol(merge_table2)]+1),cluster_rows = F,border_color = NA,scale = "row", col = cm.colors(256))

pheatmap(log(merge_table2[merge_table2$seq_Name == "chr1",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr2",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr3",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr4",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr5",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr6",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr7",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr8",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr9",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr10",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr11",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr12",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr14",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr15",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr16",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr17",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr18",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr19",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr20",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr21",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr22",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")
pheatmap(log(merge_table2[merge_table2$seq_Name == "chr23",5:ncol(merge_table2)]+1),cluster_rows = F,scale = "row")



ggplot(merge_table) + 
  scale_x_continuous(breaks = seq(0,3e+8,1e+7)) +
  facet_grid(.~seq_Name,scales="free_x", space="free_x") + 
  #coord_cartesian(ylim = c(-100,100)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  ggtitle("Difference in expression in lncRNAs between Tumour and Control samples") +
  xlab("Position across the genome") + ylab("control - tumour") + 
  theme_bw() +
  #geom_smooth(aes(x=start,y=(Epi_1-WT), color = context),method = "loess") +
  geom_point(aes(x=start,y=(log(1+K2_N)-log(1+K2_Y))), alpha = 1, size = 1, color = "red") +
  geom_point(aes(x=start,y=(log(1+K3_N)-log(1+K3_Y))), alpha = 1, size = 1, color = "blue") +
  geom_point(aes(x=start,y=(log(1+K4_N)-log(1+K4_Y))), alpha = 1, size = 1, color = "yellow") +
  geom_point(aes(x=start,y=(log(1+K8_N)-log(1+K8_Y))), alpha = 1, size = 1, color = "purple")


  
##################Scatter plot for DE between T and N##############################
#Load table as new
merge_table <- read.table("merge_table.txt", header = T, sep = "\t")
merge_table <- merge_table[,-8]

order = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", 
          "chr10", "chr11", "chr12", "chr14", "chr15", "chr16", 
          "chr17", "chr18", "chr19", "chr20", "chr21", "chr22","chr23")
merge_table$seq_Name = factor(merge_table$seq_Name, levels=order)
head(merge_table)


merge_table$K2 = log(merge_table$K2_N +1)-log(merge_table$K2_Y+1)
merge_table$K3 = log(merge_table$K3_N +1)-log(merge_table$K3_Y+1)
merge_table$K4 = log(merge_table$K4_N +1)-log(merge_table$K4_Y+1)
merge_table$K8 = log(merge_table$K8_N +1)-log(merge_table$K8_Y+1)
head(merge_table)
y = gather(merge_table, "samples", "diffexpression", 13:16)
y = select(y, -c(K2_N, K3_N, K4_N, K8_N, K2_Y, K3_Y, K4_Y, K8_Y))
head(y)

ggplot(y) + 
  scale_x_continuous(breaks = seq(0,3e+8,1e+8)) +
  facet_grid(.~seq_Name,scales="free_x")+#, space="free_x") + 
  #coord_cartesian(ylim = c(-100,100)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  ggtitle("Difference in expression in lncRNAs between Tumour and Control samples") +
  xlab("Position across the genome") + ylab("control - tumour") + 
  theme_bw() +
  geom_point(aes(x=start,y=diffexpression, color = samples), alpha = 1, size = 1)

#####################################################################################




#######################################################################################################
lncRNA_gtf <- read.table("~/Projects/GFF_Intersector/Subset_lncRNA&mRNA/selected_lncRNAgenes.txt", header = F, sep = "\t")

colnames(lncRNA_gtf)
lncRNA_gtf_2 <- lncRNA_gtf[lncRNA_gtf$V3 %in% "gene", ]
nrow(lncRNA_gtf_2)

write.table(lncRNA_gtf_2, file = "lncRNA_gene.txt",quote = F, row.names = T, col.names = T,sep = "\t")

