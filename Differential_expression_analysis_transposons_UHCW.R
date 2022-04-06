#Differential expression analysis of transponsons from UHCW

setwd("~/Projects/Transponsons_analysis/twelve/")
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
library(cowplot)

#TEtranscripts Analysis#
setwd("~/Projects/Transponsons_analysis/eight/")



#load packages
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
library(cowplot)

###################################################
data <- read.table("~/Projects/Transponsons_analysis/eight/eight_test_lp.cntTable",header=T,row.names=1)
groups <- factor(c(rep("TGroup",4),rep("CGroup",4)))
min_read <- 1
data2 <- data[apply(data,1,function(x){max(x)}) > min_read,]
sampleInfo <- data.frame(groups,row.names=colnames(data2))

dds <- DESeqDataSetFromMatrix(countData = data2, colData = sampleInfo, design = ~ groups)
dds$condition = relevel(dds$groups,"CGroup")
dds <- DESeq(dds)

#################p005#########################
#res <- results(dds,independentFiltering=F)
res <- results(dds,contrast=c("groups", "TGroup", "CGroup"))
#resSig <- res[(!is.na(res$padj) & (res$padj < 0.05) & sqrt(res$log2FoldChange^2) >= 2), ]
resSig <- res[(!is.na(res$padj) & (res$padj < 0.05)),]
summary(resSig)

write.table(res, file="~/Projects/Transponsons_analysis/eight/eight_test_005_gene_TE_analysis.txt", sep="\t",quote=F)
write.table(resSig, file="~/Projects/Transponsons_analysis/eight/eight_test_005_sigdiff_gene_TE.txt",sep="\t", quote=F)
write.table(resSig[1279:1295,], file="~/Projects/Transponsons_analysis/eight/eight_test_005_sigdiff_TE.txt", sep="\t", quote = F)
resSig_1 <- read.table("~/Projects/Transponsons_analysis/eight/eight_test_005_sigdiff_TE.txt", header=T,row.names=1)


#########################################################
all = row.names(resSig_1)
mat = matrix(0,nrow = length(all), ncol = length(colnames(assay(dds))) )
row.names(mat) = all
colnames(mat) = colnames(assay(dds))
counts = counts(dds, normalized = T)
counts = counts[row.names(counts) %in% all ,]
pheatmap(log(counts+1),scale = "row",border_color = NA)
#pheatmap(log(counts+1),border_color = NA)

resSig_1 <- read.table("~/Projects/Transponsons_analysis/eight/eight_test_p005_sigdiff_TE.txt", header=T,row.names=1)
all3 = row.names(resSig_1[1:12,])
mat3 = matrix(0,nrow = length(all3), ncol = length(colnames(assay(dds))) )
row.names(mat3) = all3
colnames(mat3) = colnames(assay(dds))
counts3 = counts(dds, normalized = T)
counts3 = counts[row.names(counts) %in% all3,]
pheatmap(log(counts3+1),scale = "row",border_color = NA)


all4 = row.names(resSig_1[13:15,])
mat4 = matrix(0,nrow = length(all4), ncol = length(colnames(assay(dds))) )
row.names(mat4) = all4
colnames(mat4) = colnames(assay(dds))
counts4 = counts(dds, normalized = T)
counts4 = counts[row.names(counts) %in% all4,]
pheatmap(log(counts4+1),scale = "row",border_color = NA)


all5 = row.names(resSig_1[16:17,])
mat5 = matrix(0,nrow = length(all5), ncol = length(colnames(assay(dds))) )
row.names(mat5) = all5
colnames(mat5) = colnames(assay(dds))
counts5 = counts(dds, normalized = T)
counts5 = counts[row.names(counts) %in% all5,]
pheatmap(log(counts5+1),scale = "row",border_color = NA)



################P0.1#########################
#clean the history

data <- read.table("~/Projects/Transponsons_analysis/eight/eight_test_lp.cntTable",header=T,row.names=1)
groups <- factor(c(rep("TGroup",4),rep("CGroup",4)))
min_read <- 1
data2 <- data[apply(data,1,function(x){max(x)}) > min_read,]
sampleInfo <- data.frame(groups,row.names=colnames(data2))

dds <- DESeqDataSetFromMatrix(countData = data2, colData = sampleInfo, design = ~ groups)
dds$condition = relevel(dds$groups,"CGroup")
dds <- DESeq(dds)


#####################p01################################
#res <- results(dds,independentFiltering=F)
res <- results(dds,contrast=c("groups", "TGroup", "CGroup"))
resSig <- res[(!is.na(res$padj) & (res$padj < 0.1)),]
summary(resSig)

write.table(res, file="~/Projects/Transponsons_analysis/eight/eight_test_01_gene_TE_analysis.txt", sep="\t",quote=F)
write.table(resSig, file="~/Projects/Transponsons_analysis/eight/eight_test_01_sigdiff_gene_TE.txt",sep="\t", quote=F)
write.table(resSig[1839:1866,], file="~/Projects/Transponsons_analysis/eight/eight_test_01_sigdiff_TE.txt", sep="\t", quote = F)

write.table(resSig[1839:1866,], file="res_DE_TE_FDR0.1.txt", sep="\t", quote = F, row.names = T, col.names = T)

all = row.names(res)
mat = matrix(0,nrow = length(all), ncol = length(colnames(assay(dds))) )
row.names(mat) = all
colnames(mat) = colnames(assay(dds))
counts = counts(dds, normalized = T)
counts = counts[row.names(counts) %in% all,]
pheatmap(log(counts+1),scale = "row",border_color = NA, cluster_rows = F)
#pheatmap(log(counts+1),border_color = NA)

resSig_1 <- read.table("~/Projects/Transponsons_analysis/eight/eight_test_p01_TE.txt", header=T,row.names=1)
all1 = row.names(resSig_1)
mat1 = matrix(0,nrow = length(all1), ncol = length(colnames(assay(dds))) )
row.names(mat1) = all1
colnames(mat1) = colnames(assay(dds))
counts1 = counts(dds, normalized = T)
counts1 = counts[row.names(counts) %in% all1, ]
pheatmap(log(counts1+1),scale = "row",border_color = NA)

T_counts <- c(5582.034558,6969.909028,6304.236456,6218.9598096)
N_counts <- c(7247.095919,8508.85641,12065.80968,12227.32486)
log_counts <- c(log(T_counts+1), log(N_counts +1))
log_counts
groups <- c("T","T","T","T","N","N","N","N")
logcounts <- c(8.627488, 8.849501, 8.749136, 8.735519, 8.888494, 9.048980, 9.398214, 9.411510)
sample_ID <- c("K3Y", "K4Y", "K2Y","K8Y", "K3N", "K8N", "K2N","K4N")
MamTip2 <- data.frame(groups,logcounts,sample_ID)

ggplot(MamTip2, aes(x = groups, y = logcounts, colour = groups, fill= groups)) +
  scale_y_continuous(limits = c(5,10))+
  geom_boxplot(alpha=0.4) +
  theme(legend.position="none")

t.test(log(T_counts+1), log(N_counts +1) )

##########################################################
####plot normalized counts for a single gene####
plot_gene = function(ID){
  
  count = counts[row.names(counts) %in% ID,]
  names = c(rep("Y", 4), rep("N",4), rep("Y", 4))
  data = cbind(count, names)
  ggplot(as.data.frame(data)) +
    geom_boxplot(aes(x = names, y = as.numeric(as.character(count)))) + 
    theme_bw() 
  
  
}

##########################################################



plot_gene("MamTip2:hAT-Tip100:DNA")
########################################################################


all2 = row.names(resSig_1[1:17,])
mat2 = matrix(0,nrow = length(all2), ncol = length(colnames(assay(dds))) )
row.names(mat2) = all2
colnames(mat2) = colnames(assay(dds))
counts2 = counts(dds, normalized = T)
counts2 = counts[row.names(counts) %in% all2, ]
pheatmap(log(counts2+1),scale = "row",border_color = NA)


all3 = row.names(resSig_1[18:22,])
mat3 = matrix(0,nrow = length(all3), ncol = length(colnames(assay(dds))) )
row.names(mat3) = all3
colnames(mat3) = colnames(assay(dds))
counts3 = counts(dds, normalized = T)
counts3 = counts[row.names(counts) %in% all3,]
pheatmap(log(counts3+1),scale = "row",border_color = NA)


all4 = row.names(resSig_1[23:28,])
mat4 = matrix(0,nrow = length(all4), ncol = length(colnames(assay(dds))) )
row.names(mat4) = all4
colnames(mat4) = colnames(assay(dds))
counts4 = counts(dds, normalized = T)
counts4 = counts[row.names(counts) %in% all4,]
pheatmap(log(counts4+1),scale = "row",border_color = NA)




##################
rld <- rlog(dds)
rld2 = rld
a <- plotPCAEx(rld2,1,2,cond = c("groups"),ntop = 1000)
b <- plotPCAEx(rld2,2,3,cond = c("groups"),ntop = 1000)
c <- plotPCAEx(rld2,1,3,cond = c("groups"),ntop = 1000)

plot_grid(a,b,c, labels = c("A", "B","C"), ncol = 1, nrow = 3, rel_widths = c(1,1), rel_heights = c(1,1))


