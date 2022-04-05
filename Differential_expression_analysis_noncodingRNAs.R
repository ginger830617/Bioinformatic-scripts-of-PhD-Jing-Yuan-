#Differential expression analysis of UHCW

rm(list=ls())
setwd("/Users/u1566273/Projects/Marcos_datasets")

###########you must alter these paths to make them useful!##################
htseqDir = "/Users/u1566273/Projects/Marcos_datasets/htseq_counting"
sampleDataFilename = "sapmletable.txt"


####################################################################
# load libraries
####################################################################
# load the required functions

library(DESeq2)
library(baySeq)
source(file='~/bins/dataFunctions.r')
library(limma)
library(ggplot2)
library(gplots)
library(reshape2)
library(pheatmap)
library(genefilter)
library(VennDiagram)
library(gridExtra)
library(gtable)
library(cowplot)


####################################################################
# Do DESeq2 -> create dds object
####################################################################
##  Read in the file identifying the different samples
sampleTable = read.table("~/Projects/Marcos_datasets/sampletable.txt",header=TRUE)

##  Read in the results from the LibiNorm analysis (the counts files)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory = htseqDir,design = ~ patientID + Tumour )

##  And perform the analysis (details in the manual)
dds<-DESeq(ddsHTSeq)

#save.image("first_save.RData")
#load("first_save.RData")

####################################################################
# run a PCA
####################################################################

rld <- rlog(dds)
plotPCAEx(rld,1,2,cond = c("Tumour"))
plotPCAEx(rld,3,3,cond = c("Tumour"))
plotPCAEx(rld,2,2,cond = c("Tumour"))
plotPCAEx(rld,1,3,cond = c("Tumour"))

################calculate the variance for each gene###############
rv <- rowVars(assay(rld))

################select the ntop genes by variance##################
select <- order(rv, decreasing=TRUE)[1:500]


###perform a PCA on the data in assay(x) for the selected genes###
pca <- prcomp(t(assay(rld)[select,]))


# analysis of PC1
hist(pca$rotation[,"PC1"])

#plotCounts(dds,"ENSG00000117322",intgroup = c("patientID"))

write.table(names(pca$rotation[,1][pca$rotation[,1] > 0.07]),quote=F, row.names = F)
write.table(names(pca$rotation[,3][pca$rotation[,3] < -0.05]),quote=F, row.names = F)



####################################################################
#scatter plots####full####
####################################################################
Tvector = c(2,3,4,6,9,12)
Nvector = c(1,5,7,8,10,11)

par(mfrow = c(6,6))

for(i in Tvector) {
  for(j in Tvector){
    plot_xy_scatter(dds,i,j)
  }
}

for(i in Nvector) {
  for(j in Nvector){
    plot_xy_scatter(dds,i,j)
  }
}

for(i in Tvector){
  for(j in Nvector){
    plot_xy_scatter(dds,i,j)
  }
}




####################################################################
# do the differential expression analysis 
####################################################################

res <- results(dds, contrast=c("Tumour","y","n"),alpha=0.001)
summary(res)
resSig = res[which(!is.na(res$padj)),]
res2 = resSig[resSig$padj < 0.001, ]
row.names(res2)
#write.table(row.names(res2), file = "Differentially_expressed_gene.txt",quote = F, row.names = F, col.names = F)
res0rdered <- res2[order(res2$padj),]
row.names(head(res0rdered,20))


###################################################################
#generate a heatmap
###################################################################

all = row.names(res2)
mat = matrix(0,nrow = length(all), ncol = length(colnames(assay(dds))) )
row.names(mat) = all
colnames(mat) = colnames(assay(dds))

counts = counts(dds, normalized = T)
counts = counts[row.names(counts) %in% all ,]
pheatmap(log(counts+1),scale = "row",border_color = NA,annotation_legend = F,fontsize = 18,show_rownames = F)
pheatmap(log(counts+1),border_color = NA,annotation_legend = F,fontsize = 18,show_rownames = F)

#for(i in 1:nrow(mat)){
  #for(j in 1:length(colnames(assay(dds)))){
   # mat[i,j] = assay(dds)[row.names(assay(dds)) == row.names(mat)[i],j]
  #}
#}

#the matrix of dds is log and scale, which means the heatmap shows normalized gene expression values
#without scale , it shows absolute gene expression values. 
#scale could provide significant features
#pheatmap(log(mat+1),scale = "row",border_color = NA,annotation_legend = F,fontsize = 18,show_rownames = F)
#pheatmap(log(counts+1),border_color = NA)


##################################################################
#######################remove K9##################################
##################################################################
sampleTable2 = sampleTable[sampleTable$Samples != "K9_n",]
sampleTable3 = sampleTable2[sampleTable2$Samples != "K9_y",]
ddsHTSeq2<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable2,directory = htseqDir,design = ~Tumour )
ddsHTSeq3<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable3,directory = htseqDir,design = ~Tumour )
dds2<-DESeq(ddsHTSeq2)
dds3<-DESeq(ddsHTSeq3)

#
res <- results(dds, contrast=c("Tumour","y","n"),alpha=0.001)
summary(res)
resSig = res[which(!is.na(res$padj)),]
res2 = resSig[resSig$padj < 0.001, ]
row.names(res2)
write.table(row.names(res2), file = "Differentially_expressed_gene2.txt",quote = F, row.names = F, col.names = F)
res0rdered <- res2[order(res2$padj),]
row.names(head(res0rdered,20))

#
rld <- rlog(dds)
plotPCAEx(rld,1,2,cond = c("Tumour"))
plotPCAEx(rld,2,3,cond = c("Tumour"))
plotPCAEx(rld,1,3,cond = c("Tumour"))
plotPCAEx(rld,3,3,cond = c("Tumour"))
plotPCAEx(rld,2,2,cond = c("Tumour"))

#
all = row.names(res2)
mat = matrix(0,nrow = length(all), ncol = length(colnames(assay(dds))) )
row.names(mat) = all
colnames(mat) = colnames(assay(dds))
for(i in 1:nrow(mat)){
  for(j in 1:length(colnames(assay(dds)))){
    mat[i,j] = assay(dds)[row.names(assay(dds)) == row.names(mat)[i],j]
  }
}

pheatmap(log(mat+1),scale = "row",border_color = NA)

#
Tvector = c(2,3,4,6,10)
Nvector = c(1,5,7,8,9)

par(mfrow = c(5,5))

for(i in Tvector) {
  for(j in Tvector){
    plot_xy_scatter(dds,i,j)
  }
}

for(i in Nvector) {
  for(j in Nvector){
    plot_xy_scatter(dds,i,j)
  }
}

for(i in Tvector){
  for(j in Nvector){
    plot_xy_scatter(dds,i,j)
  }
}





#####################################################################
########################remove K6_y##################################
sampleTable4 = sampleTable3[sampleTable3$Samples != "K6_y",]
ddsHTSeq4<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable4,directory = htseqDir,design = ~Tumour )
dds4<-DESeq(ddsHTSeq4)

#
res <- results(dds4, contrast=c("Tumour","y","n"),alpha=0.001)
summary(res)
resSig = res[which(!is.na(res$padj)),]
res2 = resSig[resSig$padj < 0.001, ]
row.names(res2)
write.table(row.names(res2), file = "Differentially_expressed_gene3.txt",quote = F, row.names = F, col.names = F)
res0rdered <- res2[order(res2$padj),]
row.names(head(res0rdered,20))

#
rld <- rlog(dds)
plotPCAEx(rld,1,2,cond = c("Tumour"))
plotPCAEx(rld,2,3,cond = c("Tumour"))
plotPCAEx(rld,1,3,cond = c("Tumour"))
plotPCAEx(rld,3,3,cond = c("Tumour"))
plotPCAEx(rld,2,2,cond = c("Tumour"))

#
all = row.names(res2)
mat = matrix(0,nrow = length(all), ncol = length(colnames(assay(dds4))) )
row.names(mat) = all
colnames(mat) = colnames(assay(dds4))
counts = counts(dds4, normalized = T)
counts = counts[row.names(counts) %in% all ,]
pheatmap(log(counts+1),scale = "row",border_color = NA)
pheatmap(log(counts+1),border_color = NA)

#

all = row.names(res2)
mat = matrix(0,nrow = length(all), ncol = length(colnames(assay(dds))) )
row.names(mat) = all
colnames(mat) = colnames(assay(dds))
counts = counts(dds, normalized = T)
counts = counts[row.names(counts) %in% all ,]
pheatmap(log(counts+1),scale = "row",border_color = NA)
pheatmap(log(counts+1),border_color = NA)


#
Tvector = c(2,3,4,9)
Nvector = c(1,5,6,7,8)

par(mfrow = c(4,4))
for(i in Tvector) {
  for(j in Tvector){
    plot_xy_scatter(dds,i,j)
  }
}

par(mfrow = c(5,5))
for(i in Nvector) {
  for(j in Nvector){
    plot_xy_scatter(dds,i,j)
  }
}

par(mfrow = c(5,4))
for(i in Tvector){
  for(j in Nvector){
    plot_xy_scatter(dds,i,j)
  }
}


##################################################################


############# Make a basic volcano plot###########################
head(res)
# Make a basic volcano plot
with(res, plot(log2FoldChange, 
               -log10(pvalue), 
               main="Tumour vs Control", 
               xlim=c(-10,10),ylim=c(0,40)))
# Add colored points: red if padj<0.05, orange of log2FC>2, green if both)
with(subset(res, padj<0.05), 
     points(log2FoldChange, 
            -log10(pvalue), 
            pch=20, col="red"))

with(subset(res, padj>0.05), 
     points(log2FoldChange, -log10(pvalue), 
            pch=20, 
            col="orange"))

with(subset(res, padj<0.001 & abs(log2FoldChange)>2), 
     points(log2FoldChange, 
            -log10(pvalue), 
            pch=20, 
            col="green"))

#abline(v=log2(2),col='blue',lwd=0)
###########################################



###################################################################
#Remove extensions in gene IDs
###################################################################
full = read.table("/Volumes/Jing.Y3/GBM_12modified_codes/Reanalysis_12_database/Differentially_expressed_gene.txt")
K9_removed = read.table("/Volumes/Jing.Y3/GBM_12modified_codes/Reanalysis_12_database/Differentially_expressed_gene2.txt")
K6_y_removed = read.table("/Volumes/Jing.Y3/GBM_12modified_codes/Reanalysis_12_database/Differentially_expressed_gene3.txt")

a<- full
a<- full[,1]
a = sub("\\.\\d*","",a)
write.table(a, file = "Differentially_expressed_genefull.txt",quote = F, row.names = F, col.names = F)

c<- K9_removed
c<- K9_removed[,1]
c = sub("\\.\\d*","",c)
write.table(c, file = "Differentially_expressed_genefullK9re.txt",quote = F, row.names = F, col.names = F)

b<- K6_y_removed
b<- K6_y_removed[,1]
b = sub("\\.\\d*","",b)
write.table(b, file = "Differentially_expressed_genefullK6yre.txt",quote = F, row.names = F, col.names = F)
###################################################################






#####################################################################
#Plot Venn Diagram
###################################################################
full= read.table("/Volumes/Jing.Y3/GBM_12modified_codes/Reanalysis_12_database/Differentially_expressed_genefull.txt")
K9_removed = read.table("/Volumes/Jing.Y3/GBM_12modified_codes/Reanalysis_12_database/Differentially_expressed_geneK9re.txt")
K6_y_removed = read.table("/Volumes/Jing.Y3/GBM_12modified_codes/Reanalysis_12_database/Differentially_expressed_geneK6yre.txt")

a<- full[,1]
b<- K9_removed[,1]
c<- K6_y_removed[,1]

gene<-sort(unique(c(as.vector(a),as.vector(b),as.vector(c))))
counts<-matrix(0,nrow=length(gene),ncol=3)
for(i in 1:length(gene)){
  counts[i,1]<-as.vector(gene[i]%in%a)
  counts[i,2]<-as.vector(gene[i]%in%b)
  counts[i,3]<-as.vector(gene[i]%in%c)
  
}

colnames(counts)<-c("full","K9_removed","K6_y_removed")
cols<-c("Red","Green","Blue")
vc<-vennCounts(counts)

venn <- vennDiagram(object = vc,circle.col=cols)
###################################################################




####################################################################
#############compared to previous full dataset######################
full= read.table("/Volumes/Jing.Y3/GBM_12modified _codes/DeSeq2/Reanalysis_12_database/Differentially_expressed_genefull.txt")
full1 = read.table("/Volumes/Jing.Y3/NGBM/Deseq2/full/Differentially_expressed_gene.txt")
a<- full[,1]
b<- full1[,1]
gene<-sort(unique(c(as.vector(a),as.vector(b))))
counts<-matrix(0,nrow=length(gene),ncol=2)
for(i in 1:length(gene)){
  counts[i,1]<-as.vector(gene[i]%in%a)
  counts[i,2]<-as.vector(gene[i]%in%b)
  
}

colnames(counts)<-c("full","full1")
cols<-c("Red","Green")
vc<-vennCounts(counts)

venn <- vennDiagram(object = vc,circle.col=cols)
#################################################################






