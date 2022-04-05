###TE transcript#######
#This package TEtranscripts estimates both gene and TE transcript abundances in RNA-seq data and conducts differential expression analysis on the resultant count tables. 


#sample table
Samples	fileName	patientID	Tumour
K2_n	WTCHG_109791_02Aligned.sortedByCoord.out.bam.htseq.counts	K2	n
K2_y	WTCHG_109791_04Aligned.sortedByCoord.out.bam.htseq.counts	K2	y
K3_y	WTCHG_109791_05Aligned.sortedByCoord.out.bam.htseq.counts	K3	y
K4_y	WTCHG_109791_06Aligned.sortedByCoord.out.bam.htseq.counts	K4	y
K4_n	WTCHG_110084_07Aligned.sortedByCoord.out.bam.htseq.counts	K4	n
K6_y	WTCHG_110084_12Aligned.sortedByCoord.out.bam.htseq.counts	K6	y
K8_n	WTCHG_110084_13Aligned.sortedByCoord.out.bam.htseq.counts	K8	n
K9_n	WTCHG_110084_14Aligned.sortedByCoord.out.bam.htseq.counts	K9	n
K9_y	WTCHG_110734_15Aligned.sortedByCoord.out.bam.htseq.counts	K9	y
K3_n	WTCHG_110734_16Aligned.sortedByCoord.out.bam.htseq.counts	K3	n
K6_n	WTCHG_110734_18Aligned.sortedByCoord.out.bam.htseq.counts	K6	n
K8_y	WTCHG_110734_19Aligned.sortedByCoord.out.bam.htseq.counts	K8	y


#TEtranscripts

usage: TEtranscripts -t treatment sample [treatment sample ...]
-c control sample [control sample ...]
--GTF genic-GTF-file
--TE TE-GTF-file
[optional arguments]

Required arguments:
  -t | --treatment [treatment sample 1 treatment sample 2...]
Sample files in group 1 (e.g. treatment/mutant), separated by space
-c | --control [control sample 1 control sample 2 ...]
Sample files in group 2 (e.g. control/wildtype), separated by space
--GTF genic-GTF-file  GTF file for gene annotations
--TE TE-GTF-file      GTF file for transposable element annotations

Optional arguments:
  
  *Input/Output options*
  --format [input file format]
Input file format: BAM or SAM. DEFAULT: BAM
--stranded [option]   Is this a stranded library? (yes, no, or reverse).
DEFAULT: yes.
--sortByPos           Input file is sorted by chromosome position.
--project [name]      Prefix used for output files (e.g. project name)
DEFAULT: TEtranscript_out


#1. Indexing (gene indexing and TE indexing)
inputs:
  alignment files(SAM/BAM)
Gene annotations (GTF)
TE annotations (GTF)

command lines:
  TEtranscripts --sortByPos --format BAM --mode multi -t RNAseq1.bam RNAseq2.bam -c CtlRNAseq1.bam CtlRNAseq.bam --project sample_sorted_test


TEtranscripts  --sortByPos --format BAM -t WTCHG_110734_19Aligned.sortedByCoord.out.bam -c WTCHG_110084_13Aligned.sortedByCoord.out.bam --GTF /home/u1566273/reference_genome/Genoma_data/gencode.v22.annotation.gtf --TE /home/u1566273/reference_genome/Genoma_data/hg38_rmsk_TE.gtf --project ../counts/TEtranscripts_k8


               
#2.Quantification (TE and Gene)
Paired-end (testdata_PE folder)

Input files: testData_treatment_rep1_PE.bam testData_treatment_rep2_PE.bam testData_control_rep1_PE.bam testData_control_rep2_PE.bam
GTF files (testdata_GTF folder): hg19_refGene.gtf, hg19_rmsk_TE.gtf
   # You will need to decompress them before use.
CMD line:
   TEtranscripts --sortByPos --mode multi --TE hg19_rmsk_TE.gtf --GTF hg19_refGene.gtf --project pairedEnd_test -t testData_treatment_rep1_PE.bam testData_treatment_rep2_PE.bam -c testData_control_rep1_PE.bam testData_control_rep2_PE.bam 2> log

The expected output files are:
  pairedEnd_test.cntTable
  pairedEnd_test_DESeq2.R
  pairedEnd_test_gene_TE_analysis.txt
  pairedEnd_test_sigdiff_gene_TE.txt


TEtranscripts  --sortByPos --format BAM -t WTCHG_110734_19Aligned.sortedByCoord.out.bam WTCHG_109791_04Aligned.sortedByCoord.out.bam WTCHG_109791_05Aligned.sortedByCoord.out.bam WTCHG_109791_06Aligned.sortedByCoord.out.bam -c WTCHG_110084_13Aligned.sortedByCoord.out.bam WTCHG_109791_02Aligned.sortedByCoord.out.bam WTCHG_110734_16Aligned.sortedByCoord.out.bam WTCHG_110084_07Aligned.sortedByCoord.out.bam WTCHG_110734_18Aligned.sortedByCoord.out.bam --GTF /home/u1566273/reference_genome/Genoma_data/gencode.v22.annotation.gtf --TE /home/u1566273/reference_genome/Genoma_data/hg38_rmsk_TE.gtf --mode multi -f 2 -p 0.001 --project /home/u1566273/RNASEQ_data/counts/nine_test_lp

TEcount  --sortByPos --format BAM --mode multi -b WTCHG_109791_02.bam --GTF /home/u1566273/reference_genome/Genoma_data/gencode.v22.annotation.gtf --TE /home/u1566273/reference_genome/Genoma_data/GRCh38_rmsk_TE.gtf --project /home/u1566273/RNASEQ_data/counts/nine_test_lp


#grep DE_TE list
grep -f my_gene_and_transcript_list.txt genes.gtf > selected_genes.gtf

grep -f 60kb_lncRNA.txt gencode.v28.long_noncoding_RNAs.gtf > selected_lncRNAgenes_60Kb.gtf
grep -f 60kb_mRNA.txt /Volumes/Jing.Y4/Backup/TCGA38\ GENOME\ HUMAN/gencode.v22.annotation.gtf > selected_mRNAgenes_60Kb.gtf
grep -f 60kb_mRNA.txt /Volumes/Jing.Y4/Backup/TCGA38\ GENOME\ HUMAN/gencode.v22.annotation.gtf > selected_mRNAgenes_60Kb.gtf


grep -w -F -f DE_TE.txt /home/u1566273/reference_genome/Genoma_data/hg38_rmsk_TE_1.gtf > selected_TE.gtf 
grep -f DE_TE.txt /home/u1566273/reference_genome/Genoma_data/hg38_rmsk_TE_1.gtf > selected_TE.txt
grep -f TE_DE1.txt /home/u1566273/reference_genome/Genoma_data/hg38_rmsk_TE_1.gtf > selected_TE1.txt


while read gene; do grep '$gene' /home/u1566273/reference_genome/Genoma_data/gencode.v22.annotation.gtf > /home/u1566273/RNASEQ_data/subset_mRNA.txt; done < /home/u1566273/RNASEQ_data/Differentially_expressed_gene1.txt

while read gene; do grep '$gene' /home/u1566273/reference_genome/Genoma_data/hg38_rmsk_TE.gtf ; done < /home/u1566273/RNASEQ_data/DE_TE.txt > subset_TE.txt

while read gene; do grep $gene /home/u1566273/reference_genome/Genoma_data/gencode.v28.long_noncoding_RNAs.gtf ; done < /home/u1566273/RNASEQ_data/Differentially_expressed_gene_lncRNA0.001.txt > subset_lncRNA.txt





