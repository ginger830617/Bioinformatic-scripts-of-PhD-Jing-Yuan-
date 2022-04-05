#eSNP_Karotyping

#load relevant packages
#install_github("BenvenLab/eSNPKaryotyping/eSNPKaryotyping")
library("devtools")
library("zoo")
library("gplots")
library("eSNPKaryotyping")


############################ R in Linux ##############################

#Preparing Genome Indexes 
#bowtie2-build /home/u1566273/reference_genome/hg38/hg38.fa /home/u1566273/reference_genome/hg38/hg38

#Build the SAMtools genome index run: 
#samtools faidx /home/u1566273/reference_genome/hg38/hg38.fa

#Build the Picard genome dictionary file run: 
#java -jar /home/u1566273/program/picard.jar CreateSequenceDictionary R=/home/u1566273/reference_genome/hg38/hg38.fa O=/home/u1566273/reference_genome/hg38/hg38.dict

#The downloaded files have to be restructured 
#Edit_dbSNP_Files(Directory= "/home/u1566273/reference_genome/hg38/dbsnp_hg38/", File_Name="", Organism= "Human")

#Alignment of the reads to the reference
#Tophat(Directory = "/home/u1566273/esnp/FASTQ/”, Library_Type = "Paired", Threads = 8, Transcripts_Annotation = “/home/u1566273/reference_genome/hg38/hg38.gtf”, Bowtie_Genome_Index = “/home/u1566273/reference_genome/hg38/hg38” ) 

#Calling SNPs from the align BAM file 
#CreateVCF(Directory= "~/esnp/BAM/“, Genome_Fa= "~/reference_genome/hg38/hg38.fa", Picard_Path= "~/program/", GATK_Path= "~/program/")

#samtools view -H accepted_hits_rg_sorted.bam | grep '@RG'

#java -jar /home/u1566273/program/picard.jar MarkDuplicates I=/home/u1566273/esnp/BAM/accepted_hits_rg_sorted.bam O=dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics

#java -jar /home/u1566273/program/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/reference_genome/hg38/hg38.fa -I split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o output.vcf


################################### R in local ###############################

#Processing and filtering the variants table, try this in both sever and local
#This function creates file with SNPs data at the BAM directory called variantTable.csv 
#table=EditVCF(Directory = ,Organism = )
# Argument: 1. Directory - The Path of to the BAM Directory
#           2. Organism - "Human" or "Mouse

#table=EditVCF(Directory="/Volumes/eSNP/BAM/", Organism="Human")
#Next, run the following lines that sort the SNP table according to the chromosomal position: 

table= read.table("/Volumes/eSNP/K2T/BAM/variantTable.csv",header = TRUE)
table$chr=as.numeric(table$chr)
table=table[order(table$chr,table$position),]
table=table[table$chr>0,]


#Run the MajorMinorCalc function which filters the SNPs according to specified parameters (see Note 9) and calculates the allelic ratio of each SNP: 
#table2=MajorMinorCalc(Table, minDP, maxDP, minAF)
# Argument: 1. Table - The variable containing the table of SNPs
#           2. minDP - minimal reading depth required from accepted SNP, usually set to 20
#           3. maxDP - maximal reading depth required from accepted SNP, usually set to very high number
#           4. minAF - minimal minor allele frequency, usually set to 0.2-0.25

table2=MajorMinorCalc(Table=table, minDP=20, maxDP=10000, minAF=0.2)


# Plot Allelic ratio along the genome for duplication detection:
#PlotGenome(orderedTable = ,Window = ,Ylim = ,PValue = ,Organism = )
# Argument: 1. orderedTable - The variable containing the output of the MajorMinorCalc function
#           2. Window - an odd number determining the window of the moving average plot, usually 151
#           3. Ylim - the plot y axes maximal limit, usually 3
#           4. PValue - Determine if to add the P-Value bar to the graph, accepts "TRUE" or "FALSE" values
#           5. Organism - "Human" or "Mouse"

plot_true=PlotGenome(orderedTable=table2, Window=151, Ylim=3, PValue="TRUE", Organism="Human")


#Before the analysis, make sure that the common SNPs files were downloaded and edited as required (see Note 3), and then run the following function that intersects the observed SNPs with the list of the common SNPs from dbSNP database: 
#tbl=DeletionTable(Directory, Table, dbSNP_Data_Directory, dbSNP_File_Name, Genome_Fa_dict, Organism) 
#intersect the observed SNPs with the common SNPs table from dbSNPs, Creates file with the LOH data called             Deletions.txt
#Del_tbl=DeletionTable(Directory, Table, dbSNP_Data_Directory, dbSNP_File_Name, Genome_Fa_dict, Organism) 
# Argument: 1. Directory - The Path of to the BAM Directory, containing the variantTable.csv file
#           2. Table - The variable containing the output of the MajorMinorCalc function
#           3. dbSNP_Data_Directory - The path for the directory where the edited dbSNP file are (created previously by                    the Edit_dbSNP_Files function)
#           4. dbSNP_File_Name - The edited dbSNP file names, without the chromosome number
#           5. Genome_Fa_dict - the path for the dictionary file created for the whole genome FASTA file.
#           6. Organism - "Human" or "Mouse"

tbl=DeletionTable(Directory="/Volumes/eSNP/K2T/BAM/", Table=table2, dbSNP_Data_Directory="/Volumes/eSNP/reference/hg38/dbsnp_hg38/", dbSNP_File_Name= "Edited_Common_chr", Genome_Fa_dict = "/Volumes/eSNP/reference/hg38/hg38.dict", Organism="Human")


#Plot each SNP, without any summarization
#Plot_Zygosity_Sinle(Table = ,Organism = )
# Argument: 1. Table - The LOH table containing the output of the DeletionTable function
#           2. Organism - "Human" or "Mouse"

Plot_Zygosity_Sinle(Table =tbl, Organism = "Human")


#Plot blocks of heterozygous and homozygous SNPs
#Plot_Zygosity_Blocks(Table = ,Window = ,Max = ,Max2 = ,Organism = )
# Argument: 1. Table - The deletion table containing the output of the DeletionTable function
#           2. window - the block size in bp, usually 1500000
#           3. Max - How many Heterozygouse SNP need to be in a block to get the full color, usually 6
#           4. Max2 - How many Homozygouse SNP need to be in a block to get the full color, usually 60
#           5. Organism - "Human" or "Mouse"

Plot_Zygosity_Blocks(Table = tbl, Window = 1500000, Max = 6, Max2 = 60, Organism = "Human")

