###RSEM quantification of RNA-seq datasets of UHCW

###directory
/Volumes/shared210/Jing
sshfs u1566273@myfiles.warwick.ac.uk:  ~/myfiles
cd ~/myfiles/Shared210--FOHN/Jing/Marcos/2_trim


###FastQC raw reads
for i in *; do fastqc -o ../fastqc_1 $i & done
cd ~/RNA_data/fastqc_report
multiqc -p .

for variable in *; do echo $variable ; done
wc -l for i in *; do echo $i & done
wc –l *
  python ~/program/fastqc_analysis.py
more all_mod_scores.csv



###trim###
for i in WTCHG_110084_07 WTCHG_110084_12 WTCHG_110084_13 WTCHG_110084_14 WTCHG_110734_15 WTCHG_110734_16 WTCHG_110734_18 WTCHG_110734_19; do java -jar /home/u1566273/program/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ${i}_1.fastq.gz ${i}_2.fastq.gz  ../2_trim/${i}_1.fastq.gz ../2_trim/${i}_1UP.fastq.gz ../2_trim/${i}_2.fastq.gz ../2_trim/${i}_2UP.fastq.gz ILLUMINACLIP:/home/u1566273/program/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:24 TRAILING:24 SLIDINGWINDOW:4:20 MINLEN:40 ; done

java -jar /home/u1566273/program/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 WTCHG_109791_02_1.fastq.gz WTCHG_109791_02_2.fastq.gz  ../3_trim/WTCHG_109791_02_1.fastq.gz ../3_trim/WTCHG_109791_02_1UP.fastq.gz ../3_trim/WTCHG_109791_02_2.fastq.gz ../3_trim/WTCHG_109791_02_2UP.fastq.gz ILLUMINACLIP:/home/u1566273/program/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:24 TRAILING:24 SLIDINGWINDOW:4:20 MINLEN:40 

trim_galore --quality 20 --phred33 --stringency 3 --length 20 --paired \
--gzip --output_dir ../3_trim \
WTCHG_109791_02_1.fastq.gz WTCHG_109791_02_2.fastq.gz


####FastQC after trimming
for i in *; do fastqc -o ../fastqc_report $i & done
cd ~/RNA_data/fastqc_report
multiqc -p .

for variable in *; do echo $variable ; done
for i in *; do gunzip $i & done
wc -l for i in *; do echo $i & done
wc –l *
  for i in *; do fastqc -o ../2_fastqc $i & done
python ~/program/fastqc_analysis.py
more all_mod_scores.csv



#########lncRNA#############
Release 28 (GRCh38.p12)
Create Mapping Indices
Before we can perform NGS read mapping, we will create the genome indices using the genome FASTA file as input. You can re-use these indices in all your future RnA-seq mapping. However, if you wish to map to a different genome build/assembly, you have to re-run this step using different genome sequences and save the indices in a different directory.
Here, we will create indices by STAR and RSEM.

mkdir GENOME_data/star

STAR --runThreadN 40 --runMode genomeGenerate --genomeDir star \
--genomeFastaFiles gencode.v28.lncRNA_transcripts.fa \
--sjdbGTFfile gencode.v28.long_noncoding_RNAs.gtf


STAR --runThreadN 4 --runMode genomeGenerate --genomeDir /home/u1566273/reference_genome/Genoma_data/star_lncRNA --genomeFastaFiles GRCh38.d1.vd1.fa --sjdbGTFfile gencode.v28.long_noncoding_RNAs.gtf

mkdir GENOME_data/rsem

rsem-prepare-reference --gtf gencode.v28.long_noncoding_RNAs.gtf GRCh38.d1.vd1.fa rsem/GRCh38.lncRNA


####2PASS MODE
Mapping with STAR :This has be to done with each sample every time
STAR --genomeDir GENOME_data/star --readFilesCommand zcat \
--readFilesIn RNASEQ_data/GM12878.rep1.R1.fastq.gz RNASEQ_data/GM12878.rep1.R2.fastq.gz \
--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within \
--twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM \
--runThreadN 20 --outFileNamePrefix "RNASEQ_data/star_GM12878_rep1/"

for i in WTCHG_110084_12 WTCHG_110084_14 WTCHG_110734_15 WTCHG_110734_18; do STAR --genomeDir /home/u1566273/reference_genome/Genoma_data/star_lncRNA --readFilesCommand zcat --readFilesIn ${i}_1.fastq.gz ${i}_2.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM \
--runThreadN 40 --outFileNamePrefix /home/u1566273/RNASEQ_data/star2nd/${i} & done

for i in WTCHG_110084_07 WTCHG_110084_12 WTCHG_110084_13 WTCHG_110734_16; do STAR --genomeDir /home/u1566273/reference_genome/Genoma_data/star_lncRNA --readFilesCommand zcat --readFilesIn ${i}_1.fastq.gz ${i}_2.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM \
--runThreadN 20 --outFileNamePrefix /home/u1566273/RNASEQ_data/star2nd/${i} & done

for i in WTCHG_109791_02 WTCHG_109791_04 WTCHG_109791_05 WTCHG_109791_06; do STAR --genomeDir /home/u1566273/reference_genome/Genoma_data/star_lncRNA --readFilesCommand zcat --readFilesIn /home/u1566273/RNASEQ_data/trimmed/${i}_1.fastq.gz /home/u1566273/RNASEQ_data/trimmed/${i}_2.fastq.gz --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --outSAMunmapped Within --twopassMode Basic --outFilterMultimapNmax 1 --quantMode TranscriptomeSAM \
--runThreadN 20 --outFileNamePrefix /home/u1566273/RNASEQ_data/star2nd/${i} & done


###Simple statistics with samtools flagstat
samtools flagstat star_GM12878_rep1/Aligned.sortedByCoord.out.bam


###Quantification with RSEM
we use RSEM to quantify the expression of genes and transcript. 

for i in in WTCHG_109791_05 WTCHG_110734_16; do rsem-calculate-expression --bam --no-bam-output -p 10 --paired-end --forward-prob 0 \ 
/home/u1566273/RNASEQ_data/star2nd/${i}Aligned.toTranscriptome.out.bam \
/home/u1566273/reference_genome/Genoma_data/rsem/GRCh38.lncRNA /home/u1566273/RNASEQ_data/rsem_K3N/${i}rsem >& \
/home/u1566273/RNASEQ_data/rsem_K3N/${i}rsem.log & done

for i in WTCHG_109791_05; do rsem-calculate-expression --bam --no-bam-output -p 10 --paired-end --forward-prob 0 /home/u1566273/RNASEQ_data/star2nd/${i}Aligned.toTranscriptome.out.bam /home/u1566273/reference_genome/Genoma_data/rsem/GRCh38.lncRNA /home/u1566273/RNASEQ_data/rsem_K3T/${i}rsem >& /home/u1566273/RNASEQ_data/rsem_K3T/${i}rsem.log & done

for i in WTCHG_110734_18; do rsem-calculate-expression --bam --no-bam-output -p 10 --paired-end --forward-prob 0 /home/u1566273/RNASEQ_data/star2nd/${i}Aligned.toTranscriptome.out.bam /home/u1566273/reference_genome/Genoma_data/rsem/GRCh38.lncRNA /home/u1566273/RNASEQ_data/rsem_K6N/${i}rsem >& /home/u1566273/RNASEQ_data/rsem_K6N/${i}rsem.log & done

for i in WTCHG_110084_12; do rsem-calculate-expression --bam --no-bam-output -p 10 --paired-end --forward-prob 0 /home/u1566273/RNASEQ_data/star2nd/${i}Aligned.toTranscriptome.out.bam /home/u1566273/reference_genome/Genoma_data/rsem/GRCh38.lncRNA /home/u1566273/RNASEQ_data/rsem_K6T/${i}rsem >& /home/u1566273/RNASEQ_data/rsem_K6T/${i}rsem.log & done

for i in WTCHG_109791_02; do rsem-calculate-expression --bam --no-bam-output -p 10 --paired-end --forward-prob 0 /home/u1566273/RNASEQ_data/star2nd/${i}Aligned.toTranscriptome.out.bam /home/u1566273/reference_genome/Genoma_data/rsem/GRCh38.lncRNA /home/u1566273/RNASEQ_data/rsem_K2N/${i}rsem >& /home/u1566273/RNASEQ_data/rsem_K2N/${i}rsem.log & done

for i in WTCHG_109791_04; do rsem-calculate-expression --bam --no-bam-output -p 10 --paired-end --forward-prob 0 /home/u1566273/RNASEQ_data/star2nd/${i}Aligned.toTranscriptome.out.bam /home/u1566273/reference_genome/Genoma_data/rsem/GRCh38.lncRNA /home/u1566273/RNASEQ_data/rsem_K2T/${i}rsem >& /home/u1566273/RNASEQ_data/rsem_K2T/${i}rsem.log & done

for i in WTCHG_109791_06; do rsem-calculate-expression --bam --no-bam-output -p 10 --paired-end --forward-prob 0 /home/u1566273/RNASEQ_data/star2nd/${i}Aligned.toTranscriptome.out.bam /home/u1566273/reference_genome/Genoma_data/rsem/GRCh38.lncRNA /home/u1566273/RNASEQ_data/rsem_K4T/${i}rsem >& /home/u1566273/RNASEQ_data/rsem_K4T/${i}rsem.log & done

for i in WTCHG_110084_07; do rsem-calculate-expression --bam --no-bam-output -p 10 --paired-end --forward-prob 0 /home/u1566273/RNASEQ_data/star2nd/${i}Aligned.toTranscriptome.out.bam /home/u1566273/reference_genome/Genoma_data/rsem/GRCh38.lncRNA /home/u1566273/RNASEQ_data/rsem_K4N/${i}rsem >& /home/u1566273/RNASEQ_data/rsem_K4N/${i}rsem.log & done

for i in WTCHG_110734_19; do rsem-calculate-expression --bam --no-bam-output -p 10 --paired-end --forward-prob 0 /home/u1566273/RNASEQ_data/star2nd/${i}Aligned.toTranscriptome.out.bam /home/u1566273/reference_genome/Genoma_data/rsem/GRCh38.lncRNA /home/u1566273/RNASEQ_data/rsem_K8T/${i}rsem >& /home/u1566273/RNASEQ_data/rsem_K8T/${i}rsem.log & done

for i in WTCHG_110084_13; do rsem-calculate-expression --bam --no-bam-output -p 10 --paired-end --forward-prob 0 /home/u1566273/RNASEQ_data/star2nd/${i}Aligned.toTranscriptome.out.bam /home/u1566273/reference_genome/Genoma_data/rsem/GRCh38.lncRNA /home/u1566273/RNASEQ_data/rsem_K8N/${i}rsem >& /home/u1566273/RNASEQ_data/rsem_K8N/${i}rsem.log & done

for i in WTCHG_110084_14; do rsem-calculate-expression --bam --no-bam-output -p 10 --paired-end --forward-prob 0 /home/u1566273/RNASEQ_data/star2nd/${i}Aligned.toTranscriptome.out.bam /home/u1566273/reference_genome/Genoma_data/rsem/GRCh38.lncRNA /home/u1566273/RNASEQ_data/rsem_K9N/${i}rsem >& /home/u1566273/RNASEQ_data/rsem_K9N/${i}rsem.log & done

for i in WTCHG_110734_15; do rsem-calculate-expression --bam --no-bam-output -p 10 --paired-end --forward-prob 0 /home/u1566273/RNASEQ_data/star2nd/${i}Aligned.toTranscriptome.out.bam /home/u1566273/reference_genome/Genoma_data/rsem/GRCh38.lncRNA /home/u1566273/RNASEQ_data/rsem_K9T/${i}rsem >& /home/u1566273/RNASEQ_data/rsem_K9T/${i}rsem.log & done


###Prepare gene-level and transcript-level expression matrices
We use paste command to join the rsem.genes.results files side-by-side, then use cut to select the columns containing the expected_count information, and place them into a final output file. Repeat the same step for isoforms.
This one-line command assumes the genes (and transcripts) in each files are in the same order. If they are not, you will have to sort the files before joining them together.
cd ~/LSLNGS2015/
  
  paste RNASEQ_data/rsem_*/rsem.genes.results | tail -n+2 | \
cut -f1,5,12,19,26 > RNASEQ_data/edgeR.genes.rsem.txt

paste RNASEQ_data/rsem_*/rsem.isoforms.results | tail -n+2 | \
cut -f1,5,13,21,29 > RNASEQ_data/edgeR.isoforms.rsem.txt

wc -l RNASEQ_data/rsem_*/*results

paste /home/u1566273/RNASEQ_data/rsem_K2N/*rsem.genes.results | tail -n+2 | cut -f1,5 > /home/u1566273/RNASEQ_data/rsem_K2N/WTCHG_109791_02.genes.rsem.counts

paste /home/u1566273/RNASEQ_data/rsem_K3N/*rsem.genes.results | tail -n+2 | cut -f1,5 > /home/u1566273/RNASEQ_data/rsem_K3N/WTCHG_110734_16.genes.rsem.counts

paste /home/u1566273/RNASEQ_data/rsem_K4N/*rsem.genes.results | tail -n+2 | cut -f1,5 > /home/u1566273/RNASEQ_data/rsem_K4N/WTCHG_110084_07.genes.rsem.counts

paste /home/u1566273/RNASEQ_data/rsem_K6N/*rsem.genes.results | tail -n+2 | cut -f1,5 > /home/u1566273/RNASEQ_data/rsem_K6N/WTCHG_110734_18.genes.rsem.counts

paste /home/u1566273/RNASEQ_data/rsem_K8N/*rsem.genes.results | tail -n+2 | cut -f1,5 > /home/u1566273/RNASEQ_data/rsem_K8N/WTCHG_110084_13.genes.rsem.counts

paste /home/u1566273/RNASEQ_data/rsem_K2T/*rsem.genes.results | tail -n+2 | cut -f1,5 > /home/u1566273/RNASEQ_data/rsem_K2T/WTCHG_109791_04.genes.rsem.counts

paste /home/u1566273/RNASEQ_data/rsem_K3T/*rsem.genes.results | tail -n+2 | cut -f1,5 > /home/u1566273/RNASEQ_data/rsem_K3T/WTCHG_109791_05.genes.rsem.counts

paste /home/u1566273/RNASEQ_data/rsem_K4T/*rsem.genes.results | tail -n+2 | cut -f1,5 > /home/u1566273/RNASEQ_data/rsem_K4T/WTCHG_109791_06.genes.rsem.counts

paste /home/u1566273/RNASEQ_data/rsem_K8T/*rsem.genes.results | tail -n+2 | cut -f1,5 > /home/u1566273/RNASEQ_data/rsem_K8T/WTCHG_110734_19.genes.rsem.counts

paste /home/u1566273/RNASEQ_data/rsem_K6T/*rsem.genes.results | tail -n+2 | cut -f1,5 > /home/u1566273/RNASEQ_data/rsem_K6T/WTCHG_110084_12.genes.rsem.counts

paste /home/u1566273/RNASEQ_data/rsem_K9T/*rsem.genes.results | tail -n+2 | cut -f1,5 > /home/u1566273/RNASEQ_data/rsem_K9T/WTCHG_110734_15.genes.rsem.counts

paste /home/u1566273/RNASEQ_data/rsem_K9N/*rsem.genes.results | tail -n+2 | cut -f1,5 > /home/u1566273/RNASEQ_data/rsem_K9N/WTCHG_110084_14.genes.rsem.counts


###Create annotations in BED format
Below, we use a combination of commands to convert annotations recorded in the GTF file into BED format.

  grep -P "\tgene\t" gencode.v28.long_noncoding_RNAs.gtf | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
awk -F $'\t' 'BEGIN { OFS=FS } { print $1,$2-1,$3,$6,".",$4,$10,$12,$14 }' | \
sort -k1,1 -k2,2n > gencode.v28.long_noncoding_RNAs_gene.bed

grep -P "\ttranscript\t" gencode.v28.long_noncoding_RNAs.gtf | cut -f1,4,5,7,9 | \
sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
awk -F $'\t' 'BEGIN { OFS=FS } { print $1,$2-1,$3,$10,".",$4,$14,$16,$18 }' | \
sort -k1,1 -k2,2n > gencode.v28.long_noncoding_RNAs.transcript.bed  

grep -f my_gene_and_transcript_list.txt genes.gtf > selected_genes.gtf
grep -f 60kb_lncRNA.txt gencode.v28.long_noncoding_RNAs.gtf > selected_lncRNAgenes_60Kb.gtf


