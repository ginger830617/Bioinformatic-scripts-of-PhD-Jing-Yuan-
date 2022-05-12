###RNA-seq analysis of UHCW###

###check the quality of sequences
Fastqc
for i in *; do fastqc -o ../2_fastqc $i & done
more all_mod_scores.csv
for i in *; do fastqc -o ../4_fastqc $i & done

Trim the sequences
for i in WTCHG_109791_02 WTCHG_109791_04 WTCHG_109791_05 WTCHG_109791_06 WTCHG_110084_07 WTCHG_110084_12 WTCHG_110084_14 WTCHG_110734_15 WTCHG_110734_16 WTCHG_110734_18; do java -jar /home/u1566273/program/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ${i}_1.fastq.gz ${i}_2.fastq.gz  ../2_trim/${i}_1.fastq.gz ../2_trim/${i}_1UP.fastq.gz ../2_trim/${i}_2.fastq.gz ../2_trim/${i}_2UP.fastq.gz ILLUMINACLIP:/home/u1566273/program/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:24 TRAILING:24 SLIDINGWINDOW:4:20 MINLEN:40 ; done


###Alignment
STEP 1: building the STAR index
~/Jing/bin/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /mnt/data/Jing/references/TCGA38 --genomeFastaFiles /mnt/data/Jing/references/TCGA38/GRCh38.d1.vd1.fa --sjdbGTFfile /mnt/data/Jing/references/TCGA38/gencode.v22.annotation.gtf --sjdbOverhang 100


STEP 2:Alignment 1st pass
##if work in shareholder, put absolute path in input files.

--readFilesIn /home/u1566273/Jing/GBM_38_12/new/3_trimmed/${i}*_1.fastq /home/u1566273/Jing/GBM_38_12/new/3_trimmed/${i}*_2.fastq
--runThreadN 8
--outFilterMultimapScoreRange 1
--outFilterMultimapNmax 20
--outFilterMismatchNmax 10
--alignIntronMax 500000
--alignMatesGapMax 1000000
--sjdbScore 2
--alignSJDBoverhangMin 1
--genomeLoad NoSharedMemory
--outFilterMatchNminOverLread 0.33
--outFilterScoreMinOverLread 0.33
--sjdbOverhang 100
--outSAMstrandField intronMotif
--outSAMtype None
--outSAMmode None
--outFileNamePrefix /home/u1566273/Jing/GBM_38_12/new/5_mapping/${i} & done

for i in WTCHG_109791_02 WTCHG_109791_04 WTCHG_109791_05 WTCHG_109791_06 WTCHG_110084_07 WTCHG_110084_12 WTCHG_110084_13 WTCHG_110084_14 WTCHG_110734_15 WTCHG_110734_16 WTCHG_110734_18 WTCHG_110734_19; do STAR --genomeDir /mnt/data/Jing/references/TCGA38 --readFilesIn ${i}*_1.fastq ${i}*_2.fastq --runThreadN 4 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif --outSAMtype None --outSAMmode None --outFileNamePrefix ${i} & done


STEP3:Intermediate index generation
STAR --runMode genomeGenerate --genomeDir /home/u1566273/reference_genome/Genoma_data/star2 --genomeFastaFiles /home/u1566273/reference_genome/Genoma_data/star/GRCh38.d1.vd1.fa --sjdbOverhang 100 --runThreadN 20 --sjdbFileChrStartEnd WTCHG_109791_02SJ.out.tab WTCHG_109791_04SJ.out.tab WTCHG_109791_05SJ.out.tab WTCHG_109791_06SJ.out.tab WTCHG_110084_07SJ.out.tab WTCHG_110084_12SJ.out.tab WTCHG_110084_13SJ.out.tab WTCHG_110084_14SJ.out.tab WTCHG_110734_15SJ.out.tab WTCHG_110734_16SJ.out.tab WTCHG_110734_18SJ.out.tab WTCHG_110734_19SJ.out.tab


STEP4: Alignment 2nd pass
--genomeDir <output_path from previous step>
  --readFilesIn <fastq_left_1>,<fastq_left2>,... <fastq_right_1>,<fastq_right_2>,...
--runThreadN <runThreadN>
  --outFilterMultimapScoreRange 1
--outFilterMultimapNmax 20
--outFilterMismatchNmax 10
--alignIntronMax 500000
--alignMatesGapMax 1000000
--sjdbScore 2
--alignSJDBoverhangMin 1
--genomeLoad NoSharedMemory
--limitBAMsortRAM 0
--readFilesCommand <bzcat|cat|zcat>
  --outFilterMatchNminOverLread 0.33
--outFilterScoreMinOverLread 0.33
--sjdbOverhang 100
--outSAMstrandField intronMotif
--outSAMattributes NH HI NM MD AS XS
--outSAMunmapped Within
--outSAMtype BAM SortedByCoordinate
--outSAMheaderHD @HD VN:1.4
--outSAMattrRGline <formatted RG line provided by wrapper>
  
  
for i in WTCHG_109791_05 WTCHG_109791_06 WTCHG_110084_07 WTCHG_110084_12 WTCHG_110084_13 WTCHG_110084_14 WTCHG_110734_15 WTCHG_110734_16 WTCHG_110734_18 WTCHG_110734_19; do STAR --genomeDir /home/u1566273/reference_genome/Genoma_data/star2 --readFilesCommand zcat --readFilesIn /home/u1566273/myfiles/Shared210--FOHN/Jing/Marcos/2_trim/${i}_1.fastq.gz /home/u1566273/myfiles/Shared210--FOHN/Jing/Marcos/2_trim/${i}_2.fastq.gz --runThreadN 10 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --limitBAMsortRAM 0 --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $i & done


STEP5: Sort bam files INDEX
samtools index <in.bam> [out.index]

for i in WTCHG_109791_02 WTCHG_109791_04 WTCHG_109791_05 WTCHG_109791_06 WTCHG_110084_07 WTCHG_110084_12 WTCHG_110084_13 WTCHG_110084_14 WTCHG_110734_15 WTCHG_110734_16 WTCHG_110734_18 WTCHG_110734_19; do samtools index ./$i*out.bam & done


STEP 6: Generate HTSeq-0.6.1p1
htseq-count \
-m intersection-nonempty \
-i gene_id \
-r name \
-s no \
- gencode.v22.annotation.gtf

Here, change “-s no” to “-s reverse”(default is “yes”, “reverse” means “yes” with reversed strand interpretation)


#install HTSeq#
tar -zxvf HTSeq-0.6.1.tar.gz
u1566273@vettel[HTSeq-0.6.1] python setup.py install --user


move sorted bam files to 6.3_sorting
for i in *;do python ~/program/HTSeq-0.6.1/scripts/htseq-count -m intersection-nonempty -f bam -i gene_id -r name -s reverse ${i} ~/reference_genome/Genoma_data/star2/gencode.v22.annotation.gtf > $i_htseq.counts & done

python ~/program/HTSeq-0.6.1/scripts/htseq-count -m intersection-nonempty -f bam -i gene_id -r name -s reverse WTCHG_109791_04Aligned.sortedByCoord.out.bam ~/reference_genome/Genoma_data/star2/gencode.v22.annotation.gtf > ./sort/WTCHG_109791_04_htseq.counts

python ~/program/HTSeq-0.6.1/scripts/htseq-count -m intersection-nonempty -f bam -i gene_id -r name -s reverse WTCHG_109791_06Aligned.sortedByCoord.out.bam ~/reference_genome/Genoma_data/star2/gencode.v22.annotation.gtf > ./sort/WTCHG_109791_06_htseq.counts

python ~/program/HTSeq-0.6.1/scripts/htseq-count -m intersection-nonempty -f bam -i gene_id -r name -s reverse WTCHG_110084_07Aligned.sortedByCoord.out.bam ~/reference_genome/Genoma_data/star2/gencode.v22.annotation.gtf > ./sort/WTCHG_110084_07_htseq.counts

python ~/program/HTSeq-0.6.1/scripts/htseq-count -m intersection-nonempty -f bam -i gene_id -r name -s reverse WTCHG_110084_13Aligned.sortedByCoord.out.bam ~/reference_genome/Genoma_data/star2/gencode.v22.annotation.gtf > ./sort/WTCHG_110084_13_htseq.counts

python ~/program/HTSeq-0.6.1/scripts/htseq-count -m intersection-nonempty -f bam -i gene_id -r name -s reverse WTCHG_110084_12Aligned.sortedByCoord.out.bam ~/reference_genome/Genoma_data/star2/gencode.v22.annotation.gtf > ./sort/WTCHG_110084_12_htseq.counts

python ~/program/HTSeq-0.6.1/scripts/htseq-count -m intersection-nonempty -f bam -i gene_id -r name -s reverse WTCHG_110734_15Aligned.sortedByCoord.out.bam ~/reference_genome/Genoma_data/star2/gencode.v22.annotation.gtf > ./sort/WTCHG_110734_15_htseq.counts

python ~/program/HTSeq-0.6.1/scripts/htseq-count -m intersection-nonempty -f bam -i gene_id -r name -s reverse WTCHG_110734_16Aligned.sortedByCoord.out.bam ~/reference_genome/Genoma_data/star2/gencode.v22.annotation.gtf > ./sort/WTCHG_110734_16_htseq.counts

python ~/program/HTSeq-0.6.1/scripts/htseq-count -m intersection-nonempty -f bam -i gene_id -r name -s reverse WTCHG_110734_18Aligned.sortedByCoord.out.bam ~/reference_genome/Genoma_data/star2/gencode.v22.annotation.gtf > ./sort/WTCHG_110734_18_htseq.counts

python ~/program/HTSeq-0.6.1/scripts/htseq-count -m intersection-nonempty -f bam -i gene_id -r name -s reverse WTCHG_110734_19Aligned.sortedByCoord.out.bam ~/reference_genome/Genoma_data/star2/gencode.v22.annotation.gtf > ./sort/WTCHG_110734_19_htseq.counts

python ~/program/HTSeq-0.6.1/scripts/htseq-count -m intersection-nonempty -f bam -i gene_id -r name -s reverse WTCHG_110084_14Aligned.sortedByCoord.out.bam ~/reference_genome/Genoma_data/star2/gencode.v22.annotation.gtf > ./sort/WTCHG_110084_14_htseq.counts


check the volume on roots
df .
du -sh 

