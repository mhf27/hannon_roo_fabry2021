#!/bin/bash

#SBATCH -n 1 # one CPU
#SBATCH -N 1 # on one node
#SBATCH -t 0-2:00 # Running time of 2hr
#SBATCH --mem 20000 # Memory request
#SBATCH --mail-type=ALL
#SBATCH --mail-user=martin.fabry91@googlemail.com

source ~/software/MyVE/bin/activate


echo 'Start' ${1}
filename=${1%%.fq}
echo 'Start' $filename;
fastx_trimmer -Q33 -f 2 -l 49 -i $1 -o $filename.trimmed;
echo 'Trimmer done ' $filename;
#2 mismatches per read, add --alignIntronMax 1 to suppress spliced reads, add --alignEndsType EndToEnd to account for difficult to map indels (makes it more relatable to bowtie)
STAR --genomeDir /Users/fabry01/genomes/drosophila/dm6/UCSC/dm6/Sequence/STAR-Index/te_fused_clean_3 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1000 --outMultimapperOrder Random --outSAMmultNmax 1 --alignIntronMax 100 --alignEndsType EndToEnd --outFilterMismatchNoverLmax 0.04 --outReadsUnmapped Fastx --readFilesIn $filename.trimmed --runThreadN 6 --outFileNamePrefix $filename.te.;        
echo 'STAR dm6 done' $filename;
samtools index $filename.te.Aligned.sortedByCoord.out.bam;
echo 'Index bai done' $filename;
STAR --genomeDir /Users/fabry01/genomes/drosophila/dm6/UCSC/dm6/Sequence/STAR-Index/ --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outMultimapperOrder Random --outSAMmultNmax 1 --alignIntronMax 100 --alignEndsType EndToEnd --outFilterMismatchNoverLmax 0.04 --readFilesIn $filename.te.Unmapped.out.mate1 --runThreadN 6 --outFileNamePrefix $filename.dm6.;
samtools index $filename.dm6.Aligned.sortedByCoord.out.bam;
echo 'Index bai done' $filename;

wait

# 16 = antisense, 0 = sense transcripts mapping to TE genome
# Use sorted bam file as input
# Generates list with read number for each TE

#Split te reads into sense (0) and antisense (16)

samtools view -f 0x10 -b $filename.te.Aligned.sortedByCoord.out.bam > $filename.te.0.bam;
samtools index $filename.te.0.bam;
samtools idxstats $filename.te.0.bam | cut -f 1,3 > $filename.te.0.chrom_reads.txt;


samtools view -F 0x10 -b $filename.te.Aligned.sortedByCoord.out.bam > $filename.te.16.bam;
samtools index $filename.te.16.bam;
samtools idxstats $filename.te.16.bam | cut -f 1,3 > $filename.te.16.chrom_reads.txt;

#Split dm6 reads into sense (0) and antisense (16)

samtools view -f 0x10 -b $filename.dm6.Aligned.sortedByCoord.out.bam > $filename.dm6.0.bam;
samtools index $filename.dm6.0.bam;
samtools idxstats $filename.dm6.0.bam | cut -f 1,3 > $filename.dm6.0.chrom_reads.txt;


samtools view -F 0x10 -b $filename.dm6.Aligned.sortedByCoord.out.bam > $filename.dm6.16.bam;
samtools index $filename.dm6.16.bam;
samtools idxstats $filename.dm6.16.bam | cut -f 1,3 > $filename.dm6.16.chrom_reads.txt;


wait

# Normalisation with scaling factor!!!

echo 'Start Normalization'

#Calculate total reads
var1=$(samtools view -c -F 260 $filename.dm6.Aligned.sortedByCoord.out.bam)

var2=$(samtools view -c -F 260 $filename.te.Aligned.sortedByCoord.out.bam)

var3=$((var1+var2))


#calculate reads for condition (te.0) / total reads

var4=$(samtools view -c -F 260 $filename.te.0.bam)

var5=$(bc <<<"scale=10; $var4 / $var3")

echo $var5 > $filename.te.0.scaling.factor.txt

bamCoverage -b $filename.te.0.bam -bs 10 -p 4 --scaleFactor $var5 --normalizeUsing CPM -o $filename.te.0.cpm.norm.bw

wait

#calculate reads for condition (te.16) / total reads

var4=$(samtools view -c -F 260 $filename.te.16.bam)

var5=$(bc <<<"scale=10; $var4 / $var3")

echo $var5 > $filename.te.16.scaling.factor.txt

bamCoverage -b $filename.te.16.bam -bs 10 -p 4 --scaleFactor $var5 --normalizeUsing CPM -o $filename.te.16.cpm.norm.bw

wait

#calculate reads for condition (dm6.0) / total reads

var4=$(samtools view -c -F 260 $filename.dm6.0.bam)

var5=$(bc <<<"scale=10; $var4 / $var3")

echo $var5 > $filename.dm6.0.scaling.factor.txt

bamCoverage -b $filename.dm6.0.bam -bs 10 -p 4 --scaleFactor $var5 --normalizeUsing CPM -o $filename.dm6.0.cpm.norm.bw

wait

#calculate reads for condition (dm6.16) / total reads

var4=$(samtools view -c -F 260 $filename.dm6.16.bam)

var5=$(bc <<<"scale=10; $var4 / $var3")

echo $var5 > $filename.dm6.16.scaling.factor.txt

bamCoverage -b $filename.dm6.16.bam -bs 10 -p 4 --scaleFactor $var5 --normalizeUsing CPM -o $filename.dm6.16.cpm.norm.bw


wait

wait

#htseq on dm6.Ali

htseq-count -s reverse -f bam -i gene_name $filename.dm6.Aligned.sortedByCoord.out.bam /Users/fabry01/genomes/drosophila/dm6/UCSC/dm6/Sequence/GTF-files/ensembl_dm6_v2_rrna_dep.gtf > $filename.count.htseq

exit
