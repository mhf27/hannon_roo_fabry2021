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
STAR --genomeDir /Users/fabry01/genomes/drosophila/dm6/UCSC/dm6/Sequence/STAR-Index/te_fused_clean_3 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1000 --outMultimapperOrder Random --outSAMmultNmax 1 --alignIntronMax 100 --alignEndsType EndToEnd --outFilterMismatchNoverLmax 0.04 --outReadsUnmapped Fastx --readFilesIn $filename.trimmed --limitBAMsortRAM 1300000000 --runThreadN 6 --outFileNamePrefix $filename.te.;        
echo 'STAR dm6 done' $filename;
samtools index $filename.te.Aligned.sortedByCoord.out.bam;
echo 'Index bai done' $filename;
STAR --genomeDir /Users/fabry01/genomes/drosophila/dm6/UCSC/dm6/Sequence/STAR-Index/ --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outMultimapperOrder Random --outSAMmultNmax 1 --alignIntronMax 100 --alignEndsType EndToEnd --outFilterMismatchNoverLmax 0.04 --readFilesIn $filename.te.Unmapped.out.mate1 --limitBAMsortRAM 1300000000 --runThreadN 6 --outFileNamePrefix $filename.dm6.;
samtools index $filename.dm6.Aligned.sortedByCoord.out.bam;
echo 'Index bai done' $filename;

wait



# Normalisation with scaling factor!!!

echo 'Start Normalization'

#Calculate total reads
var1=$(samtools view -c  $filename.dm6.Aligned.sortedByCoord.out.bam)

var2=$(samtools view -c  $filename.te.Aligned.sortedByCoord.out.bam)

var3=$((var1+var2))


#calculate reads for condition (te) / total reads

var4=$(samtools view -c $filename.te.Aligned.sortedByCoord.out.bam)

var5=$(bc <<<"scale=10; $var4 / $var3")

echo $var5 > $filename.te.scaling.factor.txt

bamCoverage -b $filename.te.Aligned.sortedByCoord.out.bam -bs 10 -p 4 --scaleFactor $var5 --normalizeUsing CPM -o $filename.te.cpm.norm.bw

wait

#calculate reads for condition (dm6) / total reads

var4=$(samtools view -c $filename.dm6.Aligned.sortedByCoord.out.bam)

var5=$(bc <<<"scale=10; $var4 / $var3")

echo $var5 > $filename.dm6.scaling.factor.txt

bamCoverage -b $filename.dm6.Aligned.sortedByCoord.out.bam -bs 10 -p 4 --scaleFactor $var5 --normalizeUsing CPM -o $filename.dm6.cpm.norm.bw

wait

exit