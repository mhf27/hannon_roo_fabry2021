#!/bin/bash

# Analysis for Clust1, alignment to TE list only

source ~/software/MyVE/bin/activate


for i in *.fq; do (
        filename=${i%%.fq}
        echo 'Start' $filename;
        fastx_clipper -a AGATCGGAAGAGCACACGTCT -l 15 -c -Q33 -i $i -o $filename.clipped;
        seqtk trimfq -b 4 -e 4 $filename.clipped > $filename.trimmed;
        ) &
done

wait

echo 'Start Part 1'

for i in *.trimmed; do (
        filename=${i%%.trimmed}
        #2 mismatches per read, add --alignIntronMax 1 to suppress spliced reads, add --alignEndsType EndToEnd to account for difficult to map indels (makes it more relatable to bowtie)
        STAR --outSAMunmapped Within --genomeDir /Users/fabry01/genomes/drosophila/dm6/UCSC/dm6/Sequence/STAR-Index/te_fused_clean_3 --limitBAMsortRAM 1098266374 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1000000000 --outMultimapperOrder Random --outSAMmultNmax 1 --alignIntronMax 100 --alignEndsType EndToEnd --outFilterMismatchNoverLmax 0.04 --readFilesIn $filename.trimmed --runThreadN 2 --outFileNamePrefix $filename.;
        echo 'STAR dm6 done' $filename;
        samtools index $filename.Aligned.sortedByCoord.out.bam;
        samtools view -hf 4 $filename.Aligned.sortedByCoord.out.bam > $filename.unmap.sam;
        samtools view -S -b $filename.unmap.sam > $filename.unmap.bam;
        bamToFastq -i $filename.unmap.bam -fq $filename.unmap.fq;
        STAR --genomeDir /Users/fabry01/genomes/drosophila/dm6/UCSC/dm6/Sequence/STAR-Index/ --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 100000 --outMultimapperOrder Random --outSAMmultNmax 1 --alignIntronMax 100 --alignEndsType EndToEnd --outFilterMismatchNoverLmax 0.04 --readFilesIn $filename.unmap.fq --runThreadN 2 --outFileNamePrefix $filename.dm6.multi.;
        samtools index $filename.dm6.multi.Aligned.sortedByCoord.out.bam;
        echo 'Index bai done' $filename;
        ) &
done

wait



for i in *.Aligned.sortedByCoord.out.bam; do (
        filename=${i%%.Aligned.sortedByCoord.out.bam}
        samtools view -h $i |awk '(length($10) > 20) && (length($10) < 22) || $1 ~ /^@/' | samtools view -bS - > $filename.21bp.Aligned.sortedByCoord.out.bam
        samtools view -h $i |awk '(length($10) > 22) && (length($10) < 30) || $1 ~ /^@/' | samtools view -bS - > $filename.23_29bp.Aligned.sortedByCoord.out.bam
        ) &
done

wait

# the flag in samtools view removes all reads not associated with a chr or TE, that mean all read counts printed can be used immediately

for i in *.Aligned.sortedByCoord.out.bam; do (
        filename=${i%%.Aligned.sortedByCoord.out.bam}
        samtools view -f 0x10 -b $i > $filename.0.bam;
        samtools index $filename.0.bam;
        samtools idxstats $filename.0.bam | cut -f 1,3 > $filename.0.te.chrom_reads.txt;
        echo $filename.0.te.chrom_reads.txt;
        cut -f 2 $filename.0.te.chrom_reads.txt | paste -sd+ - | bc;
        ) &
done

wait

for i in *.Aligned.sortedByCoord.out.bam; do (
        filename=${i%%.Aligned.sortedByCoord.out.bam}
        samtools view -F 0x10 -b $i > $filename.16.bam;
        samtools index $filename.16.bam;
        samtools idxstats $filename.16.bam | cut -f 1,3 > $filename.16.te.chrom_reads.txt;
        echo $filename.16.te.chrom_reads.txt;
        cut -f 2 $filename.16.te.chrom_reads.txt | paste -sd+ - | bc;
        ) &
done

wait

exit
