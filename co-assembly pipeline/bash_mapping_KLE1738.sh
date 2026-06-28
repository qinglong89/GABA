#!/bin/bash

WorkDir=/mnt/home1/qinglong/@Shotgun/@CDI_analysis/co-assembly_BR71-BR79
THREADS=80

for folder in $WorkDir/RawReadProcess/*
do

	cd $folder

	SampleID=`basename $folder`

	read1=`ls *R1*`
	read2=`ls *R2*`

	#note: bowtie2 database has to be built with the same version of bowtie2 used for mapping
	bowtie2 -x /mnt/home1/qinglong/@Shotgun/@Firmicutes_CAG:114_detection/KLE1738/KLE1738 \
		-1 $read1 \
		-2 $read2 \
		--al-conc-gz ''$SampleID'_bowtie2_KLE1738_mapped' \
		-S ''$SampleID'_bowtie2_KLE1738.sam' \
		-p $THREADS --very-sensitive

	#convert to bam format and exclude the unmapped reads (-F 4)
	#samtools view -F 4 -@ $THREADS -b -S ''$SampleID'_bowtie2_KLE1738.sam' > ''$SampleID'_bowtie2_KLE1738_mapped.bam'	#for single-end sequencing
	samtools view -F 4 -f 0x2 -@ $THREADS -b -S ''$SampleID'_bowtie2_KLE1738.sam' > ''$SampleID'_bowtie2_KLE1738_mapped.bam'       #include only properly paired reads for downstream

	#uniquely mapped reads with MAPQ not less than 2: if average per-read quality score is higher than Q20, then there will be less than 5 mismatches per alignment
	samtools view -q 2 -@ $THREADS -b -S ''$SampleID'_bowtie2_KLE1738_mapped.bam' > ''$SampleID'_bowtie2_KLE1738_MAPQ>=2.bam'

	#sort (by position in reference genome) and index bam file (output two files – .bam and .bam.bai)
	samtools sort -@ $THREADS ''$SampleID'_bowtie2_KLE1738_MAPQ>=2.bam' -o ''$SampleID'_bowtie2_KLE1738_MAPQ>=2_sorted.bam'
	samtools index -b -@ $THREADS ''$SampleID'_bowtie2_KLE1738_MAPQ>=2_sorted.bam'

	#visualize alignments
	#samtools tview ''$SampleID'_bowtie2_KLE1738_MAPQ>=2_sorted.bam' /mnt/home1/qinglong/@Shotgun/@Firmicutes_CAG:114_detection/KLE1738/Evtepia_gabavorous_KLE1738.fna

	#generate stats from mapping
	echo ''$SampleID':' 'total pairs of paired reads:' $(zcat $read1 | wc -l) \
	                    'MAPQ>0:' $(samtools view -@ $THREADS -c -q 0 ''$SampleID'_bowtie2_KLE1738_mapped.bam') \
	                    'MAPQ>2:' $(samtools view -@ $THREADS -c -q 2 ''$SampleID'_bowtie2_KLE1738_mapped.bam') \
	                    'MAPQ>5:' $(samtools view -@ $THREADS -c -q 5 ''$SampleID'_bowtie2_KLE1738_mapped.bam') \
	                    'MAPQ>40:' $(samtools view -@ $THREADS -c -q 40 ''$SampleID'_bowtie2_KLE1738_mapped.bam') \
	                    'MAPQ>42:' $(samtools view -@ $THREADS -c -q 42 ''$SampleID'_bowtie2_KLE1738_mapped.bam') \
			    'Coverage_mean:' $(samtools depth -a ''$SampleID'_bowtie2_KLE1738_MAPQ>=2_sorted.bam' | awk '{sum+=$3} END {print sum/NR}') \
			    'Coverage_stdev:' $(samtools depth -a ''$SampleID'_bowtie2_KLE1738_MAPQ>=2_sorted.bam' | awk '{sum+=$3; a[NR]=$3} END {for(i in a) y+=(a[i]-(sum/NR))^2; print sqrt(y)/(NR-1)}') \
			    'Coverage_breadth (%):' $(samtools depth -a ''$SampleID'_bowtie2_KLE1738_MAPQ>=2_sorted.bam' | awk '{RefCount++; if($3>0) total+=1} END {print (total/RefCount)*100}') \
			    > ''$SampleID'_mapping_stats.txt'
	#The -a flag (samtools depth) indicates that the depth must be calculated at all positions, including at those with zero depth

	rm ''$SampleID'_bowtie2_KLE1738_mapped.1' ''$SampleID'_bowtie2_KLE1738_mapped.2' ''$SampleID'_bowtie2_KLE1738.sam'
	rm ''$SampleID'_bowtie2_KLE1738_mapped.bam' ''$SampleID'_bowtie2_KLE1738_MAPQ>=2.bam'
done
