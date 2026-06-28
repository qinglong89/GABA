#!/bin/bash
WorkDir=/mnt/home1/qinglong/@Shotgun/@CDI_analysis/co-assembly_BR71-BR79

THREADS=80

#########################################################################################################################################
##################################"Step 1: sequence normalization with bbnorm"###########################################################
#########################################################################################################################################
for folder in $WorkDir/RawReadProcess/*
do
	cd $folder

	#generate prefix
	SampleName=`basename $folder`

	#host read removal
	/home/qinglong/softwares/bowtie2-binary-2.3.4.3/bowtie2 -p $THREADS -x /home/qinglong/MetaDatabases/hg19_bowtie2_index/hg19 \
	                                                        -1 ''$SampleName'_R1.fastq.gz' \
        	                                                -2 ''$SampleName'_R2.fastq.gz' \
                	                                        --un-conc ''$SampleName'_nohost'

	#normalization with bbnorm.sh
#        /home/qinglong/softwares/bbmap-38.34/bbnorm.sh in1=''$SampleName'_nohost.1' in2=''$SampleName'_nohost.2' \
#      	                                               out1=''$SampleName'_nohost_bbnorm_R1.fq' out2=''$SampleName'_nohost_bbnorm_R2.fq' \
#              	                                       target=40 min=5 threads=$THREADS
done

#########################################################################################################################################
##################################"Step 2: co-assembly with megahit"#####################################################################
#########################################################################################################################################
mkdir -p $WorkDir/co-assembly && cd $WorkDir/co-assembly

cat $WorkDir/RawReadProcess/*/*_nohost.1 > AllSamples_nohost_R1.fq
cat $WorkDir/RawReadProcess/*/*_nohost.2 > AllSamples_nohost_R2.fq

/home/qinglong/softwares/megahit-binary-1.1.4/megahit -t $THREADS -m 0.9 --min-contig-len 2000 \
                                                      -1 AllSamples_nohost_R1.fq \
                                                      -2 AllSamples_nohost_R2.fq \
                                                      -o megahit_assemblies

rm AllSamples_nohost_R1.fq AllSamples_nohost_R2.fq

#########################################################################################################################################
##################################"Step 3: read mapping for Anvi'o processes"############################################################
#########################################################################################################################################
mkdir -p $WorkDir/anvio7_process && cd $WorkDir/anvio7_process

#activate virtual environment of anvio (v7) from minibioconda
source /home/DataAnalysis/miniconda2/bin/activate /home/qinglong/Programs/anvio7

#reformat the header identifier of the non-redundant contigs and only keep contigs with a minimum length of 2000 nucleotides
anvi-script-reformat-fasta $WorkDir/co-assembly/megahit_assemblies/final.contigs.fa -o $PWD/contigs.fa -l 2000 --simplify-names

#build the index for the contigs
mkdir contigs-index
bowtie2-build --threads $THREADS --quiet contigs.fa ./contigs-index/contigs   #the prefix for the index is "contigs"

#generate BAM files for each sample with the contigs index
for folder in $WorkDir/RawReadProcess/*
do
	cd $folder

	SampleName=`basename $folder`

	bowtie2 --threads $THREADS --quiet --very-sensitive -x $WorkDir/anvio7_process/contigs-index/contigs \
		-1 ''$SampleName'_nohost.1' -2 ''$SampleName'_nohost.2' \
		--un-conc ''$SampleName'_coassembly_unmatched' \
		-S ''$SampleName'.sam'   # QC reads from step 1

	samtools view -@ $THREADS -F 4 -bS ''$SampleName'.sam' > ''$SampleName'-RAW.bam'	#exclude unmapped reads

	anvi-init-bam --num-threads $THREADS ''$SampleName'-RAW.bam' -o ''$SampleName'.bam'	#sort and index for anvio downstream

	#generate stats from mapping
	echo ''$SampleName':' 'q>0:' $(samtools view -@ $THREADS -c -q 0 ''$SampleName'.bam') \
                              'q>2:' $(samtools view -@ $THREADS -c -q 2 ''$SampleName'.bam') \
                              'q>5:' $(samtools view -@ $THREADS -c -q 5 ''$SampleName'.bam') \
                              'q>40:' $(samtools view -@ $THREADS -c -q 40 ''$SampleName'.bam') \
                              'q>42:' $(samtools view -@ $THREADS -c -q 42 ''$SampleName'.bam') \
                              > ''$SampleName'_mapping_stats.txt' &

	rm ''$SampleName'.sam' ''$SampleName'-RAW.bam'   #remove unnecessary large files to save space
done

#########################################################################################################################################
##################################"Step 4: anvi'o profiling & metagenomic binning"#######################################################
#########################################################################################################################################
#anvi'o processes for contigs database

cd $WorkDir/anvio7_process
anvi-gen-contigs-database -f contigs.fa -o contigs.db -n 'non-redundant contigs databse'
anvi-run-hmms -c contigs.db --num-threads $THREADS

#anvi-profile
for folder in $WorkDir/RawReadProcess/*
do
	cd $folder

	SampleName=`basename $folder`

	#anvio needs a name for each profile (sample), rules include ASCII letters, digits, and the underscore character ('_'); do not start with numbers
	NewName=`basename $folder | sed 's/-/_/g'`

	anvi-profile -i ''$SampleName'.bam' \
		     -c $WorkDir/anvio7_process/contigs.db \
		     --min-contig-length 2000 \
		     --num-threads $THREADS \
		     --sample-name $(echo "S_$NewName") \
		     --skip-SNV-profiling
done

#anvi-merge
anvi-merge $WorkDir/RawReadProcess/*/*ANVIO_PROFILE/PROFILE.db -o $WorkDir/anvio7_process/SAMPLES-MERGED \
	   -c $WorkDir/anvio7_process/contigs.db \
	   --skip-hierarchical-clustering

#metagenomic binning: metabat2 has high accuracy (low redundancy of bins) for large dataset and requires very less RAM (important)
anvi-cluster-contigs -p $WorkDir/anvio7_process/SAMPLES-MERGED/PROFILE.db \
		     -c $WorkDir/anvio7_process/contigs.db \
		     --collection-name metabat2 \
		     --driver metabat2 \
		     --num-threads $THREADS --just-do-it

#this command will include most contigs in the output bins: by default include >2500 bp contigs, --minCV and --minCVSum are two most important parameters to include more contigs for binning
anvi-cluster-contigs -p $WorkDir/anvio7_process/SAMPLES-MERGED/PROFILE.db \
		     -c $WorkDir/anvio7_process/contigs.db \
		     --collection-name metabat2_refined \
		     --driver metabat2 \
		     -m 1500 --minCV 0.01 --minCVSum 0.1 --minClsSize 100000 \
		     --num-threads $THREADS --just-do-it

#########################################################################################################################################
##################################"Step 5: gene-level taxonomy & gene function"##########################################################
#########################################################################################################################################
#gene-level taxonomy assignment, following Genome Biology, 2017, 18: 182.
cd $WorkDir/anvio7_process

anvi-get-sequences-for-gene-calls -c contigs.db -o gene_calls.fa

#make sure BLASTDB variable (nt nr & taxonomy) is indicated in $PATH
#requires 100 GB RAM to store nt database, time-consuming, and less genes can be annotated in this way
#/home/qinglong/softwares/ncbi-blast-binary-2.9.0/bin/blastn -task megablast \
#		-query gene_calls.fa \
#		-db /mnt/home1/qinglong/Databases/BLAST_NT_DB_June2021/nt \
#		-outfmt '6 std staxids scomnames sscinames sskingdoms' \
#		-out gene_calls.blastn \
#		-num_threads $THREADS -evalue 1e-20

#/home/qinglong/softwares/MEGAN6/tools/blast2lca --input gene_calls.blastn \
#						--topPercent 10 --format BlastTap \
#						--acc2taxa /mnt/home1/qinglong/Databases/MEGAN6_acc2tax_maps/nucl_acc2tax-Jul2019.abin \
#						--output gene_calls_BLASTN-MEGAN_taxonomy.txt


#DIAMOND alignment (output: gene_calls.daa), only requires 20 GB RAM and it is quick
/home/qinglong/Programs/diamond_v2.0.9/diamond blastx -p $THREADS \
		-d /mnt/home1/qinglong/Databases/DIAMOND_NR_DB_June2021/nr_June2021 \
		-q gene_calls.fa \
		-a gene_calls

/home/qinglong/Programs/diamond_v2.0.9/diamond view -a gene_calls.daa -o gene_calls.m8	#convert to blast tabular format

/home/qinglong/softwares/MEGAN6/tools/blast2lca --input gene_calls.daa \
						--topPercent 10 --format DAA \
						--acc2taxa /mnt/home1/qinglong/Databases/MEGAN6_acc2tax_maps/prot_acc2tax-Jul2019X1.abin \
						--output gene_calls_DIAMOND-MEGAN_taxonomy.txt


#"Here in this analysis, we found Firmcutes CAG:114 or Evtepia gabavorous from DIAMOND-MEGAN taxonomy, but not in BLASTN-MEGAN taxonomy"
#parse DIAMOND-MEGAN taxonomy output for updating anvi'o contig database
echo "gene_callers_id;t_domain;t_phylum;t_class;t_order;t_family;t_genus;t_species" > temp0
cat gene_calls_DIAMOND-MEGAN_taxonomy.txt | cut -d ';' -f 1,3,5,7,9,11,13,15 > temp1
cat temp0 temp1	| sed 's/;/\t/g' > temp2
#make sure every line has 8 fields, otherwise anvi-import-taxonomy-for-genes won't work
awk -v c=8 'BEGIN{FS=OFS="\t"} {for(i=NF+1; i<=c; i++) $i="Unknown"} 1' temp2 > temp3	#still remain one field empty since last field of each line has \t ending
cat temp3 | sed 's/\t\t/\tUnknown\t/g' > gene_calls_DIAMOND-MEGAN_taxonomy_parsed4anvio.txt	#this is clean now
rm temp*

#awk 'BEGIN{FS="\t"}{print NF}' gene_calls_DIAMOND-MEGAN_taxonomy_parsed4anvio.txt | sort -n | uniq -c | sort -nr	#check again if every lines has 8 fields

anvi-import-taxonomy-for-genes -c contigs.db -i gene_calls_DIAMOND-MEGAN_taxonomy_parsed4anvio.txt -p default_matrix


#anvi'o v7 only support eggNOG mapper v2.0.0, v2.0.1
#perform gene annotation with eggNOG mapper: https://github.com/eggnogdb/emapper-benchmark/blob/master/benchmark_analysis.ipynb
anvi-get-sequences-for-gene-calls -c contigs.db --get-aa-sequences -o gene_calls_aa.fa

source /home/DataAnalysis/miniconda2/bin/activate /home/qinglong/Programs/eggnog-mapper_v2.0.1

#download_eggnog_data.py --data_dir /mnt/home1/qinglong/Databases/EggNOG-mapper_v2.0.1_DB

emapper.py  -i gene_calls_aa.fa \
	    --data_dir /mnt/home1/qinglong/Databases/EggNOG-mapper_v2.0.1_DB \
	    -o gene_calls_aa \
	    --cpu $THREADS -m diamond

rm gene_calls_aa.emapper.seed_orthologs

#need to add prefix 'g' to the gene caller IDs, otherwise anvi'o will not take it
sed -n '4p' gene_calls_aa.emapper.annotations > temp0	#extract header line only
cat gene_calls_aa.emapper.annotations | sed '/^#/d' | sed 's/^/g/' > temp1	#add prefix to gene caller IDs
cat temp0 temp1 > gene_calls_aa.emapper.annotations_parsed4anvio
rm temp*

source /home/DataAnalysis/miniconda2/bin/deactivate /home/qinglong/Programs/eggnog-mapper_v2.0.1


source /home/DataAnalysis/miniconda2/bin/activate /home/qinglong/Programs/anvio7
#anvi-setup-ncbi-cogs --num-threads $THREADS	#EggNOGmapper class depends COG data, this one-time run
anvi-script-run-eggnog-mapper -c contigs.db \
			      --annotation gene_calls_aa.emapper.annotations_parsed4anvio \
                              --num-threads $THREADS --drop-previous-annotations \
			      --use-version 2.0.1	#import gene function

#########################################################################################################################################
##################################"Step 6: taxonomy estimation for single-copy genes (SCGs) & metagenomic bins"##########################
#########################################################################################################################################
#"perform taxonomic assignment for single-copy genes"

cd $WorkDir/anvio7_process
source /home/DataAnalysis/miniconda2/bin/activate /home/qinglong/Programs/anvio7

#Setting up anvi’o SCG taxonomy, "only need to run the set-up once"
#anvi-setup-scg-taxonomy && anvi-setup-scg-taxonomy --reset

#identify all the 22 single-copy core genes in the contig database
anvi-run-scg-taxonomy -c contigs.db --num-parallel-processes 22 --num-threads 2

#option "--update-profile-db-with-taxonomy" and "--compute-scg-coverages"
#can show "stacked taxonomic abundance plot" in the `layers` tab of your interactive interface.
anvi-estimate-scg-taxonomy -c contigs.db --metagenome-mode -p SAMPLES-MERGED/PROFILE.db \
                           --compute-scg-coverages \
                           --simplify-taxonomy-information \
                           --update-profile-db-with-taxonomy \
                           --num-threads $THREADS

#add sample groupings and the order of samples
#anvi-import-misc-data sample_groups.txt --target-data-table layers -p ./SAMPLES-MERGED/PROFILE.db

#########################################################################################################################################

#anvi-interactive -p ./SAMPLES-MERGED/PROFILE.db -c contigs.db -C metabat2 --server-only

#use "anvi-refine" to do human-guided refinement for MAG bins with high redundancy (> 10%; "gene-level taxonomy will be useful in refinement process"

#after human-guided bin refinement, perform bin-level taxonomic estimation
anvi-estimate-scg-taxonomy -c contigs.db -C metabat2 -p SAMPLES-MERGED/PROFILE.db \
                           --compute-scg-coverages \
                           --simplify-taxonomy-information \
                           --num-threads $THREADS

anvi-summarize -p SAMPLES-MERGED/PROFILE.db -c contigs.db -o SAMPLES-SUMMARY -C metabat2 --just-do-it

