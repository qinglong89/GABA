#!/bin/bash

WorkDir=/mnt/home1/qinglong/@Shotgun/@CDI_analysis/co-assembly_BR71-BR79/anvio7_process


mkdir $WorkDir/BinTaxa_VoteFrom_gene-level-taxonomy
cd $WorkDir/BinTaxa_VoteFrom_gene-level-taxonomy

for bin in $WorkDir/SAMPLES-SUMMARY/bin_by_bin/*
do
	BinName=`basename $bin`

	#extract gene caller IDs of a specific metagenomic bin
	cat $bin/''$BinName'-gene_calls.txt' | cut -d$'\t' -f 1 | sed 1d > bin_genes_list.txt

	#extract the taxa of all gene calls of a specific metagenomic bin
	grep -w -F -f bin_genes_list.txt $WorkDir/gene_calls_DIAMOND-MEGAN_taxonomy_parsed4anvio.txt \
			> gene_calls_DIAMOND-MEGAN_taxonomy_parsed4anvio_bin.txt	#exclude header

	TotalGenes=`cat gene_calls_DIAMOND-MEGAN_taxonomy_parsed4anvio_bin.txt | wc -l`

	#parse for each taxonomic level
	for column in 2 3 4 5 6 7 8
	do
		#extract the best bacterial hit
		cat gene_calls_DIAMOND-MEGAN_taxonomy_parsed4anvio_bin.txt | cut -d$'\t' -f $column | sort | uniq -c | sort -nr | sed '/Unknown/d' | head -n 1 > temp0

		Count=`cat temp0 | sed 's/d__/\td__/g' | sed 's/p__/\tp__/g' | sed 's/c__/\tc__/g' | sed 's/o__/\to__/g' \
				 | sed 's/f__/\tf__/g' | sed 's/g__/\tg__/g' | sed 's/s__/\ts__/g' | cut -d$'\t' -f 1`

		Confidence=`echo "scale=2; $Count*100/$TotalGenes" | bc`

		Taxa=`cat temp0  | sed 's/d__/\td__/g' | sed 's/p__/\tp__/g' | sed 's/c__/\tc__/g' | sed 's/o__/\to__/g' \
				 | sed 's/f__/\tf__/g' | sed 's/g__/\tg__/g' | sed 's/s__/\ts__/g' |cut -d$'\t' -f 2`

		rm temp0

		echo "$Taxa;$Confidence"

	done >> temp1

	taxonomy=`cat temp1 | tr $'\n' $'\t'`	#transpose: from lines to rows

	echo -e "$BinName\t$taxonomy" > ''$BinName'_taxonomy_from_genes.txt'

	rm temp1 gene_calls_DIAMOND-MEGAN_taxonomy_parsed4anvio_bin.txt bin_genes_list.txt
done

cat *_taxonomy_from_genes.txt > AllBins_taxonomy_from_gene-level-taxonomy.txt

rm *_taxonomy_from_genes.txt
