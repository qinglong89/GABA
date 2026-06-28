#!/bin/bash

WorkDir=/mnt/home1/qinglong/@Shotgun/@CDI_analysis/co-assembly_BR71-BR79/anvio7_process


cd $WorkDir

for bin in $WorkDir/SAMPLES-SUMMARY/bin_by_bin/*
do
	BinName=`basename $bin`

	cd $bin

	#extract gene caller IDs of a specific metagenomic bin
	cat ''$BinName'-gene_calls.txt' | cut -d$'\t' -f 1 | sed 1d > bin_genes_list.txt

	cat bin_genes_list.txt | sed "s/$/\t$BinName/" > bin_genes_list_info.txt
done

cd $WorkDir

cat ./SAMPLES-SUMMARY/bin_by_bin/*/bin_genes_list_info.txt > list_gene_bin_info.txt

rm ./SAMPLES-SUMMARY/bin_by_bin/*/bin_genes_list*.txt
