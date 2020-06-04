#!/bin/bash

# This script takes a lastal output file and makes it into a compatible file for circos

# CIRCOS TEMPLATE FILE
circos_temp="../neosexchromosome/intermediate/satsuma_warbler_gt_with_mt/circos_100kb_1_perc/conf_circos_template.50kb.labelFix.conf"

lastal="intermediate/lastal_ZF/AlaArv_ref_AlaArv_EDI/AlaArv_align_converted"
genome_fai="data/external_raw/genome/GCA_902810485.1_skylark_genome_genomic.fasta.fai"
outdir="intermediate/lastal_ZF/circos_plotting"
target="zf"
query="AlaArv"
mkdir $outdir
circos_temp="conf_circos_template.50kb.labelFix.conf"
cp code/${circos_temp} ${outdir}

# Create comparative species karyotype file
cat data/external_raw/genome/Taeniopygia_guttata.taeGut3.2.4.dna.toplevel_final.fasta.fai | awk '{ print "chr","-","zf"$1,$1,"0",$2,"red"}' > ${outdir}/${target}_karyotype.txt

#cp ./intermediate/satsuma_warbler_zf_circos/conf_circos_template.50kb.conf ${outdir}/
#cp $circos_temp ${outdir}

####################################


# Calculate prop of alignments to the ZF chromosomes (first one gives the each number of matching base pairs)
cat $lastal | awk '{print $1,$10"_CHR_"$14}' | awk '{a[$2]+=$1} END {OFS="\t"; for(i in a) print a[i], i}' | awk '$1>10000 {print}' | sed 's/_CHR_/\t/' | while read bp scaff chr ; do cat $genome_fai | awk '$1=="'"$scaff"'" {print $2}' | awk -v OFS="\t" '{print "'"$scaff"'","'"$chr"'","'"$bp"'",$1,"'"$bp"'"/$1}' ; done > ${outdir}/${query}_${target}_scaffold_chr_prop_assignment.out


# Only print those where i) matches to zebra finch is larger than 1% and scaffold match is > 100kb. 
cat ${outdir}/${query}_${target}_scaffold_chr_prop_assignment.out | awk '$3 > 100000 && $5 > 0.01 {print $0}' | sort -k2,2g > ${outdir}/${query}_${target}_scaffold_chr_prop_assignment_100kb.out


# Write to a circos links file
awk 'NR==FNR{c[$1$2]++;next};c[$10$14] > 0' ${outdir}/${query}_${target}_scaffold_chr_prop_assignment_100kb.out $lastal | awk '$1>200 {print "'"$query"'"$10,$12,$13,"'"$target"'"$14,$16,$17,$1}' > ${outdir}/${query}_${target}_circos_links.txt


# And for each chromosome separately
cd ${outdir}
awk '{print >> "'"$query"'""_""'"$target"'""_"$4"_links.txt"}' ${query}_${target}_circos_links.txt
cd ../../../


# Make circos karyotype file
cat $genome_fai | cut -f 1,2 | sed 's/Contig//' | awk '{print "chr","-","grw"$1,$1,"0",$2,"blue"}' > ${outdir}/${query}_karyotype.txt
