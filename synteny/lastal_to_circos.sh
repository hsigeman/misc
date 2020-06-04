#!/bin/bash

# This script takes a lastal output file and plots circos plots for each chromosome separately
# Filtering steps: 
## 1. At least 100kb and 1 % of scaffold must align to each chromosome to assign the scaffold to it. 
## 2. Then 1 Mb (doesn't make sense, will change it)


# Modify these
lastal="intermediate/lastal_ZF/AlaArv_ref_AlaArv_EDI/AlaArv_align_converted"
genome_fai="data/external_raw/genome/GCA_902810485.1_skylark_genome_genomic.fasta.fai"
query="AlaArv"
singleton="intermediate/freebayes_17nov2019_parallel/AlaArv_ref_AlaArv_EDI/AlaArv.singletons.bed"


lastal="intermediate/lastal_ZF/CisJun_ref_CisJun_B10K/CisJun_align_converted"
genome_fai="data/external_raw/genome/B10K-DU-002-30.genomic.fasta.fai"
query="CisJun"
singleton="intermediate/freebayes_17nov2019_parallel/CisJun_ref_CisJun_B10K/CisJun.singletons.bed"


lastal="intermediate/lastal_ZF/Sylvietta_virens_B10K_align_converted"
genome_fai="data/external_raw/genome/Sylvietta_virens_B10K.genomic.fasta.fai"
query="SylVir"
singleton="intermediate/freebayes_17nov2019_parallel/SylBra_ref_SylBra_B10K/SylBra.singletons.bed"


# Keep these the same
target="zf"
outdir="intermediate/lastal_ZF/circos_plotting"
#circos_temp="../neosexchromosome/intermediate/satsuma_warbler_gt_with_mt/circos_100kb_1_perc/conf_circos_template.50kb.labelFix.conf"
#circos_temp="../neosexchromosome/intermediate/satsuma_warbler_gt_with_mt/circos_100kb_1_perc/conf_circos_template.50kb.conf2"
mkdir $outdir

# CIRCOS TEMPLATE FILE
circos_temp="conf_circos_template.50kb.labelFix.conf"
#cp code/${circos_temp} ${outdir}

source activate neosc

# Create comparative species karyotype file
cat data/external_raw/genome/Taeniopygia_guttata.taeGut3.2.4.dna.toplevel_final.fasta.fai | awk '{ print "chr","-","zf"$1,$1,"0",$2,"red"}' > ${outdir}/${target}_karyotype.txt

#cp ./intermediate/satsuma_warbler_zf_circos/conf_circos_template.50kb.conf ${outdir}/
#cp $circos_temp ${outdir}

####################################
# Optional: count the number of singletons per 5kb and sex
bedtools makewindows -g $genome_fai -w 100000 > ${outdir}/${query}.genome.windows.bed

cat data/meta/samples_sex.tsv | while read sample sp sex ; do cat $singleton | grep -v CHROM | grep $sample | bedtools coverage -a ${outdir}/${query}.genome.windows.bed -b stdin -counts | awk '{print "'"$query"'"$0}' > ${outdir}/${query}_${sex}.singletons.bed ; done

join -t ':' -1 1 -2 1  <(awk  '{print $1$2$3 ":" $0;}' ${outdir}/${query}_female.singletons.bed | sort -t ':' -k1,1 ) <(awk '{print $1$2$3 ":" $0;}' ${outdir}/${query}_male.singletons.bed | sort -t ':' -k1,1 ) | cut -d ':' -f 2- | sed 's/:/\t/' | sed 's/ /\t/g' | cut -f 1,2,3,4,8 | awk '{print $1,$2,$3,$4-$5}' > ${outdir}/${query}_singletonDiff.out



# Calculate prop of alignments to the ZF chromosomes (first one gives the each number of matching base pairs)
cat $lastal | awk '{print $1,$10"_CHR_"$14}' | awk '{a[$2]+=$1} END {OFS="\t"; for(i in a) print a[i], i}' | awk '$1>10000 {print}' | sed 's/_CHR_/\t/' | while read bp scaff chr ; do cat $genome_fai | awk '$1=="'"$scaff"'" {print $2}' | awk -v OFS="\t" '{print "'"$scaff"'","'"$chr"'","'"$bp"'",$1,"'"$bp"'"/$1}' ; done > ${outdir}/${query}_${target}_scaffold_chr_prop_assignment.out


# Only print those where i) matches to zebra finch is larger than 1% and scaffold match is > 100kb. 
cat ${outdir}/${query}_${target}_scaffold_chr_prop_assignment.out | awk '$3 > 100000 && $5 > 0.01 {print $0}' | sort -k2,2g > ${outdir}/${query}_${target}_scaffold_chr_prop_assignment_100kb.out 


# Write to a circos links file
awk 'NR==FNR{c[$1$2]++;next};c[$10$14] > 0' ${outdir}/${query}_${target}_scaffold_chr_prop_assignment_100kb.out $lastal | awk '{print "'"$query"'"$10,$12,$13,"'"$target"'"$14,$16,$17}' > ${outdir}/${query}_${target}_circos_links_allChr.txt


# Filter for scaffolds that align with more than 1 Mb to this chromosome
cat ${outdir}/${query}_${target}_circos_links_allChr.txt | awk '{print $3-$2,$1"_CHR_"$4}' | awk '{a[$2]+=$1} END {OFS="\t"; for(i in a) print a[i], i}' | awk '$1>500000 {print $2}' | sed 's/_CHR_/\t/' > ${outdir}/${query}_${target}_link_filter.list

awk 'NR==FNR{c[$1$2]++;next};c[$1$4] > 0' ${outdir}/${query}_${target}_link_filter.list ${outdir}/${query}_${target}_circos_links_allChr.txt > ${outdir}/${query}_${target}_circos_links_allChr_filter.txt


# And for each chromosome separately
cd ${outdir}
rm ${query}_${target}_*_links.txt
awk '{print >> "'"$query"'""_""'"$target"'""_"$4"_links.txt"}' ${query}_${target}_circos_links_allChr_filter.txt
cd ../../../


# Make circos karyotype file
cat $genome_fai | cut -f 1,2 | sed 's/Contig//' | awk '{print "chr","-","'"$query"'"$1,$1,"0",$2,"blue"}' > ${outdir}/${query}_karyotype.txt


### Parallel bundlelinks - 50kbp and 5kbp
cd ${outdir}

find . -name "*_links.txt" | grep $query | grep -v bundles |  sed 's|./||'> ${query}_links.list

cat ${query}_links.list | sed 's/.txt//' | while read links ; do ~/bin/circos-tools-0.23/tools/bundlelinks/bin/bundlelinks -max_gap 50000 -min_bundle_size 50000 -strict -links $links.txt > $links.bundles.txt ; done



#############################################################
# 100kb - For automatic scaffold order and colours

ls | grep bundles.txt | grep $query | sed 's/_links.bundles.txt//' | sed "s/${query}_${target}_//" | while read chr ; do cat ${query}_${target}_${chr}_links.bundles.txt | cut -f 1 -d " " | sort | uniq | tac - | tr '\n' ';' | sed 's/;$//' | awk '{print "'"$chr"'"";"$0}' | while read line ; do cat circos_template.conf | sed -e "s/CHR/$chr/" | sed -e "s/chromosomes = LINE/chromosomes = $line/" ; done | sed "s/TARGET/${target}/" | sed "s/QUERY/${query}/"  > ${query}_${target}_${chr}.conf ; done


# Bash line for finding the middle of the largest chunk of alignment and order the grw scaffolds in the reverse order
ls | grep bundles.txt | grep $query | sed 's/_links.bundles.txt//' | sed "s/${query}_${target}_//" | while read chr ; do cat ${query}_${target}_${chr}_links.bundles.txt | sed 's/ /\t/g' | cut -f 1 | sort | uniq | while read scaff ; do cat ${query}_${target}_${chr}_links.bundles.txt | sed 's/ /\t/g' |  cut -f 1-3 | awk '$1=="'"$scaff"'" { print $1,$3-$2,$2 }' | sed 's/ /\t/g' | sort -k2,2gr | head -n 1 ; done | while read scaff length start ; do cat ${query}_${target}_${chr}_links.bundles.txt | sed 's/ /\t/g' | awk '{ if($1=="'"$scaff"'" && $2=="'"$start"'" && $6>$5) printf "%s %.0f %.0f %.0f\n", $1,$5,$6,$5+(($6-$5)/2) ; else if ($1=="'"$scaff"'" && $2=="'"$start"'" && $5>$6) printf "%s %.0f %.0f %.0f\n", $1,$5,$6,$6+(($5-$6)/2) }'; done | sed 's/ /\t/g' | sort -k4,4gr | cut -f 1 | tr -d "\n" | sed "s/${query}/;${query}/g" | awk '{print "'"$chr"'"$0}' | while read line ; do cat ${query}_${target}_${chr}.conf | sed -e "s/chromosomes_order = LINE/chromosomes_order = $line/" ; done > ${query}_${target}_${chr}.conf2 ; done


# Assign different colors to the different scaffolds - This one also takes the scaffold order into account
ls | grep bundles.txt | grep $query | sed 's/_links.bundles.txt//' | sed "s/${query}_${target}_//" | while read chr; do cat ${query}_${target}_${chr}_links.bundles.txt | sed 's/ /\t/g' | cut -f 1 | sort | uniq | while read scaff ; do cat ${query}_${target}_${chr}_links.bundles.txt | sed 's/ /\t/g' |  cut -f 1-3 | awk '$1=="'"$scaff"'" { print $1,$3-$2,$2 }' | sed 's/ /\t/g' | sort -k2,2gr | head -n 1 ; done | while read scaff length start ; do cat ${query}_${target}_${chr}_links.bundles.txt | sed 's/ /\t/g' | awk '{ if($1=="'"$scaff"'" && $2=="'"$start"'" && $6>$5) printf "%s %.0f %.0f %.0f\n", $1,$5,$6,$5+(($6-$5)/2) ; else if ($1=="'"$scaff"'" && $2=="'"$start"'" && $5>$6) printf "%s %.0f %.0f %.0f\n", $1,$5,$6,$6+(($5-$6)/2) }'; done | sed 's/ /\t/g' | sort -k4,4gr | cut -f 1 | cat -n | sed 's/ /\t/g' | while read nr scaff ; do cat ${query}_${target}_${chr}_links.bundles.txt | sed 's/ /\t/g' | awk '$1=="'"$scaff"'" { print $0"color=""'"$nr"'" }' ; done > ${query}_${target}_${chr}_links.bundles.color.txt ; done

ls | grep bundles.color.txt | grep $query | sed 's/_links.bundles.color.txt//' | sed "s/${query}_${target}_//" | while read chr ; do cat ${query}_${target}_${chr}_links.bundles.color.txt | cut -f 1 | sort | uniq | while read scaff ; do cat ${query}_karyotype.txt | awk '$3=="'"$scaff"'" {print}' ; done > ${query}_karyotype_${chr}.txt ; done



mkdir plots
cat ${target}_karyotype.txt | cut -d " " -f 3 | while read chr; do circos -conf ${query}_${target}_${chr}.conf2 ; mv circos.png plots/${query}_${target}_${chr}.png ; mv circos.svg plots/${query}_${target}_${chr}.svg ; rm circos.png ; rm circos.svg ; done



