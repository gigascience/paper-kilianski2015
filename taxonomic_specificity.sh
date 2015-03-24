#!/bin/sh

#Count the number of bases aligned to the correct genus, species, and strain, as well as the incorrect genus

count_bases(){
	fp=$1
	term=$2
	type=$3

	[[ "$type" -eq "blasr" ]] && readcol=1
	[[ "$type" -eq "last" ]] && readcol=7

	for read in $(egrep $term $fp | tr ' ' '\t' | cut -f $readcol | sort -u); do
		echo $read
		#[[ "$type" -eq "blasr" ]] && egrep $term $fp | grep $read | tr ' ' '\t' | cut -f 6,7 | awk '{print $2-$1}'|sort -u | sort -nrk1 | head -1
		#[[ "$type" -eq "last" ]] && egrep $term $fp | grep $read | tr ' ' '\t' | cut -f 12 | sort -nrk1 | head -1
	done #| tr -d '-' | awk '{sum+=$1}END{print sum}'
}

count_bases_absent(){
	fp=$1
	term=$2
	type=$3

	[[ "$type" -eq "blasr" ]] && readcol=1
	[[ "$type" -eq "last" ]] && readcol=7

	for read in $(egrep -v $term $fp | grep -v readname | tr ' ' '\t' | cut -f $readcol | sort -u); do 
		[[ "$type" -eq "blasr" ]] && egrep -v $term $fp | grep $read | tr ' ' '\t' | cut -f 6,7 | awk '{print $3-$2}' | sort -nrk1 | head -1
		[[ "$type" -eq "last" ]] && egrep -v $term $fp | grep $read | tr ' ' '\t' | cut -f 12 | sort -nrk1 | head -1
	done | tr -d '-' | awk '{sum+=$1}END{print sum}'
}

count_reads(){
	fp=$1
	term=$2
	type=$3

	[[ "$type" -eq "blasr" ]] && egrep $term $fp | tr ' ' '\t' | cut -f 1 | sort -u | wc -l
	[[ "$type" -eq "last" ]] && egrep $term $fp | tr ' ' '\t' | cut -f 7 | sort -u | wc -l 
}
count_reads_absent(){
	fp=$1
	term=$2
	type=$3

	[[ "$type" -eq "blasr" ]] && readcol=1
	[[ "$type" -eq "last" ]] && readcol=7

	egrep -v $term $fp | fgrep -v readname | tr ' ' '\t' | cut -f $readcol | sort -u | wc -l 
}
longest_alignment(){
	fp=$1
	term=$2
	type=$3

	egrep $term $fp |tr ' ' '\n' | tr '\t' '\n' | grep channel | sort -u | while read readname; do
		[[ $type == "blasr" ]] && egrep $term $fp | grep $readname | tr ' ' '\t' | cut -f 6,7 | awk '{print $2-$1}' | sort -nrk1 | head -1 | tr -d '-'
		[[ $type == "last" ]] && egrep $term $fp | grep $readname | cut -f 12 | sort -nrk1 | head -1
	done
}
longest_alignment_absent(){
	fp=$1
	term=$2
	type=$3

	egrep -v $term $fp |tr ' ' '\n' | tr '\t' '\n' | grep channel | sort -u | while read readname; do
		[[ $type == "blasr" ]] && egrep -v $term $fp | grep $readname | tr ' ' '\t' | cut -f 6,7 | awk '{print $2-$1}' | sort -nrk1 | head -1 | tr -d '-'
		[[ $type == "last" ]] && egrep -v $term $fp | grep $readname | cut -f 12 | sort -nrk1 | head -1
	done
}
read_num_and_align_len(){
	echo "$#/$(for var in "$@"; do echo $var; done | awk '{sum+=$1}END{print sum}')"
}
taxon(){
	fp=$1
	type=$4
	genus_dat=$(longest_alignment $fp $2 $type)
	genus=$(read_num_and_align_len $genus_dat)
	species_dat=$(longest_alignment $fp $3 $type)
	species=$(read_num_and_align_len $species_dat)
	nongenus_dat=$(longest_alignment_absent $fp $2 $type)
	nongenus=$(read_num_and_align_len $nongenus_dat)
	echo "$fp	$genus	$species	$nongenus"
}
for p in last blasr; do 
	echo "Dataset	GenusReads/GenusBases	SpeciesReads/SpeciesBases NonGenusReads/NonGenusBases" > taxonomic_specificity.$p.tsv
	taxon ecoli_reads.fa.$p.read.stats 'Escherichia|Shigella' coli $p >> taxonomic_specificity.$p.tsv
	taxon MVA_reads.fa.$p.read.stats 'Vaccinia|Cowpox' 'MVA|Cowpox' $p >> taxonomic_specificity.$p.tsv
	taxon lister_reads.fa.$p.read.stats Vaccinia Lister $p >> taxonomic_specificity.$p.tsv
	taxon cowpox_reads.fa.$p.read.stats Cowpox Cowpox $p >> taxonomic_specificity.$p.tsv
done
