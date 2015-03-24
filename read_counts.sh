#!/bin/bash

echo "Dataset	nReads	nBases	meanLength	alignedReads	meanAlignLen	meanIden" > read.counts.blasr.tsv

for f in *_reads.fa; do
	n=${f%_reads.fa}
	nReads=$(cat $f | grep -c ">channel")
	nBases=$(cat $f | grep -v ">channel" | wc -c)
	meanLength=$(python -c "print(round($nBases/float($nReads),1))")
	alignedReads=$(cat $n\_reads.fa.blasr.true.read.stats|sed 's/ .*//g'|sort -u|wc -l)
	meanAlignLen=$(cat $n\_reads.fa.blasr.true.read.stats|tr ' ' '\t'|cut -f 6,7|awk '{sum+=$2; sum-=$1; n+=1}END{print sum/n}')
	meanIden=$(cat $n\_reads.fa.blasr.true.read.stats|tr ' ' '\t'|awk '{sum+=$4; n+=1}END{print sum/n}')
	echo "$n	$nReads	$nBases	$meanLength	$alignedReads	$meanAlignLen	$meanIden"
done >> read.counts.blasr.tsv

echo "Dataset	nReads	nBases	meanLength	alignedReads	meanAlignLen	meanIden" > read.counts.last.tsv

for f in *_reads.fa; do
	n=${f%_reads.fa}
	nReads=$(cat $f | grep -c ">channel")
	nBases=$(cat $f | grep -v ">channel" | wc -c)
	meanLength=$(python -c "print(round($nBases/float($nReads),1))")
	alignedReads=$(cat $n\_reads.fa.last.true.read.stats|sed 's/ .*//g'|sort -u|wc -l)
	meanAlignLen=$(cat $n\_reads.fa.last.true.read.stats|cut -f 12|awk '{sum+=$1; n+=1}END{print sum/n}')
	meanIden=$(cat $n\_reads.fa.last.true.read.stats|cut -f 13 | awk '{sum+=$1; n+=1}END{print sum/n}')
	meanIden=$(python -c "print(round(100*$meanIden/$meanAlignLen,4))")
	echo "$n	$nReads	$nBases	$meanLength	$alignedReads	$meanAlignLen	$meanIden"
done >> read.counts.last.tsv

