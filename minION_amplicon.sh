#!/bin/sh

dataset="$1" #Name to describe dataset
fp=$dataset\_reads.fa #Path to FASTQ file that will be created
fast5_folder="$2" #Folder containing all fast5 files for dataset
ref="$3" #FASTA all possible hits
true_ref="$4" #FASTA true reference
fo=$fp

#Path to executables
blasr='bin/blasr'
samtools='bin/samtools'
pileup_parser='bin/pileup.parser.py'
parse_maf='bin/parse_maf.py'
#poRe must be installed via R -- Version 0.6 used.

if [ ! -s $ref ]; then echo "$ref is empty or does not exist"; exit; fi

if [ ! -s $true_ref ]; then echo "$true_ref is empty or does not exist"; exit; fi

extract_reads(){
	#Extract fastq files from fast5 directory
	echo "library(poRe)
	extract.run.fasta('$fast5_folder')"| R --no-save --no-restore --vanilla

	#Concatenate all of the reads into a single file
	cat $fast5_folder/extracted/*template.fasta > $fp
}

run_blasr(){
	#Run BLASR on the whole set and the true ref, outputting per alignment details
	$blasr -m 4 $fp $ref > $fo.blasr.read.stats 
	$blasr -m 4 $fp $true_ref > $fo.blasr.true.read.stats 
	$blasr -sam $fp $true_ref  | $samtools view -bT $ref - | $samtools sort - $fo.blasr.true.srt 
	$samtools index $fo.blasr.true.srt.bam
	$samtools mpileup -f $true_ref $fo.blasr.true.srt.bam | $pileup_parser > $fo.blasr.true.pileup
}

run_last(){
	[ ! -s ${ref%.fa}.prj} ] && lastdb ${ref%.fa} $ref
        lastal -s 2 -T 0 -Q 0 -r 1 -a 1 -b 1 ${ref%.fa} $fp > $fo.maf 
        #Per read stats
        cat $fo.maf | ./$parse_maf > $fo.last.read.stats

	[ ! -s ${true_ref%.fa}.prj} ] && lastdb ${true_ref%.fa} $true_ref
        lastal -s 2 -T 0 -Q 0 -r 1 -a 1 -b 1 ${true_ref%.fa} $fp > $fo.true.maf 
        maf-convert.py -d sam $fo.true.maf | $samtools view -bT $true_ref - | $samtools sort - $fo.last.true.srt # -d switch added to get @SQ header lines
        $samtools index $fo.last.true.srt.bam
	$samtools mpileup -f $true_ref $fo.last.true.srt.bam | $pileup_parser > $fo.last.true.pileup
        #Per read stats
        cat $fo.true.maf | ./$parse_maf > $fo.last.true.read.stats
}

analyze_output(){
#Analyze the output data
echo "library(ggplot2)
library(outliers)
pdf(file='$fo.pdf')

#Function to read in FASTQ files
read.fa<-function(fp){
	dat<-scan(fp, sep='\n',what='character')
	out<-data.frame(name=c(), nucs=c(), length=c(), stringsAsFactors=FALSE)
        while(length(dat)>0){
                #If the first line doesn't mark a new sequence, delete the first line and try again
                if(substr(dat[1],1,1)!='>'){
                        dat<-dat[-1]
                }else{
                        title=substr(dat[1],2,nchar(dat[1])) #Name of sequence
                        dat<-dat[-1]
			nucs=''
			while(substr(dat[1],1,1)!='>'){
				nucs<-paste(nucs, dat[1],sep='')
				dat<-dat[-1]
				if(length(dat)==0){break}
			}
			out<-rbind(out, data.frame(name=title, nucs=nucs, length=nchar(nucs), stringsAsFactors=FALSE))
                }
        }
	return(out)
}

read.blasr.stats<-function(fp, keep.top=TRUE){
	dat<-read.table(file=fp, sep=' ', stringsAsFactors=FALSE, header=F)
	names(dat)<-strsplit('qName tName score percentSimilarity qStrand qStart qEnd qLength tStrand tStart tEnd tLength mapQV', ' ')[[1]] #BLASR -m 4 format
	dat\$readname<-sapply(dat\$qName, function(x){strsplit(x, '/')[[1]][1]})
	if(keep.top){
		#Only keep the best alignment for each read
		dat<-do.call(rbind, lapply(unique(dat\$readname), function(r){
			t.dat<-subset(dat, readname==r)
			if(nrow(t.dat)==1){return(t.dat)}
			t.dat<-subset(t.dat, score==max(t.dat\$score))
			return(t.dat)
		}))
	}
	dat<-dat[rev(order(dat\$score)),]
	dat\$align_len<-abs(dat\$tEnd-dat\$tStart)
	#Calculate the aligned read length X percentSimilarity
	dat\$weighted_len<-dat\$align_len*dat\$percentSimilarity

	return(dat)
}

read.last.stats<-function(fp, keep.top=TRUE){
	dat<-read.table(file=fp, sep='\t', stringsAsFactors=FALSE, header=T)
	dat\$readname<-sapply(dat\$readname, function(x){strsplit(x, '/')[[1]][1]})
	if(keep.top){
		#Only keep the best alignment for each read
		dat<-do.call(rbind, lapply(unique(dat\$readname), function(r){
			t.dat<-subset(dat, readname==r)
			if(nrow(t.dat)==1){return(t.dat)}
			t.dat<-subset(t.dat, score==max(t.dat\$identities))
			return(t.dat)
		}))
	}
	dat<-dat[rev(order(dat\$identities)),]
	dat\$percentSimilarity<-100*dat\$identities/(dat\$identities+dat\$insertions+dat\$deletions+dat\$mismatches)
	dat\$qLength<-dat\$readlen
	dat\$tStart<-dat\$refstart
	dat\$tEnd<-dat\$refstart+dat\$align_len
	#Calculate the aligned read length X percentSimilarity
	dat\$weighted_len<-dat\$align_len*dat\$percentSimilarity
	dat\$tName<-dat\$refname
	
	return(dat)
}

print.aligned.reads<-function(dat, pileup.fp, type=''){
	main='Similarity to Known Reference by Read Length'
	if(type != ''){main=paste(main, type)}
	print(qplot(data=dat, align_len, percentSimilarity, main=main, geom='point',
			xlab='Aligned Read Length', ylab='Percent Similarity')+theme_bw())
	#Get the lengths of the input reads from the FASTQ
	reads<-read.fa('$fp')

	#Mark each read as one that aligned or not
	reads\$Aligned<-sapply(reads\$name, function(x){x %in% dat\$readname})

	#Plot the read lengths according to whether or not they aligned
	main='Read Length Distribution by Alignment to Reference'
	if(type != ''){main=paste(main, type)}
	print(qplot(data=reads, length, geom='histogram', xlab='Read Length', ylab='Number of Reads',
		    main=main)+theme_bw()+facet_grid(Aligned~., scales='free_y'))

	#For the reads that aligned, plot the read length against the alignment length
	main='Read Length vs. Alignment Length'
	if(type != ''){main=paste(main, type)}
	print(qplot(data=dat, qLength, align_len, xlab='Read Length', ylab='Alignment Length', 
		    main=main, geom='point')+theme_bw()) #+geom_abline(slope=1, intercept=0))

	#Plot the length and quality of each read along the true reference
	read.pos<-do.call(rbind, lapply(1:nrow(dat), function(q){
		data.frame(Read=dat\$readname[q], Percent_Identity=dat\$percentSimilarity[q], Position=dat\$tStart[q]:dat\$tEnd[q])
	}))
	#Don't plot more than 40 reads, so that the labels remain visible
	max.reads<-40
	if(length(unique(read.pos\$Read))>max.reads){read.pos<-subset(read.pos, Read %in% sample(unique(read.pos\$Read), max.reads))}
	main='Read Identity and Position'
	if(type != ''){main=paste(main, type)}
	print(qplot(data=read.pos, Position, Read, geom='line', color=Percent_Identity, size=5, xlab='Amplicon Position', main=main)
		+theme_bw()+guides(size=FALSE)+scale_color_continuous(name='Percent Identity'))

	#Show the coverage along the length of the reference
	pileup<-read.table(file=pileup.fp, sep='\t', fill=T, stringsAsFactors=FALSE,header=T)
	main='Coverage of Reference'
	if(type != ''){main=paste(main, type)}
	print(qplot(data=pileup, pos, cov, xlab='Position', ylab='Depth', geom='bar', stat='identity', main=main)+theme_bw())

}

## ALIGNMENT TO KNOWN TRUE REFERENCE

#What is the true accuracy of each read?
#Show the reads aligned to the true reference by BLASR
print.aligned.reads(read.blasr.stats('$fo.blasr.true.read.stats'), '$fo.blasr.true.pileup', type='(BLASR)')
#Show the reads aligned to the true reference by LAST
print.aligned.reads(read.last.stats('$fo.last.true.read.stats'), '$fo.last.true.pileup', type='(LAST)')


## ALIGNMENT TO LARGE DATABASE OF REFFERENCES

print.reference.scores<-function(dat, type=''){
	#Reformat the names (because greengenes gave it something much too long)
	dat\$tName<-gsub('\"', '', dat\$tName)
	dat\$tName<-gsub(\"'\", '', dat\$tName)
	dat\$tName<-sapply(dat\$tName, function(x){x<-strsplit(x, '_')[[1]];if(substr(x[2],1,2)=='gi'){x<-x[-2]};return(paste(x, collapse='_'))})
	dat\$tName<-sapply(dat\$tName, function(x){y<-strsplit(x, '_')[[1]]; if(length(y)<10){return(x)}; y<-y[7:(length(y)-3)];return(paste(y, collapse=' '))})
	dat\$tName<-sapply(dat\$tName, function(x){if(substr(x,1,1) %in% as.character(c(1:10))){y<-strsplit(x, '_')[[1]]; return(paste(y[-1], collapse='_'))};return(x)})
	dat\$tName<-gsub('_', ' ', dat\$tName)

	#Keep a single alignment for each read against each reference, choosing the alignment with the highest weighted length (align_len * percentSimilarity)
	dat<-do.call(rbind, lapply(unique(dat\$readname), function(r){
		t.dat<-subset(dat, readname==r, stringsAsFactors=FALSE)
		if(nrow(t.dat)==1){return(t.dat)}
		return(do.call(rbind, lapply(unique(t.dat\$tName), function(Ref){
			t.t.dat<-subset(t.dat, tName==Ref, stringsAsFactors=FALSE)
			if(nrow(t.t.dat)==1){return(t.t.dat)}
			return(subset(t.t.dat, weighted_len==max(t.t.dat\$weighted_len), stringsAsFactors=FALSE))
		})))
	}))

	#Count the number of aligned reads for each reference
	ref<-data.frame(name=unique(dat\$tName), stringsAsFactors=FALSE)
	ref\$nreads<-sapply(ref\$name, function(r){sum(dat\$tName==r)})
	ref\$alignedBases<-sapply(ref\$name, function(r){sum(dat\$align_len[dat\$tName==r])})
	#Score based on average percent identity of the aligned reads
	ref\$percIden<-sapply(ref\$name, function(r){sum(dat\$percentSimilarity[dat\$tName==r])/sum(dat\$tName==r)})
	#Score based on aligned reads and identity
	ref\$score<-ref\$nreads * ref\$percIden
	#ref\$score<-sapply(ref\$name, function(r){sum(dat\$weighted_len[dat\$tName==r])})

	print(qplot(data=ref, nreads, alignedBases, geom='point', 
			main=paste('Number of Reads and Number of Bases Aligned by Reference', type), 
			xlab='Number of Reads Aligned', ylab='Number of Bases Aligned')+theme_bw())
	ref.hit.plot(data=ref[,c('name','nreads')], main=paste('Number of Reads Aligned', type), xlab='', ylab='# Reads')
	ref.hit.plot(data=ref[,c('name','alignedBases')], main=paste('Number of Bases Aligned', type), xlab='', ylab='# Bases')
	ref.hit.plot(data=ref[,c('name','score')], main=paste('Weighted Alignment Score', type), xlab='', ylab='Score')
	ref.hit.plot(data=ref[,c('name','percIden')], main=paste('Mean Percent Identity', type), xlab='', ylab='Percent Identity')
}

ref.hit.plot<-function(data='', xlab='', ylab='', main=''){
	names(data)<-c('name','score')
	#Order the table by that score
	data<-data[rev(order(data\$score)),]
	data<-data[min(nrow(data), 20):1,]

	print(data)
	data\$name<-factor(data\$name, levels=data\$name, ordered=TRUE)

	print(qplot(data=data, name, score, geom='bar', stat='identity', main=main, xlab=xlab, ylab=ylab)+coord_flip()+theme_bw())
}

#Read in the stats on the aligned reads
dat<-read.blasr.stats('$fo.blasr.read.stats', keep.top=FALSE)
print.reference.scores(dat, type='(BLASR)')
dat<-read.last.stats('$fo.last.read.stats', keep.top=FALSE)
print.reference.scores(dat, type='(LAST)')

dev.off(); q()
" | R --no-save --no-restore
}

extract_reads
run_blasr
run_last
analyze_output
