#!/usr/bin/python

#Parse the MAF output for per-read alignment statistics

#Output will have the columns:
#score refname	ref_strand	reflen	refstart	refend	readname	readlen	readstrand	readstart	readend	align_len	identities	mismatches	insertions	deletions

import sys

ref=''
read=''

def print_read_stat(ref, read, score):
	if len(ref)==0: return()
	refname=ref[1]
	ref_strand=ref[4]
	reflen=ref[5]
	refstart=ref[2]
	if ref_strand=='+': refend=str(int(refstart)+int(reflen))
	if ref_strand=='-': refend=str(int(refstart)-int(reflen))
	readname=read[1]
	readlen=read[5]
	read_strand=read[4]
	readstart=read[2]
	if read_strand=='+': readend=str(int(readstart)+len(ref[6]))
	if read_strand=='-': readend=str(int(readstart)-len(ref[6]))
	#Now focus on the actual alignment
	ref=ref[6]
	read=read[6]
	align_len=len(ref)
	identities=str(sum([1 for q in range(0,align_len) if ref[q] == read[q]]))
	mismatches=str(sum([1 for q in range(0,align_len) if ref[q] != read[q] and ref[q]!='-' and read[q]!='-']))
	insertions=str(sum([1 for q in range(0,align_len) if ref[q] == '-']))
	deletions=str(sum([1 for q in range(0,align_len) if read[q] == '-']))
	align_len=str(align_len)
	print('\t'.join([score, refname, ref_strand, reflen, refstart, refend, readname, readlen, read_strand, readstart, readend, align_len, identities, mismatches, insertions, deletions]))

#Header
print '\t'.join(['score','refname','ref_strand','reflen','refstart','refend','readname','readlen','readstrand','readstart','readend','align_len','identities','mismatches','insertions','deletions'])

for line in sys.stdin:
	if line[0] == '#': continue
	line=line.strip()
	if line.startswith('a score='):
		score=line[8:]
		print_read_stat(ref, read, score)
		ref=''; read=''; score=''
		continue
	if len(line)==0: continue
	line=line.split(' ')
	line=[x for x in line if len(x)>0] #Remove empty positions
	if len(ref)==0: ref=line; continue
	if len(read)==0: read=line; continue
