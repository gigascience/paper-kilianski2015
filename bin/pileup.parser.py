#! /usr/bin/env python
###############################################################################
# parses pileup base string and returns the counts for all possible alleles
# for each position
# reads input (mpileup output) from sys.stdin
###############################################################################

import os
import sys

def is_int(s):
    try:
	int(s)
	return True
    except ValueError:
	return False

class parseString(object):
    
    def __init__(self, ref, string):
        self.ref = ref.upper()
        self.string = string.upper()
        self.types = {'A':0,'G':0,'C':0,'T':0,'-':[],'*':0,'+':[],'X':[]}
        self.process()
        
    def process(self):
        # remove end of read character
        self.string = self.string.replace('$','')
        while self.string != '':
            if self.string[0] == '^':
                # skip two characters when encountering '^' as it indicates
                # a read start mark and the read mapping quality
                self.string = self.string[2:]
            elif self.string[0] == '*':
                self.types['*'] += 1
                # skip to next character
                self.string = self.string[1:]
            elif self.string[0] in ['.',',']:
                self.types[self.ref] += 1
                self.string = self.string[1:]
            elif self.string[0] == '+':
                if is_int(self.string[2]): #Two digit insertion
                    insertionLength = int(self.string[1:3])
                    insertionSeq = self.string[3:3+ insertionLength]
                    self.string = self.string[3+insertionLength:]
                else:
                    insertionLength = int(self.string[1])
                    insertionSeq = self.string[2:2+ insertionLength]
                    self.string = self.string[2+insertionLength:]
                self.types['+'].append(insertionSeq)
            elif self.string[0] == '-':
                if is_int(self.string[2]): #Two digit deletion
                    deletionLength = int(self.string[1:3])
                    deletionSeq = self.string[3:3+deletionLength]
                    self.string = self.string[3+deletionLength:]
                else:
                    deletionLength = int(self.string[1])
                    deletionSeq = self.string[2:2+deletionLength]
                    self.string = self.string[2+deletionLength:]
                self.types['-'].append(deletionSeq)
            else:
                # unrecognized character or substitution
		if self.types.has_key(self.string[0]):
	                self.types[self.string[0]] += 1
		else:
	                self.types['X'].append(self.string[0])
                self.string = self.string[1:]
        return
    def __repr__(self):
        types = self.types
        return '\t'.join(map(str,[types['A'], types['C'], types['G'],types['T'],\
                                  types['*']]) +\
                         map(','.join, [types['-'],types['+'],types['X']]))
        

def main():
    print >>sys.stdout, "chrom\tpos\tref\tcov\tA\tC\tG\tT\t*\t-\t+\tX"
    for line in sys.stdin:
        toks = line.strip('\n').split('\t')
        ref = toks[2].upper()
        cov = toks[3]
        print >>sys.stdout, '\t'.join([toks[0], toks[1],ref, cov]) + '\t' + \
            parseString(ref, toks[4]).__repr__()

if __name__ == '__main__':
    main()
