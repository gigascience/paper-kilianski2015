#!/bin/sh

#Run the analysis for every set of samples in this folder, using another script, minION_amplicon.sh, for the heavy lifting

./minION_amplicon.sh MVA fast5/MVA poxvirus_ref.fa MVA_ref.fa
./minION_amplicon.sh cowpox fast5/cowpox poxvirus_ref.fa cowpox_ref.fa
./minION_amplicon.sh lister fast5/lister poxvirus_ref.fa lister_ref.fa
./minION_amplicon.sh ecoli fast5/ecoli 16S_ref.fa ecoli_ref.fa
