#!/bin/python3

import sys
import csv

IN_FREQ_FILE=sys.argv[1]
OUT_FILE=sys.argv[2]
out_no_trans=csv.writer(open(OUT_FILE,'w'), delimiter='\t')
transversion=[{'A','C'},{'A','T'},{'G','C'},{'G','T'}]

with open(IN_FREQ_FILE) as tsv:
	for line in csv.reader(tsv, dialect="excel", delimiter='\t'):
		alleles=set(line[2]+line[3])
		if alleles in transversion:
			out_no_trans.writerow(line)
