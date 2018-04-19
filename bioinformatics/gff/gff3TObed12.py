#!/usr/bin/env python
import csv
import collections
import Bio.SeqIO
import argparse
import sys

#just adapted this one to my needs
##initial script my Michael Seidel

parser = argparse.ArgumentParser(description='Converts GFF3 to BED12 (12 column BED format). See https://genome.ucsc.edu/FAQ/FAQformat.html#format1 and https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md for details.')
parser.add_argument('-i','--in', dest='gff3',type=str, help='GFF3 input file name',required=True)
parser.add_argument('-o', '--out', dest='bed12',type=str,  help='BED12 output file name',required=True)

args = parser.parse_args()


fhd = open(args.gff3)

gff = {}

c = csv.reader(filter(lambda r: not r.startswith('#'), fhd), delimiter='\t')
rec = collections.namedtuple("record", "feature chrom start end strand")

ids = set()

for row in c:
	if row[2] =="gene":
		continue

	try:
		info = dict((x.split('=', 1) for x in filter(None, row[-1].split(';')) if '=' in x ))
	except ValueError:
		print(row)
	k = info['ID'] if row[2] == 'mRNA' else info['Parent']
	gff.setdefault(k, []).append(rec(row[2], row[0], int(row[3]), int(row[4]), row[6]))
	ids.add(k)

def get_line(gid):
	t = gff[gid]
	mrna = list(filter(lambda x: x.feature == 'mRNA', t))[0]
	# mrna.chrom, mrna.start-1, mrna.end, mrna.strand

	exons = list(filter(lambda x: x.feature == 'exon', t))
	blockSizes = [x.end-x.start+1 for x in exons]
	blockStarts = [x.start-mrna.start for x in exons]

	cds = [(x.start, x.end) for x in list(filter(lambda x: x.feature == 'CDS', t))]
	cds = sorted([x for y in cds for x in y])
	try:
		thickStart = cds[0]-1
		thickEnd = cds[-1]
	except IndexError:
		print(gid)
		raise IndexError

	return [mrna.chrom, mrna.start-1, mrna.end, gid.split('.')[0], 1000, mrna.strand, thickStart, thickEnd,
			0, len(blockSizes), ','.join(list(map(str, blockSizes))), ','.join(list(map(str, blockStarts)))]



with open(args.bed12, 'w') as out:
	for gid in ids:
		out.write('\t'.join(list(map(str, get_line(gid)))) + '\n')
