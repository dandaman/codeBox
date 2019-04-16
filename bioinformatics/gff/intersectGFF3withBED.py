#!/usr/bin/env python

from bx.intervals.intersection import IntervalTree
import pandas as pd
import sys
import argparse

parser=argparse.ArgumentParser(description='Intersects a GFF3 file with regions given in a BED feature file')
parser.add_argument("GFF3", metavar="GFF3_file", type=str, help="path to GFF3 file")
parser.add_argument("BED", metavar="BED_file", type=str, help="path to BED file")
parser.add_argument("output", metavar="OUTPUT_file", type=str, help="path to filtered output GFF3 file",nargs="?",default=sys.stdout)
parser.add_argument("-t","--type", dest="tag",metavar="type", type=str, help="feature type i.e primary tag in GFF3 to return", default=None)
parser.add_argument("-a","--attribute",dest="attribute", metavar="ATTRIBUTE" ,help="Do not return GFF3. Return value of specified attribute filed of 8th column",default=None)
args = parser.parse_args()

def index_gff3(gff3_file_path):
	#following an example from https://malariageninformatics.wordpress.com/2011/07/07/using-interval-trees-to-query-genome-annotations-by-position/
	 # dictionary mapping chromosome names to interval trees
	genome = dict()

	# parse the annotations file (GFF3) and build the interval trees
	gff=pd.read_csv(gff3_file_path,sep="\t",header=None,comment="#")
	for idx,row in gff.iterrows():
		if args.tag is not None and row[2] != args.tag:
			continue
		seqid = row[0]
		start = int(row[3])
		end = int(row[4])
		tree = None
		# one interval tree per chromosome
		if seqid in genome:
			tree = genome[seqid]
		else:
			# first time we've encountered this chromosome, create an interval tree
			tree = IntervalTree()
			genome[seqid] = tree
		# index the feature
		if args.attribute is None:
			tree.add(start, end, row)
		else :
			attr=row[8].split(";")
			o=list()
			for n in attr:
				k,v=n.split("=")
				if k == args.attribute:
					o.append(v)
			o=",".join(o)
			tree.add(start, end, o)
	return genome

genome = index_gff3(args.GFF3)
d=pd.read_csv(args.BED,sep="\t",header=None)

D=d.apply(lambda x: genome[x[0]].find(x[1]+1,x[2]),axis=1)

o=list()
for g in D.tolist():
	o.append(pd.DataFrame(g))
o=pd.concat(o)
o.to_csv(args.output,sep="\t",index=False,header=False)
