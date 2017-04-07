#!/bin/env python
from Bio import Phylo
from Bio.Phylo import PhyloXMLIO
from Bio import SeqIO
from optparse import OptionParser
import sys

usage	= "%prog [options] input_phyloxml_file" 
parser = OptionParser(usage=usage)

#parser.add_option("-f","--format",dest="format",help="seqio format to export sequences into", default="fasta",metavar="SEQIOFORMAT")
parser.add_option("-i","--id",dest="id",help="use id instead of node.name for display id", default=False,action="store_true")
(options, args) = parser.parse_args()

xml	= PhyloXMLIO.read(args[0])

for t in xml.phylogenies:
	for n in t.get_terminals():
		seq=n.sequences[0].to_seqrecord()
		seq.description=""
		if not options.id:
			seq.id=n.name
		seq.name=n.name
		SeqIO.write(seq,sys.stdout,"fasta")
