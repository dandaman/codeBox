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
parser.add_option("-d","--dom2jalview",dest="domains",help="write domain architecture as jalview feature annotation file", default=None,metavar="JALVIEW_FEAUTURE_FILE")
(options, args) = parser.parse_args()

xml	= PhyloXMLIO.read(args[0])

if options.domains:
    domfile=open(options.domains,"w")
for t in xml.phylogenies:
    for n in t.get_terminals():
        seq=n.sequences[0].to_seqrecord()
        seq.description=""
        if not options.id:
            seq.id=n.name
        seq.name=n.name
        SeqIO.write(seq,sys.stdout,"fasta")
        if options.domains and n.sequences[0].domain_architecture:
            for d in n.sequences[0].domain_architecture.domains:
                domfile.write("\t".join(["from phyloxml",n.name,str(-1),str(d.start),str(d.end),d.value])+"\n")
