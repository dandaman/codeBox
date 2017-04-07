#!/usr/bin/env python
from Bio import SeqIO
import sys

out_cds=open(sys.argv[2], "w")  
out_protein=open(sys.argv[3], "w")  

for index,record in enumerate(SeqIO.parse(sys.argv[1],"genbank")):
	for f in record.features:
		if f.type =="CDS":
			cds=f.extract(record)
			if "locus_tag" in f.qualifiers:
				cds.id=f.qualifiers["locus_tag"][0]
			elif "gene" in f.qualifiers:
				cds.id=f.qualifiers["gene"][0]
			cds.description=f.qualifiers["product"][0]
			SeqIO.write([cds],out_cds,"fasta")
			protein=cds
			protein.seq=cds.seq.translate()
			SeqIO.write([protein],out_protein,"fasta")
			
			

