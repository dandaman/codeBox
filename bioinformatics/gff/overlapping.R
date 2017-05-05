#!/usr/bin/env R 
library("GenomicRanges")
require("rtracklayer")



args 	= commandArgs(trailingOnly=TRUE)
ref		= as.character(args[1]) 
query	= as.character(args[2])
mincov	= as.numeric(args[3])

ref		= import.gff(ref,format="gff3")
query	= import.gff(query,format="gff3")

hits 	= findOverlaps(ref, query)
overlaps= pintersect(ref[queryHits(hits)], query[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(query[subjectHits(hits)])
hits <- hits[percentOverlap >= mincov]

match_hit <- data.frame(ref@elementMetadata$ID[queryHits(hits)], query@elementMetadata$ID[subjectHits(hits)],stringsAsFactors=F)
write.table(match_hit,file="",sep="\t", quote=FALSE,row.names=FALSE,col.names=FALSE)
