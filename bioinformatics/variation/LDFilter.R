#!/usr/bin/env Rscript
library(gdsfmt)
library(SNPRelate)

#args: vcffile ld.threshold outprefix [referencename]
args <- commandArgs(trailingOnly = TRUE)

#Marlene's birthday
set.seed(27042012)


gds.file<-gsub("\\.vcf",".gds",args[1])

cat(sprintf("Processing vcf %s gds %s %f to %s %s\n\n", args[1],gds.file, as.numeric(args[2]), args[3],args[4] ))

if (! file.exists(gds.file) ) {
	snpgdsVCF2GDS(args[1], gds.file, method="biallelic.only")
#	snpgdsVCF2GDS(args[1], gds.file, method="copy.num.of.ref")
}

genofile <- snpgdsOpen(gds.file)

snpgdsSummary(genofile)

snpset <- snpgdsLDpruning(genofile, ld.threshold=as.numeric(args[2]))
snpset.id<-unlist(snpset)

gds2fasta <- function (gdsobj, pos.fn, snp.id = NULL, ref=NULL) {

    total.snp.ids <- read.gdsn(index.gdsn(gdsobj, "snp.id"))
    snp.ids <- total.snp.ids
    if (!is.null(snp.id)) {
        n.tmp <- length(snp.id)
        snp.id <- snp.ids %in% snp.id
        n.snp <- sum(snp.id)
        if (n.snp != n.tmp) 
            stop("Some of snp.id do not exist!")
        if (n.snp <= 0) 
            stop("No SNP in the working dataset.")
        snp.ids <- snp.ids[snp.id]
    }
    snp.idx <- match(snp.ids, total.snp.ids)

    rep.genotype <- read.gdsn(index.gdsn(gdsobj, "genotype"))[,snp.idx]
    rep.allele <- do.call(rbind, strsplit(read.gdsn(index.gdsn(gdsobj, "snp.allele"))[snp.idx], "/", fixed = TRUE))
    sample.id <- read.gdsn(index.gdsn(gdsobj, "sample.id"))

    id_file.name <- paste(pos.fn, ".id.txt", sep = "")
    cat(paste(read.gdsn(index.gdsn(gdsobj, "snp.rs.id"))[snp.idx], collapse = "\n"), "\n", file = id_file.name)

    seq.len <- length(snp.idx)
    file.name <- paste(pos.fn, ".fasta", sep = "")
    cat("", file = file.name) # Make a new empty file
    if (! is.null(ref)) {
        seq <- character(seq.len)
        for (j in 1:seq.len) {
		seq[j] <- rep.allele[j,1]
        }
	cat(">", ref, "\n", file = file.name, sep = "", append = TRUE)
        cat(seq, "\n", file = file.name, sep = "", append = TRUE)
    }

    for (i in 1:length(sample.id)) {
        seq <- character(seq.len)
        for (j in 1:seq.len) {
            if        (rep.genotype[i,j] == 0) {
                seq[j] <- rep.allele[j,1]
            } else if (rep.genotype[i,j]>0 ) {
                seq[j] <- rep.allele[j,(rep.genotype[i,j]+1)]
            } else {
                seq[j] <- "N"
            }
        }
        cat(">", sample.id[i], "\n", file = file.name, sep = "", append = TRUE)
        cat(seq, "\n", file = file.name, sep = "", append = TRUE)
    }
}
if (length(snpset.id)>0) { 
	if (length(args)==4) {
		gds2fasta(genofile,args[3],snp.id=snpset.id,ref=args[4])
	} else { 
		gds2fasta(genofile,args[3],snp.id=snpset.id)
	}
} else {
	cat(sprintf("No SNPs above threshold %f\n", as.numeric(args[2])))
}
