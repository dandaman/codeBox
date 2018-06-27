#!/home/pgsb/daniel.lang/anaconda3/bin/Rscript

args = commandArgs(trailingOnly=TRUE)


D=args[1]
o=basename(args[2])

d=gsub("^input/","",D)
d=gsub("_sets/?","",d)
dir.create(sprintf("analysis/%s",d))
print(sprintf("Processing %s",D))
for (O in o) {
	a=NULL
	dir.create(sprintf("analysis/%s/%s",d,O))
	dir.create(sprintf("analysis/%s/%s/sets",d,O))
	dir.create(sprintf("analysis/%s/%s/out",d,O))

	file.symlink(sprintf("../../../annotation/Triticum_aestivum_TGAC_PGSB_%sA_by_orthology.gaf",O), sprintf("analysis/%s/%s/input.gaf",d,O))
	file.symlink(sprintf("../../../Ontologies/%s.obo",O), sprintf("analysis/%s/%s/ontology.obo",d,O))
	file.symlink(sprintf("../../../Snakefile.sub",O), sprintf("analysis/%s/%s/Snakefile",d,O))

	if (file.exists(sprintf("./Annotation/%s.RData",O))) {
		load(sprintf("./Annotation/%s.RData",O))
	} else {
		a<-read.table(paste("./annotation/Triticum_aestivum_TGAC_PGSB_",O,"A_by_orthology.gaf",sep=""),sep="\t",stringsAsFactors=FALSE)
		a<-unique(a[,2])
		save(a,file=sprintf("./Annotation/%s.RData",O))
	}

	FF=unique(unlist(sapply(list.files(D,pattern="*.set"),function(f) {
		F=sprintf("%s/%s",D,f)
		G=unique(read.table(F,sep="\t",stringsAsFactors=FALSE,header=FALSE)[,1])

		dname=sprintf("analysis/%s/%s/sets",d,O)
		oname=sprintf("analysis/%s/%s/out",d,O)

		dir.create(dname)
		dir.create(oname)
		y=G[G %in% a] 
		if (length(y)>0) {
			x=gsub("\\.set","",f)
			write.table(y,file=paste(dname,paste(x,"set",sep="."),sep="/"),sep="\t", row.names=FALSE,col.names=FALSE,quote=FALSE)
			return(y)
		}
	})))
	write.table(unlist(FF),file=sprintf("analysis/%s/%s/pop",d,O), row.names=FALSE,col.names=FALSE,quote=FALSE)
}

