#!/home/pgsb/daniel.lang/anaconda3/bin/Rscript
library(stringr)
library(ontologyIndex)
library(ontologySimilarity)

args = commandArgs(trailingOnly=TRUE)

to=get_ontology("/nfs/pgsb/data/evograph/annotation/ontologies/to.obo", propagate_relationships=get_relation_names("/nfs/pgsb/data/evograph/annotation/ontologies/to.obo"))
po=get_ontology("/nfs/pgsb/data/evograph/annotation/ontologies/po.obo", propagate_relationships=get_relation_names("/nfs/pgsb/data/evograph/annotation/ontologies/po.obo"))
go=get_ontology("/nfs/pgsb/data/evograph/annotation/ontologies/go.obo", propagate_relationships=get_relation_names("/nfs/pgsb/data/evograph/annotation/ontologies/go.obo")[-c(4)])

toic=descendants_IC(to)
poic=descendants_IC(po)
goic=descendants_IC(go)


analysis_path=args[1]
method=args[2]
cutoff=as.numeric(args[3])
output=args[4]

pat=sprintf("*-Parent-Child-Union-%s.txt",method)
O=do.call(rbind.data.frame,lapply(dir(analysis_path,pattern="^[^.]") , function(sub) {
	DIR=paste(analysis_path,sub,sep="/")
	o=do.call(rbind.data.frame,lapply(dir(DIR,pattern="^[^.]"),function(ontology) {
		DIRR=paste(DIR,ontology,sep="/")
		o=do.call(rbind.data.frame,lapply(list.files(DIRR,pattern=pat,recursive=TRUE),function(f) {
			d=read.delim(file=paste(DIRR,f,sep="/"),sep="\t",stringsAsFactors=FALSE)
			m=str_match(f,"table-(.+)-Parent")
			sset=m[2]
			if (!is.na(sset)) {
				d=subset(d,p.adjusted<cutoff)
				f=gsub("table-","anno-",f)
				a=read.delim(paste(DIRR,f,sep="/"),stringsAsFactors=FALSE,fill=TRUE,header=FALSE)
				a=subset(a,!is.na(a[,3]) & ! a[,3] =="")
				s=read.table(sprintf("%s/sets/%s.set",DIRR,sset),stringsAsFactors=FALSE)[,1]
				a=subset(a,a[,1] %in% s)
				a=do.call(rbind.data.frame, lapply(1:nrow(a), function(x) data.frame(locus=a[x,1],ID=unlist(strsplit(gsub("\\}","",gsub("[^={ ,]+=\\{","",a[x,3])),",")))))
				a=as.data.frame(aggregate(a$locus,by=list(a$ID),paste,collapse=","))
				names(a)=c("ID","loci")
				d=merge(d,a,by="ID")
				if (nrow(d)>0) {
					data.frame(group=sub,set=sset,method=method,cutoff=cutoff,ontology=ontology,d)
				}
			}
		}))
		o
	}))
	o
}))

OO=list(GO=go,PO=po,TO=to)
Oic=list(GO=goic,PO=poic,TO=toic)

n=grep("^(loci|name)$",names(O),perl=TRUE)
OK=O[,n]
O=O[,-n]

O=transform(O,ID=as.character(O$ID))

O$depth=sapply(1:nrow(O),function(i) length(OO[[O[i,"ontology"]]]$ancestors[[O[i,"ID"]]])-1)
O$ic=sapply(1:nrow(O),function(i) Oic[[O[i,"ontology"]]][[O[i,"ID"]]])
O$partition=sapply(1:nrow(O),function(i) {AO=OO[[O[i,"ontology"]]]; a=as.character(AO$ancestors[[O[i,"ID"]]]); ifelse(length(a)>=1, AO$name[[a[1]]][1],"NA")})
O$level3=sapply(1:nrow(O),function(i) {AO=OO[[O[i,"ontology"]]]; a=as.character(AO$ancestors[[O[i,"ID"]]]); a=a[-length(a)]; ifelse(length(a)>=3, AO$name[[a[3]]][1],"NA")})

O=data.frame(O,OK)

write.table(O,file=output,sep="\t",row.names=FALSE,col.names=TRUE,quote=TRUE)

