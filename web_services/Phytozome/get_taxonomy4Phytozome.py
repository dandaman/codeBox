from intermine.webservice import Service

from ete3 import NCBITaxa
ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()

service = Service("https://phytozome.jgi.doe.gov/phytomine/service")
query = service.new_query("Organism")
query.add_view(
    "annotationVersion", "assemblyVersion", "commonName", "genus", "name",
    "proteomeId", "shortName", "species", "taxonId", "version"
)
k=["annotationVersion", "assemblyVersion", "commonName", "genus", "name", "proteomeId", "shortName", "species", "taxonId", "version"]
t=["superkingdom", "kingdom", "phylum", "class", "subclass", "order", "family", "genus","species"]
print("\t".join(k+t+["full_lineage"]))

def filterRanks(L):
	subset={ncbi.get_rank([x])[x]:x for x in L}
	#return([if x in subset: ncbi.get_taxid_translator([x])[x] else: "NA" for x in t])
	return([list(ncbi.get_taxid_translator([subset[x]]).values())[0] if x in subset else 'NA' for x in t])


for row in query.rows():
	lineage=ncbi.get_lineage(row["taxonId"])
	names=ncbi.get_taxid_translator(lineage)
	print ("\t".join([str(row[x]) for x in k] + filterRanks(lineage) + [",".join(names[z] for z in lineage)]))
