
bgee = read.csv("Downloads/CTD_exposure_events.csv"
                , skip=5000)



cgd = read.table("Downloads/curated_gene_disease_associations.tsv"
                 , sep = "\t"
                 ,skip = 1000)


.pKG.edges = read.csv("Downloads/edges.csv")
dim(.pKG.edges)
head(.pKG.edges)

set.seed(1)
.indx = sample(1:nrow(.pKG.edges), size=1e5, replace=F)
.tester = .pKG.edges[.indx,]
head(.tester)

apply(.tester, 2, function(a) length(unique(a)))
nrow(.tester)


library(readr)
.kg = read_csv("Downloads/kg_giant.csv")

head(.kg)