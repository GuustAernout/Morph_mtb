# Character vector of installed packages
pkg <- installed.packages()[, "Package"]
# Check if rentrez package is installed
'rentrez' %in% pkg
# If package is not installed --> install!
if(!('rentrez' %in% pkg)) {install.packages("rentrez")}
library("rentrez")

# get organism of given gene IDs
Organism <- function(id){
  id = id
  geneOrg = c()
  for (i in id){
    res <- entrez_search(db="gene", term=paste("(", i, "[GENE] AND Mycobacterium tuberculosis[ORGN] NOT discontinued[Properties])"))
    res$count
    res$ids
    esums <- entrez_summary(db="gene", id=res$ids)
    geneOrg <- c(geneOrg, esums$organism$scientificname)
  }
  return(geneOrg)
}

# get gene ID from given gene names
GeneID <- function(names){
  id = names
  geneIDs = c()
  for (i in id){
    res <- entrez_search(db="gene", term=paste("(", i, "[GENE] AND Mycobacterium tuberculosis H37Rv [ORGN] NOT discontinued[Properties])"))
    res$count
    if (res$count == 0){
      geneIDs <- c(geneIDs, "FALSE")
    } else {
    res$ids
    esums <- entrez_summary(db="gene", id=res$ids)
    geneIDs <- c(geneIDs, strsplit(as.character(esums$otheraliases), ",")[[1]][1])
    }
}
  return(geneIDs)
}


#Get gene names from given gene IDs
GeneName <- function(id){
  id = id
  geneName = c()
    for (i in id){
      res <- entrez_search(db="gene", term=paste("(", i, "[GENE] AND Mycobacterium tuberculosis[ORGN] NOT discontinued[Properties])"))
      res$count
      res$ids
      esums <- entrez_summary(db="gene", id=res$ids)
      geneName <- c(geneName, esums$name)
    }
  return(geneName)
}

# Get description from given gene IDs
Description <- function(id){
  id = id
  desc = c()
    for (i in id){
      res <- entrez_search(db="gene", term=paste("(", i, "[GENE] AND Mycobacterium tuberculosis[ORGN] NOT discontinued[Properties])"))
      res$count
      res$ids
      esums <- entrez_summary(db="gene", id=res$ids)
      extracteddesc <- extract_from_esummary(esums, "description")
      desc <- c(desc, extracteddesc)
    }
  return(desc)
}

