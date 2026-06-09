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
    tryCatch({
      res <- entrez_search(db="gene", term=paste("(", i, "[GENE] AND Mycobacterium tuberculosis[ORGN] NOT discontinued[Properties])"))
      esums <- entrez_summary(db="gene", id=res$ids)
      geneOrg <- c(geneOrg, esums$organism$scientificname)
    }, error = function(e) {
      geneOrg <<- c(geneOrg, NA)
    })
  }
  return(geneOrg)
}

# get gene ID from given gene names
GeneID <- function(names){
  id = names
  geneIDs = c()
  for (i in id){
    tryCatch({
      res <- entrez_search(db="gene", term=paste("(", i, "[GENE] AND Mycobacterium tuberculosis H37Rv [ORGN] NOT discontinued[Properties])"))
      if (res$count == 0){
        geneIDs <- c(geneIDs, i)  # gene not found — keep original ID
      } else {
        esums <- entrez_summary(db="gene", id=res$ids)
        # otheraliases can be NULL if NCBI has no aliases for this gene,
        # which causes strsplit() to return list() and [[1]] to crash.
        # Guard: fall back to the original input ID on any lookup failure.
        aliases_raw <- esums$otheraliases
        if (!is.null(aliases_raw) && length(aliases_raw) > 0 && nchar(as.character(aliases_raw)[1]) > 0) {
          alias_parts <- strsplit(as.character(aliases_raw)[1], ",")
          gene_alias <- if (length(alias_parts) > 0 && length(alias_parts[[1]]) > 0)
                          trimws(alias_parts[[1]][1])
                        else i
        } else {
          gene_alias <- i  # no aliases available — keep original ID
        }
        geneIDs <- c(geneIDs, gene_alias)
      }
    }, error = function(e) {
      # API unavailable, rate-limited, or unexpected format — keep original ID
      geneIDs <<- c(geneIDs, i)
    })
  }
  return(geneIDs)
}


#Get gene names from given gene IDs
GeneName <- function(id){
  id = id
  geneName = c()
  for (i in id){
    tryCatch({
      res <- entrez_search(db="gene", term=paste("(", i, "[GENE] AND Mycobacterium tuberculosis[ORGN] NOT discontinued[Properties])"))
      esums <- entrez_summary(db="gene", id=res$ids)
      geneName <- c(geneName, esums$name)
    }, error = function(e) {
      geneName <<- c(geneName, NA)
    })
  }
  return(geneName)
}

# Get description from given gene IDs
Description <- function(id){
  id = id
  desc = c()
  for (i in id){
    tryCatch({
      res <- entrez_search(db="gene", term=paste("(", i, "[GENE] AND Mycobacterium tuberculosis[ORGN] NOT discontinued[Properties])"))
      esums <- entrez_summary(db="gene", id=res$ids)
      extracteddesc <- extract_from_esummary(esums, "description")
      desc <- c(desc, extracteddesc)
    }, error = function(e) {
      desc <<- c(desc, NA)
    })
  }
  return(desc)
}

