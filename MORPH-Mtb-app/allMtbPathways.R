##############################################################################################################
## allMtbPathways with MORPH-Mtb
##############################################################################################################

# Install and load required packages ----
required_packages <- c(
  "shiny", 
  "BiocManager",
  "shinyjs", 
  "waiter", 
  "shinythemes", 
  "readr", 
  "writexl", 
  "WriteXLS", 
  "readxl", 
  "limma", 
  "edgeR", 
  "plyr", 
  "dplyr", 
  "data.table", 
  "purrr", 
  "rsample", 
  "amap", 
  "factoextra", 
  "cluster", 
  "dtwclust", 
  "proxy", 
  "kohonen", 
  "caroline", 
  "DBI", 
  "dtw", 
  "ggplot2",
  "RMySQL",
  "MASS",
  "htmlTable"
)

# Check for missing packages and install them
missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if(length(missing_packages)) install.packages(missing_packages)
#update.packages(ask = FALSE)

# Character vector of installed packages
pkg <- installed.packages()[, "Package"]
# Check if packages are installed
# If package is not installed --> install!
if(!('limma' %in% pkg)) {BiocManager::install("limma", force=TRUE)}
if(!('edgeR' %in% pkg)) {BiocManager::install("edgeR", force=TRUE)}

# Load the required packages
lapply(required_packages, library, character.only = TRUE)

###############################################################################################

## MORPH ALGORITHM ##
# Functions
## Reading data functions
### [1.1] getClusteringInformation
getClusteringInformation = function(ClusterFile)
{
  Cluster = read.delim(ClusterFile, sep="\t", header=FALSE, row.names = 1,stringsAsFactors=F) #Read tab-delimited file.
  colnames(Cluster) = ClusterFile #Change cluster column name to file-name for future identification.
  return(Cluster)
}

### [1.2] Get Pathways
GetGOIs = function(GOI_File)
{
  GOI = readLines(GOI_File, n=1) #Reads the first (and only) line from the pathway-genes file.
  GOI = strsplit(GOI,"\t")[[1]] #Splits names into a vector by tabs.
  GOI = sub("^rv|^RV", "Rv", GOI)
  return(GOI)
}


### [1.3] getGeneExpression
getGeneExpression = function(InputGE)
{
  GE = (read.delim(InputGE, sep="\t", header=TRUE,row.names=1)) #Reads the gene-expression data from the tab-delimited file.
}


## MORPH Algorithm
### [2.1] prepare Morph Object From Files {Input;Configuration file (data and cluster solution), Output:MORPH object}
prepareMorphObjectFromFiles = function(InputConfig= NULL,InputGOI = NULL) {
  Config = read.delim(InputConfig, sep = "\t", header=FALSE) #Reads the configs.txt file.
  List_GE = as.character(Config[,1]) #Reads the first column (containing paths to gene-expression data files)
  List_C = as.character(Config[,2]) #Reads the second column (containing paths to clustering solution files)
  G = c() #Initialize the vector to contain names of pathway-genes.
  if (! is.null(InputGOI)){ #Is the path to the pathway-genes file defined?
    G = GetGOIs(InputGOI) #Reads file into variable.
  }
  ge_datasets = list() #Initialize list to contain all gene-expression matrices.
  clustering_solutions = list() #Initialize list to contain all clustering solutions.
  #Go over all the datasets and solutions pairs, and insert them to the data structures.
  for (I in (1:length(List_C))) {
    cl = List_C[[I]]
    ge = List_GE[[I]]
    #Read the dataset only if it still doesn't exists in ge_datasets.
    #As each dataset can appear many time in the config file, we insert it only once to the data
    #structure (preformance reasons)
    if (is.null(ge_datasets[[ge]])) {
      ge_data = getGeneExpression(ge) #Read data into variable.
      ge_datasets[[ge]] = ge_data #Add gene-expression dataset to the list.
    }
    #Insert the current solution to the clustering solutions data stucture.
    #The processing of the solution is against the current dataset.
    #Clustering solution should be unique to each dataset (recommended).
    raw_cl =  getClusteringInformation(cl) #Read data into variable.
    genes_to_keep = intersect(rownames(raw_cl),rownames(ge_data)) #Keep only the genes that appear both in the gene-expression data AND the clustering solution.
    raw_cl = as.matrix(raw_cl[genes_to_keep,],ncol=1) #Extract the kept-genes' lines as a matrix.
    rownames(raw_cl) = genes_to_keep #Rename matrix's rows accordingly.
    raw_cl_genes = rownames(raw_cl) #Extract names genes in clustering solution.
    raw_cl_clusters = raw_cl[,1] #Extract parent clusters of genes in clustering solution.
    unclustered = setdiff(rownames(ge_data),raw_cl_genes) #Detect unclustered genes (those who appear in gene-expression but NOT in clustering solution).
    if (length(unclustered) == 0){ #Have we met NO unclustered genes?
      clustering_solutions[[cl]] = raw_cl #Add clustering solution data to the list.
      next #Move to next clustering solution.
    }
    #Have we met unclustered genes?
    uncl = rep ("unclustered",length(unclustered)) #Set their parent clusters as "unclustered".
    clustering = matrix(c(raw_cl_clusters,uncl),ncol=1) #Create a matrix with the unclustered genes at the bottom.
    rownames(clustering) = c(raw_cl_genes,unclustered) #Rename the rows according to gene names.
    clustering_solutions[[cl]] = clustering #Add clustering solution data to the list.
  }									                          
  #Create and return a MORPH object (see function description).
  morph_obj = list(config = Config, clustering_solutions = clustering_solutions, ge_datasets = ge_datasets, pathway_genes = G)
  class(morph_obj)<-"morph_data" # equvallent to Java class
  return (morph_obj)}

### [2.2] rankGenes {Input;G-pathway, C-cluster solution, GE-dataset, Output; Scores}
rankGenes = function(G,C,GE,corrs_mat = NULL) #Ranks all candidates-vs-GOIs for a "data-couple" (GE data & Clustering solution)
{if (is.null(corrs_mat)){ #Have we NOT received a coexpression (covariance) matrix (is it defined)?
  corrs_mat = cor(t(GE),t(GE[G,])) #Calculate the covariance matrix from scratch.
  print(corrs_mat)
}
  print(corrs_mat)
  G_Clusters = matrix(C[G,1], length(G),1) #Detect the parent cluster for all pathway-genes.
  rownames(G_Clusters) = G #Change the row names according to the pathway-genes.
  RelevantClusters =na.omit(unique(G_Clusters[,1]))  #Keep only the clusters that contain pathway-genes.
  Scores = list() #Initialize list to contain results.
  Scores[["Scored"]] = c() #Initialize the "Scored" field that will contain all valid calculated scores.
  for (Cluster in RelevantClusters) #Go through each cluster within those who contain pathway-genes.
  { #Acquire cluster-specific ("current") pathway-genes and candidates.
    CurrentG = na.omit(names(G_Clusters[(G_Clusters[,1]==Cluster),]))  #Acquire current pathway-genes
    CurrentB = rownames(C)[which(C[,1]==Cluster)] #Acquire current background genes.
    Candidates = setdiff(CurrentB,CurrentG) #setdiff function filters Candidates are genes which are present at the background but NOT in our pathway-genes group.
    
    #Create an ordered vector containing the normalized covariance scores for all candidates vs. pathway-genes.
    NormalizedCorrelations = getNormalizedCorrelations(corrs_mat, CurrentG, Candidates) #Calculate normalized scores.
    NormalizedCorrelations = NormalizedCorrelations[order(NormalizedCorrelations,decreasing=TRUE)] #Sort list in descending order.
    Scores[["Scored"]] = c(Scores[["Scored"]],NormalizedCorrelations) #Add normalized scores to the "Scored" field of our list.
  }
  # All genes that were not clustered with any pathway gene are placed at the bottom of the ranking (normalized score = NA)
  NonRelevantClusters = setdiff(unique(C[,1]),RelevantClusters) #Detect clusters that do NOT contain pathway-genes (non-relevant clusters).
  NonCandidates = character(0) #Initialize variable to contain names of genes who are not clustered with pathway-genes (non-candidates).
  for (NRC in NonRelevantClusters){ #Go through all non-relevant clusters.
    NC = rownames(C)[which(C[,1]==NRC)] #Extract names of all genes in current cluster.
    NonCandidates = c(NonCandidates,NC) #Add the names to the non-candidates list.
  }
  LastScores = rep(NA,length(NonCandidates)) #Set scores of all non-candidates as NA.
  names(LastScores) = NonCandidates #Rename vector according to non-candidates' names.
  Scores[["Rejected"]] = LastScores #Add non-candidates to the "Rejected" field of our list.
  #Sort all scores in a descending order.
  scores_vector = Scores[["Scored"]]
  scores_vector = scores_vector[order(scores_vector,decreasing=T)]
  Scores[["Scored"]] = scores_vector
  return(Scores)}


### [2.3] Normalized correlation {Input; correlation matrix, pathway, genes, Outputl;Ordered vector}

#Acquire cluster-specific ("current") pathway-genes and candidates.
getNormalizedCorrelations = function(corrs_mat, CurrentG, Candidates)
{AverageCorrelations = c() #Initialize variable to contain average correlations.
if (length(CurrentG)==1){ #Is there only ONE pathway gene in this cluster?
  AverageCorrelations = corrs_mat[Candidates,CurrentG] #Set average correlation as it's coexpression (covariance) score from our matrix.
}
else{ #Are there MORE than one pathway genes in this cluster?
  AverageCorrelations = apply(corrs_mat[Candidates,CurrentG],1,mean) #Calculate the mean of their coexpression (covariance) scores.
}
#Normalize and return the scores.
names(AverageCorrelations)<-Candidates #Rename the average correlations according to the appropriate candidates' names.
Mean = mean(AverageCorrelations) #Calculate the mean of all average scores.
SD = sd(AverageCorrelations) #Calculate the standard deviation of all average scores.
NormalizedCorrelations = ((AverageCorrelations - Mean) / SD) #Normalize scores according to mean and standard deviation.
return(NormalizedCorrelations)}


### [2.4] LOOCV (Leave-One-Out Cross Validation) {Input; G, C, GE, Output;Score}
LOOCV = function(G,C,GE,K,corrs_mat = NULL, NameGE = NULL, NameC = NULL)
{if (is.null(corrs_mat)){ #Have we NOT received a coexpression (covariance) matrix (is it defined)?
  corrs_mat = cor(t(GE),t(GE[G,])) #Calculate the covariance matrix from scratch.
}
  SelfRanks = numeric(0) #Initialize vector to contain self-ranks.
  for (V in G) #Go through all pathway-genes and leave one out at a time.
  { S = setdiff(G,V) #Reduce pathway genes group to contain all pathway genes EXCEPT the one left out.
  Scores = rankGenes(S,C,GE,corrs_mat) #Run rankGenes function with the reduced pathway-genes.
  Rank = which(names(Scores$Scored)==V) #Acquire the location of the left-out gene among the others.
  if (length(Rank)==0 || is.nan(Rank)){Rank = (2*K)} #Was the left-out genes NOT detected at all? Set score to 2K as penalty.
  SelfRanks=c(SelfRanks,Rank) #Add self-rank to our vector.
  }
  return(AUSR(SelfRanks,K,NameGE,NameC))}

### [2.5] AUSR (Area Under the curve of the Self-Ranked genes) {Input; Self-ranks, K-threshold, Output;Score}
AUSR = function(SelfRanks, K, NameGE = NULL, NameC = NULL)
{nGOIs = length(SelfRanks) #Count the number of self-ranks in the vector.
Fractions = numeric(0) #Initialize vector to contain fractions (i.e. fractions of pathway genes with self rank below threshold)
AUC = 0 #Initialize accumulative variable to contain Area-Under-the-Curve score.
PlotFlag = 0 #Initialize flag to indicate whether plotting should be executed.
if ((!is.null(NameGE)) & (!is.null(NameC)))	{ #Have we received both the names of gene-expression data file AND the clustering solution file?
  PlotFlag = 1 #User requested for plotting, set flag to 1.
}
for (Index in 1:K) { #Go through threshold calculation-range (1 to K by 1)
  Fraction = ((length(SelfRanks[SelfRanks<=Index])) / nGOIs) #Calculate the fraction of genes that have a self-rank smaller or equal to K of current iteration.
  Fractions = c(Fractions,Fraction) #Add the fraction score to our vector.
  AUC = AUC + Fraction #Accumulate fraction score as AUC score.
}
#Calculate and return the AUSR score (i.e. AUC divided by threshold value).
return((AUC / K))}


### [2.6] MORPH {Input; MORPH object, k, output;Ranks, Score, AUSR, C, GE, }
MORPH = function(morph_obj, K=1000, view = FALSE)
{#Acquire data from input MORPH object
  G = morph_obj$pathway_genes #Pathway-genes.
  List_C = as.character(morph_obj$config[,2]) #Paths to clustering solution files.
  List_GE = as.character(morph_obj$config[,1]) #Paths to gene-expression data files.
  ge_datasets = morph_obj$ge_datasets #Gene-expression matrices.
  cl_solutions = morph_obj$clustering_solutions #Clustering solutions.
  #Initialize global variables to contain - 
  Gene_Scores = list() #Normalized gene scores.
  AUSR_Scores = numeric(0) #AUSR scores.
  Final_Scores = list() #MORPH Results object to contain all approriate information (see function description for more details).
  #Calculate the coexpression (covariance) matrices for all gene-expression matrices.
  corr_matrices = list() #Initialize list to contain all covariance matrices.
  ges = unique(List_GE) #Eliminate duplicates from paths to gene-expression data files.
  for (ge in ges){ #Go through all gene-expression matrices
    GE = ge_datasets[[ge]] #Extract gene-expression matrix.
    P = intersect(G,rownames(GE)) #Acquire relevant genes (those which appear both in the pathway-genes AND the gene-expression data).
    corr_matrices[[ge]] = cor(t(GE),t(GE[P,])) #Calculate and save coexpression (covariance) matrix into out list.
  }
  for (I in (1:length(List_C))){ #Go through all configurations (clustering solution vs. gene-expression dataset)
    #Acquire necessary information
    C = cl_solutions[[(List_C[I])]] #Clustering solution.
    curr_corrs = corr_matrices[[List_GE[I]]] #Coexpression (covariance) matrix.
    Pathway = removeAbscentGOIs(G,C,rownames(curr_corrs)) #Eliminate all pathway genes that do not appear in this configuration.
    
    Ranking = rankGenes(Pathway,C,GE=NULL,curr_corrs) #Rank candidates vs. pathway-genes.
    AUC = numeric(0) #Initialize variable to contain AUSR score.
    
    if (view == 0){ #Is the "view" parameter set to FALSE?
      AUC = LOOCV(Pathway,C,GE=NULL,K,curr_corrs) #Proceed to execute LOOCV as usual.
    } else { #Is the "view" parameter set to TRUE?
      AUC = LOOCV(Pathway,C,GE=NULL,K,curr_corrs, List_GE[I],List_C[I]) #Pass names of configuration files for plotting purposes.
    }
    #Add approriate information to the list (see function description for more details).
    Score = list(C=List_C[I],GE=List_GE[I],AUSR = AUC, Ranking = Ranking) #Arrange information into a single data block.
    Final_Scores[[I]] = Score #Add data block to the list.
  }
  if (view == TRUE){ #Is the "view" parameter set to TRUE?
    dev.off() #Turn off console presentation in order to write into output file and not to screen.
  }
  #Set final scores class as "morph_results", and return it.
  class(Final_Scores)<-"morph_results"
  return(Final_Scores)}

### [2.7] remove abscent GOIs {Input; G, C, Output; Pathway genes appeared in both}
removeAbscentGOIs = function(G,C,GENames)
{int = intersect(G,rownames(C)) #Keep only pathway-genes that appear in clustering solution.
int = intersect(int,GENames) #Keep only pathway-genes that appear in gene-expression data.
return (int)}

### [2.8] LOOCV_MORPH {Input; Morph object, k, output; AUSR score}
LOOCV_MORPH = function(morph_obj,K=1000)
{SelfRanks = numeric(0) #Initialize variable to contain all self-ranks.
for (V in morph_obj$pathway_genes) #Go through all pathway-genes and leave one out at a time.
{# Create a new MORPH input object but EXCLUDE the tested gene.
  new_input = list(config = morph_obj$config,
                   clustering_solutions = morph_obj$clustering_solutions,
                   ge_datasets = morph_obj$ge_datasets, 
                   pathway_genes = setdiff(morph_obj$pathway_genes,V))
  class(new_input)<-"morph_data"
  morph_run_results = MORPH(morph_obj=new_input,K) #Run the MORPH algorithm on the current setting.
  Scores = getMorphPredictions(morph_run_results) #Get the normalized correlation scores of current run.
  Rank = which(names(Scores)==V) #Acquire the location of the left-out gene among the others.
  if (length(Rank)==0 || is.nan(Rank)){Rank = (2*K)} #Was the left-out genes NOT detected at all? Set score to 2K as penalty.
  SelfRanks=c(SelfRanks,Rank) #Add self-rank to our vector.
}
#Calculate and return AUSR score for the MORPH process.
return(AUSR(SelfRanks,K))}

### [2.9] getMorphResultBestConfig {Input; MORPH-result object, Output; Ranks, c, GE, AUSR, Score}
getMorphResultBestConfig = function(morph_res_obj){ # list of lists [like PERL hush table]
  #Initialize variables to contain - 
  best_ausr = -1 #Best AUSR score.
  best_ind = -1 #Index of best-AUSR configuration.
  clustering_solution = "" #Name of clustering solution file in best configuration.
  gene_expression = "" #Name of gene-expression data file in best configuration.
  for (i in 1:length(morph_res_obj)){ #Go through all configurations in MORPH Results object.
    curr_res = morph_res_obj[[i]] #Extract current results data block.
    ausr = curr_res$AUSR #Acquire AUSR score.
    if (ausr>best_ausr){ #Is the AUSR score better than the last recorded best?
      best_ind = i; best_ausr = ausr #Save both the score and the configuration index.
      clustering_solution = morph_res_obj[[i]]$C #Acquire the clustering solution file name.
      gene_expression = morph_res_obj[[i]]$GE #Acquire the gene-expression data file name.
    }
  }
  # Create and return a list containing information regarding best configuration (for more information see function description).
  morph_best_res = list(AUSR = best_ausr, Ranking = morph_res_obj[[best_ind]]$Ranking, C = clustering_solution, GE = gene_expression)
  return (morph_best_res)}

### [2.10] get Morph Predictions {Input; Morph-result object, output;Ranking score}
getMorphPredictions = function(morph_res_obj){
  return (getMorphResultBestConfig(morph_res_obj)$Ranking$Scored) #Acquire and return scores of best configuration.}
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
###############################################################################################################
# Set working directory -- needs to be personalized
setwd("C:/Users/elise/Documents/Schuul/Stage/Project/Data")

# Define the path to your Excel file
file_path <- "PathwaysMtb.xlsx"

# Get the names of all sheets in the Excel file
sheet_names <- excel_sheets(file_path)

# List pathways
listPathways <- c()

# Loop through each sheet and perform operations
for (sheet in sheet_names) {
  # Add name sheet to list
  listPathways <- c(listPathways, sheet)
  # Read the data from the current sheet
  pathway <- read_excel(file_path, sheet = sheet, col_names=FALSE)
  # Get data 
  pathway <- pathway$...1
  # Print the sheet name
  print(paste("Sheet:", sheet))
  # Print the first few rows of the data 
  print(head(pathway))
  # Write files for each pathway
  writeLines(pathway, paste0(sheet, ".txt"), sep="\t")
}
print(listPathways)

#Solutions
dfpathway <- c()
dfscore <- c()
for (pathwayName in listPathways) {
  ## MORPH input
  InputConfig = "Configs.txt"
  InputGOI = paste0(pathwayName, ".txt")
  #InputGOI = "Pathway1.txt"
  morph_input = prepareMorphObjectFromFiles(InputConfig, InputGOI)
  Scores = MORPH(morph_input)
  ## Removing Absent genes
  G = morph_input$pathway_genes #Get pathway genes
  C = (morph_input$clustering_solution)[[1]] #Get clustering solution
  GE = (morph_input$ge_datasets)[[1]] #Get gene expression dataset
  GENames = rownames(GE) #Get names of genes in dataset
  Intersection = removeAbscentGOIs(G,C,GENames)
  print(Intersection)
  ## AUSR
  BestConfig <- getMorphResultBestConfig(Scores)
  AUSRscore <- as.character(BestConfig$AUSR)
  print(BestConfig$AUSR) 
  ## MORPH gene scores
  Predictions <- getMorphPredictions(Scores)
  Candidates <- as.matrix(head(Predictions))
  CandidatesGenes <- rownames(Candidates)
  CandidatesScores <- as.character(Candidates)
  Candidates <- cbind(CandidatesGenes, CandidatesScores)
  colnames(Candidates) <- NULL
  print(head(Predictions))
  # Write file with AUSR score and candidate genes
  writeLines("AUSR score:", paste0(pathwayName, "_MORPHMtb.txt")) #Create file
  fileConn <- file(paste0(pathwayName, "_MORPHMtb.txt"),"a")#Open connection to append
  writeLines(AUSRscore, fileConn) #Append
  writeLines("\n", fileConn)
  writeLines("Top candidate genes:", fileConn)
  write.matrix(Candidates, fileConn, sep="\t")
  close(fileConn) # Close connection
  dfpathway <- c(dfpathway, pathwayName)
  dfscore <- c(dfscore, AUSRscore)
}

# Make dataframe from pathways with their AUSR scores
dfBestAUSR <- data.frame(Pathway = dfpathway, Score = dfscore)
# Make variable with pathway names from 3 pathways with the highest AUSR score
BestAUSRpathways <- head(dfBestAUSR[order(dfBestAUSR$Score, decreasing = TRUE), 1],3)
# Print name of top 3 AUSR-scored pathways
print(BestAUSRpathways)

#Get scores and descriptions of top 3 AUSR-scored pathways
bestcand <- list()
bestdescr <- list()
bestscore <- list()
for (bestPathway in BestAUSRpathways) {
  ## MORPH input
  InputConfig = "Configs.txt"
  InputGOI = paste0(bestPathway, ".txt")
  morph_input = prepareMorphObjectFromFiles(InputConfig, InputGOI)
  Scores = MORPH(morph_input)
  ## MORPH gene scores
  Predictions <- getMorphPredictions(Scores)
  Candidates <- as.matrix(head(Predictions))
  CandidatesGenes <- rownames(Candidates)
  CandidatesScores <- round(Candidates,6)
  descr <- Description(CandidatesGenes)
  bestcand <- append(bestcand, list(CandidatesGenes))
  bestdescr <- append(bestdescr, list(descr))
  bestscore <- append(bestscore, list(CandidatesScores))
}

#length <- length(bestcand[[1]])

# zip candidate genes and their descriptions together
zipped <- Map(function(x,y,z,w) list(x,y,z,w), bestcand, bestdescr, bestscore, BestAUSRpathways)

#make data frames from genes with their description
bestpath <- list()
for (i in seq_along(zipped)) {
  bestpath[[i]] <- data.frame(Candidate = zipped[[i]][[1]], 
                              Description = zipped[[i]][[2]], 
                              Score = zipped[[i]][[3]],
                              Pathway = zipped[[i]][[4]])
  # get number of rows in the data frame
  num_rows <- nrow(bestpath[[i]])
  # set the row names to numbers
  rownames(bestpath[[i]])<- seq_len(num_rows)
}
# Print bestpath
print(bestpath)


#get the candidate genes with description "hypothetical protein" 
hypotheticalProt <- list()
for (i in seq_along(bestpath)) {
      hypotheticalProt[[i]] <- subset(bestpath[[i]], Description == "hypothetical protein")
}
# Add pathwayname to dataframe with "hypothetical protein" candidate genes
for (i in seq_along(BestAUSRpathways)) {
  hypotheticalProt[[i]]$Pathway <- BestAUSRpathways[i]
}

# Print hypotheticalProt
print(hypotheticalProt)


# get file with candidate genes for top 3 AUSR-scored pathways 
for (i in seq_along(bestpath)) {
  # make data frame matrix
  bestpathfile <- as.matrix(bestpath[[i]][, c(1,3,2)])
  #open connection to append
  fileConn <- file("topAUSRscoredPathway.txt","a")
  writeLines("\n", fileConn)
  writeLines(paste0(i,"."),fileConn)#Append
  writeLines(BestAUSRpathways[i], fileConn)
  write.table(bestpathfile, fileConn, sep="\t\t", quote=FALSE, col.names=FALSE)
  writeLines("\n", fileConn)
  writeLines("Hypothetical proteins: ", fileConn)
  write.table(hypotheticalProt[[i]][, c(1,3,2)], fileConn, sep="\t\t", quote=FALSE, col.names=FALSE)
  writeLines("---------------------------------------------------------------------------------------------------------------", fileConn)
  # Close connection
  close(fileConn) 
}

###############################################################################################################
##FARIDA##

## Use gene ontology to get gene IDs and calculate pathways instead of giving the input pathways
## Maybe first use an if-statement to see what is given as input; excel file or gene ontology 
## Get gene IDs
## Use MORPH algorithm to get AUSR scores and candidate genes
## get top ... pathways with best AUSR scores
## get their candidate genes
## highlight the genes with description "hypothetical protein"
## you can use the code above as a guideline; once you have the gene IDs, you can basically use the same code

