# Functions

## Normalization and Filtering after normalization
Filtering <- function(data){
  a <- DGEList(data, group = NULL) # create DGEList object
  b <- calcNormFactors(a, method = "TMM") # perform TMM normalization
  norm_data <- cpm(b, log=FALSE) # retrieve normalized counts
  sd_expr <- apply(norm_data, 1, sd)   # SD for each gene
  threshold <- 1    # threshold
  norm_data_filtered <- norm_data[sd_expr >= threshold, ] # remove gene with sd<1
  return(norm_data_filtered)
}

log <- function(norm_data_filtered){
  log_data <- round(log2(norm_data_filtered + 1),2)
  return(log_data)
}

geneids <- function(data){
  gene_ids <- data$Gene_ID # Extract gene identifiers
  a <- DGEList(data, group = NULL) # create DGEList object
  b <- calcNormFactors(a, method = "TMM") # perform TMM normalization
  norm_data <- cpm(b, log=FALSE) # retrieve normalized counts
  sd_expr <- apply(norm_data, 1, sd)   # SD for each gene
  threshold <- 1    # threshold
  geneids <- gene_ids[sd_expr >= threshold]
  return(geneids) 
}

########################################################################################################################
## Clustering
wsskmeans <- function(logdata){
  # Elbow Method
  ## K-means
  kmax <- 10
  wssk <- sapply(2:kmax, function(k){
    kmeans(logdata, k, nstart=50,iter.max = 40)$tot.withinss})
  wssk
}
clusterK <- function(wssk, elbowK){
  kmax <- 10
  plot(2:kmax, wssk, type="b", pch =10, frame = FALSE, 
       xlab="Number of clusters",
       ylab="Total Within sum of square",
       main="K-means",
       cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
  abline(v = elbowK, col = "red", lwd = 2)
}

weightSOM <- function(data){
  # Elbow Method
  ## SOM
  som_grid <- somgrid(xdim = 376, ydim = 10, topo = "hexagonal")   # depend with the size of the dataset (3760)
  gene_som <- som(data, grid = som_grid, rlen = 100, alpha = c(0.05, 0.01), keep.data = TRUE)
  weight<- getCodes(gene_som)
  weight}
wsssom <- function(data){
  som_grid <- somgrid(xdim = 376, ydim = 10, topo = "hexagonal")   # depend with the size of the dataset (3760)
  gene_som <- som(data, grid = som_grid, rlen = 100, alpha = c(0.05, 0.01), keep.data = TRUE)
  weight<- getCodes(gene_som)
  kmax <- 10
  wsss <- sapply(2:kmax, function(k){
    kmeans(weight, k, nstart=50,iter.max = 40)$tot.withinss})
  wsss
}
clusterS <- function(wsss, elbowS){
  kmax <- 10
  plot(2:kmax, wsss, type="b", pch =10, frame = FALSE, 
       xlab="Number of clusters",
       ylab="Total Within sum of square",
       main="SOM",
       cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
  abline(v = elbowS, col = "red", lwd = 2) 
}

## Define clusters 
## k-means
kmc <- function(data, elbowK){
  kmeans(data, centers = elbowK, iter.max=40,nstart=50)
}
## SOM
SOMc <- function(weight, elbowS){
  kmeans(weight, centers = elbowS, iter.max=40,nstart=50)
} 

##########################################################################################################
## Data Preparation for MORPH ##
# Preparing 

## K-means clusters
kmeansclusters <- function(gene_ids, kmc){
  genes <- gene_ids
  cluster <- kmc$cluster
  kmeans_cluster <- data.frame(genes, cluster)
  return(kmeans_cluster)
}


## SOM clusters
somclustersss <- function(gene_ids, somc){
  gene_names <- gene_ids
  cluster_assignments <- somc$cluster
  gene_clusters <- data.frame(genes = gene_names, cluster = cluster_assignments)
  return(gene_clusters)
}


########################################################################################################################
## MORPH ALGORITHM ##
# Functions
## Reading data functions
### [1.1] getClusteringInformation
getClusteringInformation <- function(ClusterFile)
{
  Cluster = read.delim(ClusterFile, sep="\t", header=FALSE, row.names = 1,stringsAsFactors=F) #Read tab-delimited file.
  colnames(Cluster) = ClusterFile #Change cluster column name to file-name for future identification.
  return(Cluster)
}

### [1.2] Get Pathways
GetGOIs <- function(InputGOI)
{
  GOI = as.character(strsplit(InputGOI, "\t")) #Splits names into a vector by tabs.
  return(GOI)
}


### [1.3] getGeneExpression
getGeneExpression <- function(InputGE)
{
  GE = (read.delim(InputGE, sep="\t", header=TRUE,row.names=1)) #Reads the gene-expression data from the tab-delimited file.
}


## MORPH Algorithm
### [2.1] prepare Morph Object From Files {Input;Configuration file (data and cluster solution), Output:MORPH object}
prepareMorphObjectFromFiles <- function(InputGOI = NULL, ...) {
  #Config = read.delim(InputConfig, sep = "\t", header=FALSE) #Reads the configs.txt file.
  Config = data.frame(V1=c("methylation_dataset.txt","methylation_dataset.txt","drug_cholesterol_toxin.txt","drug_cholesterol_toxin.txt","ESX_WT_mutant.txt","ESX_WT_mutant.txt",
                           "MTBc_all_lineages.txt","MTBc_all_lineages.txt","primary_drug.txt","primary_drug.txt","timecourse_nitricacid.txt","timecourse_nitricacid.txt"),
                      V2=c("kmeansmethylation.txt","sommethylation.txt","kmeansdrug.txt","somdrug.txt","kmeansESX.txt","somESX.txt",
                           "kmeansMTBc_lineages.txt","somMTBc_lineages.txt","kmeansprimary.txt","somprimary.txt","kmeanstimecourse.txt","somtimecourse.txt"))
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
  class(morph_obj)<-"morph_data" # equivalent to Java class
  return (morph_obj)}

prepareMorphObjectFromFiles2 <- function(InputGOI = NULL, ...) {
  #Config = read.delim(InputConfig, sep = "\t", header=FALSE) #Reads the configs.txt file.
  Config = data.frame(V1=c("methylation_dataset.txt","methylation_dataset.txt","drug_cholesterol_toxin.txt","drug_cholesterol_toxin.txt","ESX_WT_mutant.txt","ESX_WT_mutant.txt",
                           "MTBc_all_lineages.txt","MTBc_all_lineages.txt","primary_drug.txt","primary_drug.txt","timecourse_nitricacid.txt","timecourse_nitricacid.txt", "ExpressionData.txt", "ExpressionData.txt"),
                      V2=c("kmeansmethylation.txt","sommethylation.txt","kmeansdrug.txt","somdrug.txt","kmeansESX.txt","somESX.txt",
                           "kmeansMTBc_lineages.txt","somMTBc_lineages.txt","kmeansprimary.txt","somprimary.txt","kmeanstimecourse.txt","somtimecourse.txt", "kmeansexpdata.txt", "somexpdata.txt"))
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
  class(morph_obj)<-"morph_data" # equivalent to Java class
  return (morph_obj)}


### [2.2] rankGenes {Input;G-pathway, C-cluster solution, GE-dataset, Output; Scores}
rankGenes <- function(G,C,GE,corrs_mat = NULL) #Ranks all candidates-vs-GOIs for a "data-couple" (GE data & Clustering solution)
{if (is.null(corrs_mat)){ #Have we NOT received a coexpression (covariance) matrix (is it defined)?
  corrs_mat = cor(t(GE),t(GE[G,])) #Calculate the covariance matrix from scratch.
  print(corrs_mat)
}
  print(corrs_mat)
  G_Clusters = matrix(C[G,1], length(G),1) #Detect the parent cluster for all pathway-genes.
  rownames(G_Clusters) = G #Change the row names according to the pathway-genes.
  RelevantClusters = na.omit(unique(G_Clusters[,1]))  #Keep only the clusters that contain pathway-genes.
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
getNormalizedCorrelations <- function(corrs_mat, CurrentG, Candidates)
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
LOOCV <- function(G,C,GE,K,corrs_mat = NULL, NameGE = NULL, NameC = NULL)
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
AUSR <- function(SelfRanks, K, NameGE = NULL, NameC = NULL)
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
if (PlotFlag == 1) { #Is the plotting flag "up"?
  #Set X-axis to calculation-range and plot the data.
  X = (1:K)
  plot(X,Fractions, type="l", col="blue",xlab="Self-rank threshold", ylab="Fraction of pathway genes with self rank below threshold",main=paste(NameGE,"with",NameC))
}
#Calculate and return the AUSR score (i.e. AUC divided by threshold value).
return((AUC / K))}

### [2.6] MORPH {Input; MORPH object, k, output;Ranks, Score, AUSR, C, GE, }
MORPH <- function(morph_obj, view=FALSE)
{#Acquire data from input MORPH object
  K = 1000
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
  if (view == TRUE){ #Is the "view" parameter set to TRUE?
    pdf("Threshold_Plots.pdf", onefile=TRUE) #Open PDF formatted file to hold all self-ranking plots.
  }
  for (I in (1:length(List_C))){ #Go through all configurations (clustering solution vs. gene-expression dataset)
    #Acquire necessary information
    C = cl_solutions[[(List_C[I])]] #Clustering solution.
    curr_corrs = corr_matrices[[List_GE[I]]] #Coexpression (covariance) matrix.
    rnames = rownames(curr_corrs)
    Pathway = removeAbscentGOIs(G,C,rnames) #Eliminate all pathway genes that do not appear in this configuration.
    
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
  return(Final_Scores)
  }

### [2.7] remove abscent GOIs {Input; G, C, Output; Pathway genes appeared in both}
removeAbscentGOIs <- function(G,C,GENames)
{int = intersect(G,rownames(C)) #Keep only pathway-genes that appear in clustering solution.
int = intersect(int,GENames) #Keep only pathway-genes that appear in gene-expression data.
return (int)}

### [2.8] LOOCV_MORPH {Input; Morph object, k, output; AUSR score}
LOOCV_MORPH <- function(morph_obj,K=1000)
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
  
  ### [2.11] get Scores Distribution Plots  {Input:Scores,path for saving,  output;plot}
  #getScoresDistributionPlots = function(Scores, OutputFile="ScoresDistribution.pdf", Color="blue", Type = "b")
  #{pdf(OutputFile, onefile=TRUE) #Open the output file for writing in PDF format.
  #  for (i in 1:length(Scores))	{ #Go through all MORPH Results object data blocks
  #    Ranks = Scores[[i]] #Acquire current data block.
  #    Y = as.numeric((Ranks$Ranking)$Scored) #Extract normalized correlation scores into a numeric vector.
  #    X = 1:(length(Y)) #Create an X-axis in the appropriate length.
  #    plot(X,Y,col=Color,type=Type, xlim=c(0,1), xlab="Self-rank threshold", ylab="Fraction of pathway genes with self rank below threshold",main=paste("Clustering:",Ranks$C,"; Gene expression:",Ranks$GE)) #Plot data.
  #  }
  #  dev.off() #Turn off console presentation in order to write into output file and not to screen.
  #}
  }

###################################################################################################################################################################################################################"

## Sampling for Real pathway recognition
# Get random genes as random pathway with length same as input pathway
uniform_sample <- function(vector, lengthPatway) {
  index <- sample(length(vector), lengthPatway) 
  return(vector[index])}

getScoresRandomPathway <- function(g, random, lengthPathway){
  # Make empty vector to store Scores in
  Score <- c()
  # Variable B which stores the amount of random pathways needed to be calculated
  B <- random
  lapply(1:B, function(i){
    InputGOI <- unlist(uniform_sample(g, lengthPathway))
    morph_input <- prepareMorphObjectFromFiles(InputGOI)
    print(names(morph_input)) 
    G <- morph_input$pathway_genes
    NameC <- names(morph_input$clustering_solutions)[1] 
    C <- (morph_input$clustering_solutions)[[NameC]] #Get first clustering solution
    NameGE <- names(morph_input$ge_datasets)[1] 
    GE <- (morph_input$ge_datasets)[[NameGE]] #Get first gene expression dataset
    morph_input <- prepareMorphObjectFromFiles(InputGOI)
    Scores <- MORPH(morph_input, view = TRUE)
    # Append Scores to vector Score
    Score <- c(Score, Scores)
    return(Score)
  })
}

