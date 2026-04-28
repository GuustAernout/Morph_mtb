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
  "RMySQL"
)
lapply(required_packages, library, character.only = TRUE)


data <-read.csv(file="GSE292408_genes_x_samples.txt", 
                sep='\t', 
                header = TRUE
                )
rownames(data) <- data[,1]
data[,1] <- NULL
  a <- DGEList(data, group = NULL) # create DGEList object
  b <- calcNormFactors(a, method = "TMM") # perform TMM normalization
  norm_data <- cpm(b, log=FALSE) # retrieve normalized counts
  sd_expr <- apply(norm_data, 1, sd)   # SD for each gene
  threshold <- 1    # threshold
  norm_data_filtered <- norm_data[sd_expr >= threshold,]  # remove gene with sd<1




log_data <- round(log2(norm_data_filtered + 1),2)


geneids <- data
  gene_ids <- rownames(data)  # ← Use row names instead
  a <- DGEList(data, group = NULL)
  b <- calcNormFactors(a, method = "TMM")
  norm_data <- cpm(b, log=FALSE)
  sd_expr <- apply(norm_data, 1, sd)
  threshold <- 1
  geneids <- gene_ids[sd_expr >= threshold]



########################################################################################################################
## Clustering

kmeans_data <-  kmeans(log_data, 3, nstart=50,iter.max = 40)

kmeans_file_data <- data.frame(geneids,kmeans_data$cluster)
output_file_path <- "kmeansProtonPump.txt"
write.table(kmeans_file_data, file = output_file_path, sep = "\t", row.names = FALSE,col.names = FALSE, quote = FALSE)
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
  n_genes <- nrow(data)
  xdim <- floor(n_genes / 10)  # floor ensures nodes <= genes so som() can initialise
  som_grid <- somgrid(xdim = xdim, ydim = 10, topo = "hexagonal")
  gene_som <- som(data, grid = som_grid, rlen = 100, alpha = c(0.05, 0.01), keep.data = TRUE)
  weight <- getCodes(gene_som)
  # Return both the weight matrix (for kmeans on nodes) and unit.classif
  # (which SOM node each gene belongs to, needed to map node-clusters back to genes)
  list(weights = weight, unit_classif = gene_som$unit.classif)
}
wsssom <- function(data){
  n_genes <- nrow(data)
  xdim <- floor(n_genes / 10)  # floor ensures nodes <= genes
  som_grid <- somgrid(xdim = xdim, ydim = 10, topo = "hexagonal")
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