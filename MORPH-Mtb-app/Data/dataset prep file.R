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
output_file_path <- "kmeansProtonPump.txt" # adjust file name to suit your dataset

write.table(kmeans_file_data, file = output_file_path, sep = "\t", row.names = FALSE,col.names = FALSE, quote = FALSE)




## SOM clustering
# Use log_data (filtered + log-transformed), same as the app
n_genes <- nrow(log_data)
xdim <- floor(n_genes / 10)  # floor ensures SOM nodes <= genes so som() can initialise
som_grid <- somgrid(xdim = xdim, ydim = 10, topo = "hexagonal")
gene_som <- som(log_data, grid = som_grid, rlen = 100, alpha = c(0.05, 0.01), keep.data = TRUE)

# Weight matrix: one row per SOM node
weight <- getCodes(gene_som)

# unit.classif: one entry per gene, recording which SOM node it was assigned to
unit_classif <- gene_som$unit.classif

# Choose number of SOM clusters (adjust elbowS to match your dataset)
elbowS <- 3

# Run kmeans on SOM node weights to get node-level clusters
somc <- kmeans(weight, centers = elbowS, iter.max = 40, nstart = 50)

# Map each gene to its node's cluster via unit.classif
gene_cluster_assignments <- somc$cluster[unit_classif]

# Build output: gene ID + cluster number, one row per gene
som_file_data <- data.frame(geneids, gene_cluster_assignments)
output_file_path_som <- "somProtonPump.txt"  # adjust file name to suit your dataset

write.table(som_file_data, file = output_file_path_som, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
