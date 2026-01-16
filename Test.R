# BiocManager::install('SingleCellExperiment')
# BiocManager::install(c('scuttle', 'scran', 'scater', 'uwot', 'rtracklayer'))
library(SingleCellExperiment)
library(scuttle)


mat <- read.delim(file.path("~/Téléchargements/E-MTAB-5522/counts_Calero_20160113.tsv"),
                  header=TRUE, row.names=1, check.names=FALSE)

# Only considering endogenous genes for now.
spike.mat <- mat[grepl("^ERCC-", rownames(mat)),] 
mat <- mat[grepl("^ENSMUSG", rownames(mat)),] 

# Splitting off the gene length column.
gene.length <- mat[,1]
mat <- as.matrix(mat[,-1]) 

dim(mat)

sce <- SingleCellExperiment(assays = list(counts = mat))
sce
mat2 <- counts(sce)

sce <- scuttle::logNormCounts(sce)
sce

lun.sdrf <- file.path("~/Téléchargements/E-MTAB-5522(1)/E-MTAB-5522.sdrf.txt")
coldata <- read.delim(lun.sdrf, check.names=FALSE)

# Only keeping the cells involved in the count matrix in 'mat'.
coldata <- coldata[coldata[,"Derived Array Data File"]=="counts_Calero_20160113.tsv",]

# Only keeping interesting columns, and setting the library names as the row names.
coldata <- DataFrame(
  genotype=coldata[,"Characteristics[genotype]"],
  phenotype=coldata[,"Characteristics[phenotype]"],
  spike_in=coldata[,"Factor Value[spike-in addition]"],
  row.names=coldata[,"Source Name"]
)

coldata


sce <- SingleCellExperiment(assays = list(counts=mat), colData=coldata)
sce

colData(sce)

sce <- SingleCellExperiment(list(counts=mat))
colData(sce) <- coldata
sce

stopifnot(identical(rownames(coldata), colnames(mat)))







first.10 <- sce[,1:10]
ncol(counts(first.10)) # only 10 columns.
colData(first.10) # only 10 rows.
wt.only <- sce[, sce$phenotype == "wild type phenotype"]
ncol(counts(wt.only))
colData(wt.only)


sce <- scater::logNormCounts(sce)
sce <- scater::runPCA(sce)
dim(reducedDim(sce, "PCA"))
sce <- scater::runTSNE(sce, perplexity = 0.1)
head(reducedDim(sce, "TSNE"))

reducedDims(sce)
reducedDim(sce, "TSNE")


u <- uwot::umap(t(logcounts(sce)), n_neighbors = 2)
reducedDim(sce, "UMAP_uwot") <- u
reducedDims(sce) # Now stored in the object.

colLabels(sce) <- scran::clusterCells(sce, use.dimred="PCA")
table(colLabels(sce))


saveRDS(sce, file = "sce_test.rds")
sce <- readRDS("sce_test.rds")
reducedDim(sce)
reducedDim(sce, "UMAP_uwot")
emb <- as.data.frame(reducedDim(sce, "UMAP_uwot"))
colnames(emb)[1:2] <- c("Dim1", "Dim2")
exprs <- assay(sce, "logcounts")["ENSMUSG00000064341", , drop=FALSE]
transposed.exprs <- as.data.frame(t(exprs))
# Long format
df <- cbind(
  emb,
  t(exprs)
)
df_long <- tidyr::pivot_longer(
  df,
  cols = all_of("ENSMUSG00000064341"),
  names_to = "gene",
  values_to = "expression"
)

library(Seurat)
seurat_obj <- readRDS("~/Documents/Programmation/Heimdall/Single_cell/mon_objet_seurat_umap.rds")
seurat2sce <- as.SingleCellExperiment(seurat_obj)
saveRDS(seurat2sce, file = "seurat2sce.rds")
seurat2sce <- readRDS("~/Documents/Heimdall/seurat2sce.rds")










colnames_metadata_list <- colnames(colData(seurat2sce))
obj <- seurat2sce

nested_list_metadata <- setNames(
  lapply(colnames_metadata_list, 
         function(col) {  unique(obj[[col]])  }   ),
  colnames_metadata_list
  )

  
  
  
  
  lapply(colnames_metadata_list, function(col) {
  list(unique_values = unique(colnames_metadata_list[[col]]))
})
names(result) <- names(df)
