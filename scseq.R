install.packages('BiocManager')
BiocManager::install('multtest')
BiocManager::install("tximport")
BiocManager::install('biomaRt')
BiocManager::install('org.Hs.eg.db')
BiocManager::install("EnsDb.Hsapiens.v79")
BiocManager::install("AnnotationDbi")
install.packages('Seurat')
install.packages("RColorBrewer")
library(Seurat)
library(tximport)
library(fishpond)
library(AnnotationDbi)
library(EnsDb.Hsapiens.v79)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
#devtools::install_github(repo = "mikelove/tximport", dependencies = T)
devtools::install_github("hrbrmstr/waffle")
library(waffle)
library(magrittr)
library(hrbrthemes)
library(ggplot2)
library(dplyr)
library(waffle)
hrbrthemes::import_roboto_condensed()

alevin1 <- tximport(files = "/projectnb/bf528/users/group4/project4/code/salmon_04_out/alevin/quants_mat.gz", type = "alevin")
alevin2 <- tximport(files = "/projectnb/bf528/users/group4/project4/code/salmon_05_out/alevin/quants_mat.gz", type = "alevin")
alevin3 <- tximport(files = "/projectnb/bf528/users/group4/project4/code/salmon_06_out/alevin/quants_mat.gz", type = "alevin")

ids <- alevin3$counts@Dimnames[1]
ids <- unlist(ids, use.names=FALSE)
ids <- gsub("\\.[0-9]*$", "", ids)
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ids, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
id_list <- as.data.frame(ids)
mergeids <- merge(x = id_list, y = geneIDs, by.x = "ids", by.y = "GENEID", all.x = TRUE)
mergeids <-  mergeids[match(id_list$ids, mergeids$ids), ]
mergeids <- mergeids$SYMBOL
mergeids[is.na(mergeids)] <- "NA"
mergeids <- list(mergeids)
alevin3$counts@Dimnames[1] <- mergeids

my_data1 <- CreateSeuratObject(counts = alevin1$counts, project = "proj4", min.cells = 2, min.features = 200)
my_data2 <- CreateSeuratObject(counts = alevin2$counts, project = "proj4", min.cells = 2, min.features = 200)
my_data3 <- CreateSeuratObject(counts = alevin3$counts, project = "proj4", min.cells = 2, min.features = 200)

data_merge <- merge(my_data1, y = c(my_data2, my_data3), project = "proj4")


data_merge[["percent.mt"]] <- PercentageFeatureSet(data_merge, pattern = "^MT-")
VlnPlot(data_merge, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
merge1 <- subset(data_merge, subset = nFeature_RNA > 200 & nFeature_RNA < 3800 & percent.mt <20)
merge1 <- NormalizeData(merge1, normalization.method = "LogNormalize", scale.factor = 10000)
merge1 <- FindVariableFeatures(merge1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(merge1)
merge2 <- ScaleData(merge1, features = all.genes)
merge2 <- RunPCA(merge2, features = VariableFeatures(object = merge2))
ElbowPlot(merge2)
merge2 <- FindNeighbors(merge2, dims = 1:10)
merge2 <- FindClusters(merge2, resolution = 0.85)

saveRDS(merge2, file="/projectnb/bf528/users/group4/project4/code/seurat.rda")
load_test <- readRDS("/projectnb/bf528/project_4_scrnaseq/GSM2230760_seurat.rda")
load_test2 <- readRDS("/projectnb/bf528/users/group4/project4/code/seurat.rda")

clusts <- Idents(merge2)
clusts_df <- as.data.frame(table(clusts))
cluster <- c(1:15)
freq <-  clusts_df$Freq
clusts_df <- data.frame(cluster,freq)

clust_waffle <- c(
  `Cluster 1\n(441)` = 441,`Cluster 2\n(416)` = 416,`Cluster 3\n(395)` = 395,`Cluster 4\n(360)` = 360,
  `Cluster 5\n(327)` = 327,`Cluster 6\n(308)` = 308,`Cluster 7\n(270)` = 270,`Cluster 8\n(246)` = 246,
  `Cluster 9\n(240)` = 240,`Cluster 10\n(150)` = 150,`Cluster 11\n(136)` = 136,`Cluster 12\n(102)` = 102,
  `Cluster 13\n(101)` = 101,`Cluster 14\n(47)` = 47,`Cluster 15\n(44)` = 44
)

waffle_test <- waffle(clust_waffle/10, row=10, size = 0.5,  
               colors = c("#FFFF00", "#1879bf", "#00FF00","#FF00FF", "#FF0000", "#990099","#006600", "#FF8000", 
                          "#FFFF00", "#1879bf", "#00FF00","#FF00FF", "#FF0000", "#990099","#006600", "#C0C0C0"), 
                          legend_pos = "bottom", title = "Proportion of Cells by Cluster")
waffle_test
