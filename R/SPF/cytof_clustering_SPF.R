# R Script for Paper Publication
# Author: Oliver Hunewald
# Date: 24.Aug.2023

# Acknowledgment: This workflow is inspired by the work of Nowicka et al.
# The original workflow can be found at https://f1000research.com/articles/6-748/v2

library(readxl)
library(tidyverse)
library(flowCore)
library(matrixStats)
library(ggplot2)
library(reshape2)
library(dplyr)
library(limma)
library(ggrepel)
library(FlowSOM)
library(umap)
library(cowplot)
library(ConsensusClusterPlus)

source("../helper_functions.R")
# --------------------------------------------------------------

metadata_filename <- "metadata.xlsx"
md <- read_excel(metadata_filename)

head(md)
## Make sure condition variables are factors with the right levels
md$diet <- factor(md$diet)
md$treatment <- factor(md$treatment)
head(data.frame(md))

## Define colors for conditions
color_conditions <- c("#6A3D9A", "#FF7F00", "#1965B0", "#7BAFDE", "#008000", "#800000")
names(color_conditions) <- levels(md$diet)

panel_filename <- "panel.csv"
panel <- read.csv(panel_filename)

setwd("PATH_TO_YOUR_DATA")
fcs_raw <- read.flowSet(md$file_name, transformation = FALSE,
                        truncate_max_range = FALSE)
fcs_raw


head(data.frame(panel))

markers <- panel$Marker[panel$Type %in% c("type","state")]

# Replace problematic characters
# panel$antigen <- gsub("-", "_", panel$Antigen)

panel_fcs <- pData(parameters(fcs_raw[[1]]))
head(panel_fcs)

# Spot checks
all(markers %in% unname(panel_fcs$desc))

## transform all markers
fcs <- fsApply(fcs_raw, function(x, cofactor=5){

  colnames(x) <- panel_fcs$desc
  expr <- exprs(x)
  # use unique function because markers can be lineage AND functional
  expr <- asinh(expr / cofactor)
  exprs(x) <- expr
  x
})
fcs

## Extract expression
expr <- fsApply(fcs, exprs)
dim(expr)

rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1

## Generate sample IDs corresponding to each cell in the 'expr' matrix
sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))

ggdf <- data.frame(sample_id = sample_ids, expr, check.names = FALSE)
ggdf <- melt(ggdf, id.var = "sample_id",
             value.name = "expression", variable.name = "antigen")
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- md$diet[mm]
## ----------------------------------------------------------------------------
## FlowSOM
## ----------------------------------------------------------------------------
fsom <- ReadInput(fcs, transform = FALSE, scale = FALSE)

set.seed(1234)
som <- BuildSOM(fsom, colsToUse = markers, xdim=20, ydim=20, rlen=50)
saveRDS(som, "<SET_YOUR_PATH>/som_all_400.rds")
## Get the cell clustering into 100 SOM codes
cell_clustering_som <- som$map$mapping[,1]
codes <- som$map$codes

cell_clustering1 <- code_clustering1[cell_clustering_som]

## NOT USING ANY META-CLUSTERING
cell_clustering1 <- cell_clustering_som

## ----------------------------------------------------------------------------
## Heatmap of the median marker intensities
## ----------------------------------------------------------------------------
plot_clustering_heatmap_wrapper(expr = expr[, markers],
                                expr01 = expr01[, markers],
                                cell_clustering = cell_clustering1,
                                color_clusters = color_clusters)

## ----------------------------------------------------------------------------
## Save the heatmap data as table
## ----------------------------------------------------------------------------
# Calculate the median expression
e = expr[, markers]
expr_median <- data.frame(e, cell_clustering = cell_clustering1) %>%
  group_by(cell_clustering) %>% summarize_all(funs(median))

clustering_table <- as.numeric(table(cell_clustering1))
clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)

expr_median$proportion <- clustering_prop

# write.csv(expr_median, "../HM_data_400_median_expr.csv")


e01 = expr01[, markers]
expr01_median <- data.frame(e01, cell_clustering = cell_clustering1) %>%
  group_by(cell_clustering) %>% summarize_all(funs(median))

expr01_median$proportion <- clustering_prop

# write.csv(expr01_median, "../HM_data_400_0_to_1.csv")

## ----------------------------------------------------------------------------
## UMAP
## ----------------------------------------------------------------------------
## Find and skip duplicates
dups <- which(!duplicated(expr[, markers]))
## Data subsampling: create indices by sample
inds <- split(1:length(sample_ids), sample_ids)
## How many cells to downsample per-sample
tsne_ncells <- pmin(table(sample_ids), 4000)
## Get subsampled indices
set.seed(1234)
tsne_inds <- lapply(names(inds), function(i){
  s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)
  intersect(s, dups)
})
tsne_inds <- unlist(tsne_inds)
tsne_expr <- expr[tsne_inds, markers]

X<-umap(tsne_expr)
head(X$layout)
head(X$data)

dr_umap <- data.frame(tSNE1 = X$layout[, 1], tSNE2 = X$layout[, 2],
                      X$data,check.names = FALSE)

dr_umap$sample_id <- sample_ids[tsne_inds]
mm <- match(dr_umap$sample_id, md$sample_id)
dr_umap$condition <- md$diet[mm]
dr_umap$cell_clustering1 <- factor(cell_clustering1[tsne_inds])
dr_umap$batch <- factor(md$Batch[mm])
# dr_umap$microbiota <- factor(md$Microbiota[mm])
dr_umap$treatment <- factor(md$treatment[mm])
head(dr_umap)

dr_umap$condition
## Plot t-SNE colored by clusters
ggp_umap <- ggplot(dr_umap, aes(x = tSNE1, y = tSNE2, color = sample_id)) +
  geom_point(size = 0, alpha=0.6) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
ggp_umap

ggp_umap  + facet_wrap(~ condition)

## check if there is any batch effect
ggp_umap <- ggplot(dr_umap, aes(x = tSNE1, y = tSNE2, color = batch)) +
  geom_point(size = 0, alpha=0.6) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
ggp_umap
ggp_umap  + facet_wrap(~ batch)

## ----------------------------------------------------------------------------
## Save and/or load work-space after initial clustering
## ----------------------------------------------------------------------------
# save.image(file = "<SET_YOUR_PATH>/data_clustered_400_and_UMAP.RData")

# load(file = "<SET_YOUR_PATH>/data_clustered_400_and_UMAP.RData")

## make sure the 400 clusters are loaded
# cell_clustering1 <- cell_clustering_som
## ----------------------------------------------------------------------------
## Main cluster merging:
## ----------------------------------------------------------------------------
cluster_merging1_filename <- "../400_clusters_MB.xlsx"

cluster_merging1 <- read_excel(cluster_merging1_filename)
head(cluster_merging1)

cluster_merging1$population <- factor(cluster_merging1$population)

## New clustering1m
mm <- match(cell_clustering1, cluster_merging1$cell_clustering)
cell_clustering1m <- cluster_merging1$population[mm]

## ----------------------------------------------------------------------------
## 2nd merging for umap visualization
## ----------------------------------------------------------------------------
mm <- match(cell_clustering1, cluster_merging1$cell_clustering)
cell_clustering2m <- cluster_merging1$`umap cluster`[mm]

## ----------------------------------------------------------------------------
## Update the t-SNE plot with the new annotated cell populations.
## ----------------------------------------------------------------------------
dr_umap$cell_clustering1m <- cell_clustering1m[tsne_inds]
ggp_merged <- ggplot(dr_umap, aes(x = tSNE1, y = tSNE2, color = cell_clustering1m)) +
  geom_point(size = 0.1) +
  theme_bw() +
  theme(legend.text=element_text(size=8)) +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggp_merged

## Between conditions
ggp_merged + facet_wrap(~ condition)

## ----------------------------------------------------------------------------
## Update the umap with reduced populations
## ----------------------------------------------------------------------------
dr_umap$cell_clustering2m <- cell_clustering2m[tsne_inds]
ggp_merged <- ggplot(dr_umap, aes(x = tSNE1, y = tSNE2, color = cell_clustering2m)) +
  geom_point(size = 0.2) +
  theme_bw() +
  theme(legend.text=element_text(size=8)) +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggp_merged

## Between conditions
ggp_merged + facet_wrap(~ condition)

head(dr_umap)
head(md)
dr_umap$treatment[1]
## ----------------------------------------------------------------------------
# dr_umap_sub <- dr_umap[dr_umap$treatment == "OVA CTX", ]
# dr_umap_sub <- dr_umap[dr_umap$treatment == "PBS", ]
dr_umap_sub <- dr_umap[dr_umap$treatment == "CTX", ]
ggp_merged <- ggplot(dr_umap_sub, aes(x = tSNE1, y = tSNE2, color = cell_clustering2m)) +
  geom_point(size = 0.4) +
  theme_bw() +
  theme(legend.text=element_text(size=8)) +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggp_merged

## Between conditions
ggp_merged + facet_wrap(~ condition)

ggplot(dr_umap, aes(x = tSNE1, y = tSNE2, color = `Tbet`)) +
  geom_point(size = 0.2) +
  theme_bw() +
  scale_color_gradientn("Tbet",
                        colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))

## ----------------------------------------------------------------------------
## Update heatmap
## ----------------------------------------------------------------------------
plot_clustering_heatmap_wrapper(expr = expr[, markers],
                                expr01 = expr01[, markers],
                                cell_clustering = cell_clustering1m,
                                color_clusters = color_clusters)
## ----------------------------------------------------------------------------
## Export the heatmap data
## ----------------------------------------------------------------------------
e = expr[, markers]
expr_median <- data.frame(e, cell_clustering = cell_clustering1m) %>%
  group_by(cell_clustering) %>% summarize_all(funs(median))

clustering_table <- as.numeric(table(cell_clustering1m))
clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)

expr_median$proportion <- clustering_prop

write.csv(expr_median, "../figures/HM_data_annotated_median_expr.csv")

e01 = expr01[, markers]
expr01_median <- data.frame(e01, cell_clustering = cell_clustering1m) %>%
  group_by(cell_clustering) %>% summarize_all(funs(median))

expr01_median$proportion <- clustering_prop

write.csv(expr01_median, "../figures/HM_data_annotated_0_to_1.csv")


## ----------------------------------------------------------------------------
## DE marker expression
## ----------------------------------------------------------------------------
counts_table <- table(cell_clustering1m, sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

write.csv(props, "../sample_proportions.csv")












