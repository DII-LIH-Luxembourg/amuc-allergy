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

metadata_filename <- "LP_metadata.xlsx"
md <- read_excel(metadata_filename)
## Make sure condition variables are factors with the right levels
md$Diet <- factor(md$Diet)
md$Microbiota <- factor(md$Microbiota)

## Define colors for conditions
color_conditions <- c("#6A3D9A", "#FF7F00", "#1965B0", "#7BAFDE", "#008000", "#800000")
names(color_conditions) <- levels(md$Diet)

panel_filename <- "Panel_V2.csv"
panel <- read.csv(panel_filename)

setwd("PATH_TO_YOUR_DATA")
fcs_raw <- read.flowSet(md$`File name`, transformation = FALSE,
                        truncate_max_range = FALSE)
fcs_raw

# use both, type and state markers for large clustering
markers <- panel$Marker[panel$Type %in% c("type","state")]

# Replace problematic characters
# panel$antigen <- gsub("-", "_", panel$Antigen)

panel_fcs <- pData(parameters(fcs_raw[[1]]))
head(panel_fcs)
# Spot checks
all(markers %in% panel_fcs$desc)

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
sample_ids <- rep(md$`Mouse ID`, fsApply(fcs_raw, nrow))

ggdf <- data.frame(sample_id = sample_ids, expr, check.names = FALSE)
ggdf <- melt(ggdf, id.var = "sample_id",
             value.name = "expression", variable.name = "antigen")
mm <- match(ggdf$sample_id, md$`Mouse ID`)
ggdf$condition <- md$Diet[mm]
## ----------------------------------------------------------------------------
## FlowSOM
## ----------------------------------------------------------------------------
fsom <- ReadInput(fcs, transform = FALSE, scale = FALSE)
# set.seed(566961715)
set.seed(1234)
som <- BuildSOM(fsom, colsToUse = markers, xdim=20, ydim=20, rlen=50)
saveRDS(som, "../../results/data/som_all_400.rds")
## Get the cell clustering into 400 SOM codes
cell_clustering_som <- som$map$mapping[,1]

# Not using any meta-clustering
cell_clustering1 <- cell_clustering_som

## ----------------------------------------------------------------------------
## Heatmap of the median marker intensities of the lineage markers
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

write.csv(expr_median, "../../results/data/HM_data_400_median_expr.csv")


e01 = expr01[, markers]
expr01_median <- data.frame(e01, cell_clustering = cell_clustering1) %>%
  group_by(cell_clustering) %>% summarize_all(funs(median))

expr01_median$proportion <- clustering_prop

write.csv(expr01_median, "../../results/data/HM_data_400_0_to_1.csv")

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
mm <- match(dr_umap$sample_id, md$`Mouse ID`)
dr_umap$condition <- md$Diet[mm]
dr_umap$cell_clustering1 <- factor(cell_clustering1[tsne_inds])
dr_umap$batch <- factor(md$Batch[mm])
dr_umap$microbiota <- factor(md$Microbiota[mm])
head(dr_umap)

## Plot t-SNE colored by clusters
ggp_umap <- ggplot(dr_umap, aes(x = tSNE1, y = tSNE2, color = cell_clustering1)) +
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
# save.image(file = "<PATH_TO_YOUR_DATA>/data/data_clustered_400_and_UMAP.RData")

load(file = "<PATH_TO_YOUR_DATA>/data/data_clustered_400_and_UMAP.RData")

## make sure to use the 400 clusters
cell_clustering1 <- cell_clustering_som
## ----------------------------------------------------------------------------
## Main cluster merging:
## ----------------------------------------------------------------------------
cluster_merging1_filename <- "../../results/data/400_clusters_MB.xlsx"

cluster_merging1 <- read_excel(cluster_merging1_filename)
head(cluster_merging1)

cluster_merging1$population <- factor(cluster_merging1$population)
cluster_merging1$umap <- factor(cluster_merging1$umap)
cluster_merging1$color <- factor(cluster_merging1$color)

## New clustering1m
mm <- match(cell_clustering1, cluster_merging1$cluster)
cell_clustering1m <- cluster_merging1$population[mm]
cell_clustering2m <- cluster_merging1$umap[mm]
# mm <- match(code_clustering1, cluster_merging1$cell_clustering)
# code_clustering1m <- cluster_merging1$population[mm]

## ----------------------------------------------------------------------------
## Update the plot with the new annotated cell populations.
## ----------------------------------------------------------------------------
dr_umap$cell_clustering2m <- cell_clustering2m[tsne_inds]
ggp_merged <- ggplot(dr_umap, aes(x = tSNE1, y = tSNE2, color = cell_clustering2m)) +
  geom_point(size = 0.1) +
  theme_bw() +
  theme(legend.text=element_text(size=8)) +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggp_merged

## ----------------------------------------------------------------------------
## Annotate according to conditions
## ----------------------------------------------------------------------------
## subset to SPF FF vs FR
# dr_umap_sub <- dr_umap[dr_umap$microbiota == "SPF", ]
dr_umap_sub <- dr_umap[dr_umap$microbiota == "14SM", ]
# dr_umap_sub <- dr_umap[dr_umap$microbiota == "GF", ]
dr_umap_sub <- dr_umap_sub[dr_umap_sub$condition %in% c("FF", "FR1"), ]
head(dr_umap_sub)
tail(dr_umap_sub)
ggp_merged <- ggplot(dr_umap_sub, aes(x = tSNE1, y = tSNE2, color = cell_clustering2m)) +
  geom_point(size = 0.4) +
  theme_bw() +
  theme(legend.text=element_text(size=8)) +
  scale_color_manual(values = color_clusters) +
  scale_x_reverse() +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggp_merged

ggp_merged <- ggplot(dr_umap_sub, aes(x = tSNE2, y = tSNE1, color = cell_clustering2m)) +
  geom_point(size = 0.4) +
  theme_bw() +
  theme(legend.text=element_text(size=8)) +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggp_merged

## Between conditions
ggp_merged + facet_wrap(~ condition)

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

ggplot(dr_umap, aes(x = tSNE1, y = tSNE2, color = `CD19`)) +
    geom_point(size = 0.2) +
    theme_bw() +
    scale_color_gradientn("CD19",
                          colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))

ggplot(dr_umap, aes(x = tSNE1, y = tSNE2, color = `CD45R`)) +
  geom_point(size = 0.2) +
  theme_bw() +
  scale_color_gradientn("CD45R",
                        colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))

## ----------------------------------------------------------------------------
mm <- cell_clustering1m != "Others"
# mm <- cell_clustering1m[cell_clustering1m == "Others"]
# mm <- cell_clustering1m[cell_clustering1m != "Others"]
sub_ids <- sample_ids[mm]
## remove other cells
expr_sub <- expr[mm, markers]
expr01_sub <- expr01[mm, markers]
cell_clustering1m_sub <- cell_clustering1m[mm]
dim(expr_sub)
dim(expr01_sub)

## ----------------------------------------------------------------------------
## Save and/or load work-space after initial clustering
## ----------------------------------------------------------------------------
save.image(file = "../../results/data/data_annotated.RData")

load(file = "../../results/data/data_annotated.RData")
## ----------------------------------------------------------------------------
## Update heatmap
dim(expr)
dim(expr_sub)
dim(expr01_sub)
length(cell_clustering1m_sub)
## ----------------------------------------------------------------------------
plot_clustering_heatmap_wrapper(expr = expr_sub,
                                expr01 = expr01_sub,
                                cell_clustering = cell_clustering1m_sub,
                                color_clusters = color_clusters)
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
colnames(expr) <- unname(colnames(expr))
e = expr[, markers]
head(e)
colnames(e)
colnames(expr)
expr_median <- data.frame(e, cell_clustering = cell_clustering1m) %>%
  group_by(cell_clustering) %>% summarize_all(funs(median))

clustering_table <- as.numeric(table(cell_clustering1m))
clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)

expr_median$proportion <- clustering_prop

# write.csv(expr_median, "./data/HM_data_annotated_median_expr.csv")

e01 = expr01[, markers]
expr01_median <- data.frame(e01, cell_clustering = cell_clustering1m) %>%
  group_by(cell_clustering) %>% summarize_all(funs(median))

expr01_median$proportion <- clustering_prop

# write.csv(expr01_median, "./data/HM_data_annotated_0_to_1.csv")

## ----------------------------------------------------------------------------
## cell_clustering2m
## for custom selection of phenotypes for umap
## ----------------------------------------------------------------------------
cell_clustering2m <- cluster_merging1$cat[mm]
head(cell_clustering2m)
tail(cell_clustering2m)
## ----------------------------------------------------------------------------
mm <- match(code_clustering1, cluster_merging1$cell_clustering)
code_clustering1m <- cluster_merging1$population[mm]

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
dr_umap$cell_clustering1m <- cell_clustering1m[tsne_inds]
## ----------------------------------------------------------------------------
## Update the plot with the new annotated cell populations.
## ----------------------------------------------------------------------------
dr_umap$cell_clustering2m <- cell_clustering2m[tsne_inds]
## subset to SPF FF vs FR
dr_umap_sub <- dr_umap[dr_umap$microbiota == "SPF", ]
dr_umap_sub <- dr_umap_sub[dr_umap_sub$condition %in% c("FF", "FR1"), ]
head(dr_umap_sub)
tail(dr_umap_sub)
ggp_merged <- ggplot(dr_umap_sub, aes(x = tSNE1, y = tSNE2, color = cell_clustering1m)) +
  geom_point(size = 0.4) +
  theme_bw() +
  theme(legend.text=element_text(size=8)) +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggp_merged

## Between conditions
ggp_merged + facet_wrap(~ condition)

## ----------------------------------------------------------------------------
## only umap SPF, FF vs FR
## ----------------------------------------------------------------------------
head(expr)
mb_ids <- md$`Mouse ID`[md$Microbiota == "SPF" & md$Diet %in% c("FF","FR1")]

expr_sub <- expr[sample_ids %in% mb_ids, ]
sample_ids_sub <- sample_ids[sample_ids %in% mb_ids]
cc_sub <- cell_clustering1m[sample_ids %in% mb_ids]
dim(expr_sub)
length(sample_ids_sub)
length(cc_sub)

# head(cell_clustering1m
sub_pop <- cc_sub %in% c("NK cells", "Mast cells", "ILC2", "Lti-like ILC3", "Th2 (CD44+ CD69+)", "Th1 cells (CD62L+CD44+)", "Th17 (CD44+)")
expr_sub <- expr_sub[sub_pop, ]
sample_ids_sub <- sample_ids_sub[sub_pop]
cc_sub <- cc_sub[sub_pop]

# rng <- colQuantiles(sub_expr, probs = c(0.01, 0.99))
# sub_expr01 <- t((t(sub_expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
# sub_expr01[sub_expr01 < 0] <- 0
# sub_expr01[sub_expr01 > 1] <- 1

## Find and skip duplicates
dups <- which(!duplicated(expr_sub[, markers]))
## Data subsampling: create indices by sample
inds <- split(1:length(sample_ids_sub), sample_ids_sub)
## How many cells to downsample per-sample
tsne_ncells <- pmin(table(sample_ids_sub), 1000)
## Get subsampled indices
set.seed(1234)
tsne_inds <- lapply(names(inds), function(i){
  s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)
  intersect(s, dups)
})
tsne_inds <- unlist(tsne_inds)
tsne_expr <- expr_sub[tsne_inds, markers]

X<-umap(tsne_expr)
# X<-umap(expr_sub)
head(X$layout)
head(X$data)

dr_umap <- data.frame(tSNE1 = X$layout[, 1], tSNE2 = X$layout[, 2],
                      X$data,check.names = FALSE)

dr_umap$sample_id <- sample_ids_sub[tsne_inds]
# dr_umap$sample_id <- sample_ids_sub
mm <- match(dr_umap$sample_id, md$`Mouse ID`)
dr_umap$condition <- md$Diet[mm]
dr_umap$cell_clustering1 <- factor(cc_sub[tsne_inds])
dr_umap$batch <- factor(md$Batch[mm])
dr_umap$microbiota <- factor(md$Microbiota[mm])
head(dr_umap)
tail(dr_umap)
dim(dr_umap)
## Plot t-SNE colored by clusters
ggp_umap <- ggplot(dr_umap, aes(x = tSNE1, y = tSNE2, color = cell_clustering1)) +
  geom_point(size = 0.1) +
  theme_bw() +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
ggp_umap

ggp_umap  + facet_wrap(~ condition)

###################################
mycolors <- read_excel("../../sub_clustering_V2/400_clusters/400_clusters_colors.xlsx")

# mm <- match(cell_clustering1m, mycolors$population)
#
# custom_colors <- mycolors$color[mm]

ll <- levels(cc_sub)
sort(ll)
# mm <- match(levels(cell_clustering1m), mycolors$population)
mm <- match(ll, mycolors$population)
custom_colors <- mycolors$color[mm]

mycol <- mycolors$color[mycolors$population]

head(dr_umap)
head(cell_clustering1m)
ggp_merged <- ggplot(dr_umap, aes(x = tSNE1, y = tSNE2, color = cell_clustering1)) +
  geom_point(size = 1.0) +
  theme_bw() +
  theme(legend.text=element_text(size=6)) +
  scale_color_manual(values = custom_colors) +
  guides(color = guide_legend(override.aes = list(size = 3)))


ggp_merged
ggp_merged + facet_wrap(~ condition)




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
expr_median <- data.frame(e, cell_clustering = cell_clustering1) %>%
  group_by(cell_clustering) %>% summarize_all(funs(median))

clustering_table <- as.numeric(table(cell_clustering1))
clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)

expr_median$proportion <- clustering_prop

write.csv(expr_median, "../../sub_clustering_V2/400_clusters/HM_data_annotated_median_expr.csv")

e01 = expr01[, markers]
expr01_median <- data.frame(e01, cell_clustering = cell_clustering1m) %>%
  group_by(cell_clustering) %>% summarize_all(funs(median))

expr01_median$proportion <- clustering_prop

write.csv(expr01_median, "../../sub_clustering_V2/400_clusters/HM_data_annotated_0_to_1.csv")


## ----------------------------------------------------------------------------
## DE marker expression
## ----------------------------------------------------------------------------
counts_table <- table(cell_clustering1m, sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

write.csv(props, "../../results/data/sample_proportions.csv")


