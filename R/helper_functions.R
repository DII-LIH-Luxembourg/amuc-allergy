library(RColorBrewer)
library(pheatmap)
# Define cluster colors (here there are 50 colors)
color_clusters <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00",
                    "#000000", "#0000ff", "#800080", "#ffb6c1", "#003366",
                    "#00ff00", "#666666", "#b0e0e6", "#c39797", "#66cdaa",
                    "#ff6666", "#ffc3a0", "#ff00ff", "#333333", "#cccccc",
                    "#088da5", "#c0d6e4", "#8b0000", "#660066", "#ff7f50",
                    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00",
                    "#000000", "#0000ff", "#800080", "#ffb6c1", "#003366",
                    "#00ff00", "#666666", "#b0e0e6", "#c39797", "#66cdaa",
                    "#ff6666", "#ffc3a0", "#ff00ff", "#333333", "#cccccc",
                    "#088da5", "#c0d6e4", "#8b0000", "#660066", "#ff7f50")


plot_clustering_heatmap_wrapper <- function(expr, expr01,
                                            cell_clustering, color_clusters, cluster_merging = NULL){
  # Calculate the median expression
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  ## ****************************************************************************
  ## "X" in front of each antigen
  ## TODO: locate where X comes from, or workaround by removing X
  ## ****************************************************************************
  colnames(expr_median) <- gsub('^X', '', colnames(expr_median))
  colnames(expr01_median) <- gsub('^X', '', colnames(expr01_median))
  
  
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  # Colors for the heatmap
  color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  labels_row <- paste0(expr01_median$cell_clustering, " (", clustering_prop ,"%)")
  # Annotation for the original clusters
  annotation_row <- data.frame(Cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  # Annotation for the merged clusters
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Cluster_merging <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Cluster_merging <- color_clusters2
  }
  
  pheatmap(expr_heat, color = color_heat, cluster_cols = FALSE,
           cluster_rows = cluster_rows, labels_row = labels_row,
           display_numbers = TRUE, number_color = "black",
           fontsize = 8, fontsize_number = 6, legend_breaks = legend_breaks,
           annotation_row = annotation_row, annotation_colors = annotation_colors)
}


plot_clustering_SOM_heatmap_wrapper <- function(expr, expr01,som,
                                            cell_clustering, color_clusters, cluster_merging = NULL){
  # Calculate the median expression
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  ## ****************************************************************************
  ## "X" in front of each antigen
  ## TODO: locate where X comes from, or workaround by removing X
  ## ****************************************************************************
  colnames(expr_median) <- gsub('^X', '', colnames(expr_median))
  colnames(expr01_median) <- gsub('^X', '', colnames(expr01_median))
  
  
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  # Colors for the heatmap
  color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  labels_row <- paste0(expr01_median$cell_clustering, " (", clustering_prop ,"%)")
  # Annotation for the original clusters
  annotation_row <- data.frame(Cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  color_clusters1 <- color_clusters[1:nlevels(annotation_row$Cluster)]
  names(color_clusters1) <- levels(annotation_row$Cluster)
  annotation_colors <- list(Cluster = color_clusters1)
  # Annotation for the merged clusters
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$Cluster_merging <- cluster_merging$new_cluster
    color_clusters2 <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters2) <- levels(cluster_merging$new_cluster)
    annotation_colors$Cluster_merging <- color_clusters2
  }
  
  pheatmap(expr_heat, color = color_heat, cluster_cols = FALSE,
           cluster_rows = cluster_rows, labels_row = labels_row,
           display_numbers = TRUE, number_color = "black",
           fontsize = 8, fontsize_number = 6, legend_breaks = legend_breaks,
           annotation_row = annotation_row, annotation_colors = annotation_colors)
}

# install.packages("ggridges")
library(ggridges)
plot_clustering_distr_wrapper <- function(expr, cell_clustering){
  # Calculate the median expression
  cell_clustering <- factor(cell_clustering)
  expr_median <- data.frame(expr, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  ## ****************************************************************************
  ## "X" in front of each antigen
  ## TODO: locate where X comes from, or workaround by removing X
  ## ****************************************************************************
  colnames(expr_median) <- gsub('^X', '', colnames(expr_median))
  
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  # Calculate cluster frequencies
  freq_clust <- table(cell_clustering)
  freq_clust <- round(as.numeric(freq_clust)/sum(freq_clust)*100, 2)
  cell_clustering <- factor(cell_clustering,
                            labels = paste0(levels(cell_clustering), " (", freq_clust, "%)"))
  ### Data organized per cluster
  ggd <- melt(data.frame(cluster = cell_clustering, expr,check.names = FALSE),
              id.vars = "cluster", value.name = "expression",
              variable.name = "antigen")
  ## ****************************************************************************
  ggd$antigen <- gsub('^X', '', ggd$antigen)
  
  ggd$antigen <- factor(ggd$antigen, levels = colnames(expr))
  ggd$reference <- "no"
  ### The reference data
  ggd_bg <- ggd
  ggd_bg$cluster <- "reference"
  ggd_bg$reference <- "yes"
  ggd_plot <- rbind(ggd, ggd_bg)
  ggd_plot$cluster <- factor(ggd_plot$cluster,
                             levels = c(levels(cell_clustering)[rev(cluster_rows$order)], "reference"))
  ggplot() +
    geom_density_ridges(data = ggd_plot, aes(x = expression, y = cluster,
                                             color = reference, fill = reference), alpha = 0.3) +
    facet_wrap( ~ antigen, scales = "free_x", nrow = 2) +
    theme_ridges() +
    theme(axis.text = element_text(size = 7),
          strip.text = element_text(size = 7), legend.position = "none")
}


# source("https://bioconductor.org/biocLite.R")
# biocLite("ComplexHeatmap")
library(ComplexHeatmap)
plot_clustering_heatmap_wrapper2 <- function(expr, expr01,
                                             lineage_markers, functional_markers = NULL, sample_ids = NULL,
                                             cell_clustering, color_clusters, cluster_merging = NULL,
                                             plot_cluster_annotation = TRUE){
  # Calculate the median expression of lineage markers
  expr_median <- data.frame(expr[, lineage_markers],
                            cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  
  expr01_median <- data.frame(expr01[, lineage_markers],
                              cell_clustering = cell_clustering, check.names = FALSE) %>%
    group_by(cell_clustering) %>% summarize_all(funs(median))
  
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  clustering_prop <- round(clustering_table / sum(clustering_table) * 100, 2)
  # Sort the cell clusters with hierarchical clustering
  d <- dist(expr_median[, lineage_markers], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  expr_heat <- as.matrix(expr01_median[, lineage_markers])
  # Median expression of functional markers in each sample per cluster
  expr_median_sample_cluster_tbl <- data.frame(expr01[, functional_markers,
                                                      drop = FALSE], sample_id = sample_ids, cluster = cell_clustering, check.names = FALSE) %>%
    group_by(sample_id, cluster) %>% summarize_all(funs(median))
  # Colors for the heatmap
  color_heat <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  legend_breaks = seq(from = 0, to = 1, by = 0.2)
  labels_row <- paste0(expr01_median$cell_clustering, " (", clustering_prop , "%)")
  ### Annotation for the original clusters
  annotation_row1 <- data.frame(Cluster = factor(expr01_median$cell_clustering))
  color_clusters1 <- color_clusters[1:nlevels(annotation_row1$Cluster)]
  names(color_clusters1) <- levels(annotation_row1$Cluster)
  ### Annotation for the merged clusters
  if(!is.null(cluster_merging)){
    mm <- match(annotation_row1$Cluster, cluster_merging$original_cluster)
    annotation_row2 <- data.frame(Cluster_merging =
                                    factor(cluster_merging$new_cluster[mm]))
    color_clusters2 <- color_clusters[1:nlevels(annotation_row2$Cluster_merging)]
    names(color_clusters2) <- levels(annotation_row2$Cluster_merging)
  }
  ### Heatmap annotation for the original clusters
  ha1 <- Heatmap(annotation_row1, name = "Cluster",
                 col = color_clusters1, cluster_columns = FALSE,
                 cluster_rows = cluster_rows, row_dend_reorder = FALSE,
                 show_row_names = FALSE, width = unit(0.5, "cm"),
                 rect_gp = gpar(col = "grey"))
  
  ### Heatmap annotation for the merged clusters
  if(!is.null(cluster_merging)){
    ha2 <- Heatmap(annotation_row2, name = "Cluster \nmerging",
                   col = color_clusters2, cluster_columns = FALSE,
                   cluster_rows = cluster_rows, row_dend_reorder = FALSE,
                   show_row_names = FALSE, width = unit(0.5, "cm"),
                   rect_gp = gpar(col = "grey"))
  }
  ### Cluster names and sizes - text
  ha_text <- rowAnnotation(text = row_anno_text(labels_row,
                                                gp = gpar(fontsize = 6)), width = max_text_width(labels_row))
  ### Cluster sizes - barplot
  ha_bar <- rowAnnotation("Frequency (%)" = row_anno_barplot(
    x = clustering_prop, border = FALSE, axis = TRUE,
    axis_gp = gpar(fontsize = 5), gp = gpar(fill = "#696969", col = "#696969"),
    bar_width = 0.9), width = unit(0.7, "cm"), show_annotation_name = TRUE,
    annotation_name_rot = 0, annotation_name_offset = unit(5, "mm"),
    annotation_name_gp = gpar(fontsize = 7))
  ### Heatmap for the lineage markers
  ht1 <- Heatmap(expr_heat, name = "Expr", column_title = "Lineage markers",
                 col = color_heat, cluster_columns = FALSE, cluster_rows = cluster_rows,
                 row_dend_reorder = FALSE, heatmap_legend_param = list(at = legend_breaks,
                                                                       labels = legend_breaks, color_bar = "continuous"),
                 show_row_names = FALSE, row_dend_width = unit(2, "cm"),
                 rect_gp = gpar(col = "grey"), column_names_gp = gpar(fontsize = 8))
  if(plot_cluster_annotation){
    draw_out <- ha1
  }else{
    draw_out <- NULL
  }
  if(!is.null(cluster_merging)){
    draw_out <- draw_out + ha2 + ht1 + ha_bar + ha_text
  }else{
    draw_out <- draw_out + ht1 + ha_bar + ha_text
  }
  ### Heatmaps for the signaling markers
  if(!is.null(functional_markers)){
    for(i in 1:length(functional_markers)){
      ## Rearange so the rows represent clusters
      expr_heat_fun <- as.matrix(dcast(expr_median_sample_cluster_tbl[,
                                                                      c("sample_id", "cluster", functional_markers[i])],
                                       cluster ~ sample_id, value.var = functional_markers[i])[, -1])
      draw_out <- draw_out + Heatmap(expr_heat_fun,
                                     column_title = functional_markers[i], col = color_heat,
                                     cluster_columns = FALSE, cluster_rows = cluster_rows,
                                     row_dend_reorder = FALSE, show_heatmap_legend = FALSE,
                                     show_row_names = FALSE, rect_gp = gpar(col = "grey"),
                                     column_names_gp = gpar(fontsize = 8))
    }
  }
  draw(draw_out, row_dend_side = "left")
}


## ----------------------------------------------------------------------------
differential_expression_wrapper <- function(expr_median, md, model = "lmer", formula, K){
  
  ## Fit LMM or LM for each marker separately
  fit_gaussian <- lapply(1:nrow(expr_median), function(i){
    
    data_tmp <- data.frame(y = as.numeric(expr_median[i, md$sample_id]), md)
    
    switch(model, 
           lmer = {
             fit_tmp <- lmer(formula, data = data_tmp)
           },
           lm = {
             fit_tmp <- lm(formula, data = data_tmp)
           })
    
    ## Fit contrasts one by one
    out <- apply(K, 1, function(k){
      contr_tmp <- glht(fit_tmp, linfct = matrix(k, 1))
      summ_tmp <- summary(contr_tmp)
      out <- c(summ_tmp$test$coefficients, summ_tmp$test$pvalues)
      names(out) <- c("coeff", "pval")
      return(out)
    })
    
    return(t(out))
  })
  
  ### Extract fitted contrast coefficients
  coeffs <- lapply(fit_gaussian, function(x){
    x[, "coeff"]
  })
  coeffs <- do.call(rbind, coeffs)
  colnames(coeffs) <- paste0("coeff_", contrast_names)
  rownames(coeffs) <- rownames(expr_median)
  
  ### Extract p-values
  pvals <- lapply(fit_gaussian, function(x){
    x[, "pval"]
  })
  pvals <- do.call(rbind, pvals)
  colnames(pvals) <- paste0("pval_", contrast_names)
  rownames(pvals) <- rownames(expr_median)
  
  ### Adjust the p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", contrast_names)
  
  return(list(coeffs = coeffs, pvals = pvals, adjp = adjp))
}

plot_differential_heatmap_wrapper <- function(expr_norm, sign_adjp, 
                                              condition, color_conditions, th = 2.5){
  ## Order samples by condition
  oo <- order(condition)
  condition <- condition[oo]
  expr_norm <- expr_norm[, oo, drop = FALSE]
  
  ## Create the row labels with adj p-values and other objects for pheatmap
  labels_row <- paste0(rownames(expr_norm), " (", sprintf( "%.02e", sign_adjp), ")")
  labels_col <- colnames(expr_norm)
  annotation_col <- data.frame(condition = factor(condition))
  rownames(annotation_col) <- colnames(expr_norm)
  annotation_colors <- list(condition = color_conditions)
  color <- colorRampPalette(c("#87CEFA", "#56B4E9", "#56B4E9", "#0072B2", 
                              "#000000", "#D55E00", "#E69F00", "#E69F00", "#FFD700"))(100)
  breaks = seq(from = -th, to = th, length.out = 101)
  legend_breaks = seq(from = -round(th), to = round(th), by = 1)
  gaps_col <- as.numeric(table(annotation_col$condition))
  
  pheatmap(expr_norm, color = color, breaks = breaks, 
           legend_breaks = legend_breaks, cluster_cols = FALSE, cluster_rows = FALSE, 
           labels_col = labels_col, labels_row = labels_row, gaps_col = gaps_col, 
           annotation_col = annotation_col, annotation_colors = annotation_colors, 
           fontsize = 8)
}


##
normalization_wrapper <- function(expr, th = 2.5){
  expr_norm <- apply(expr, 1, function(x){ 
    sdx <- sd(x, na.rm = TRUE)
    if(sdx == 0){
      x <- (x - mean(x, na.rm = TRUE))
    }else{ 
      x <- (x - mean(x, na.rm = TRUE)) / sdx
    }
    x[x > th] <- th
    x[x < -th] <- -th
    return(x)
  })
  expr_norm <- t(expr_norm)
}


differential_abundance_wrapper <- function(counts, md, formula, K){
  ## Fit the GLMM for each cluster separately
  ntot <- colSums(counts)
  fit_binomial <- lapply(1:nrow(counts), function(i){
    
    data_tmp <- data.frame(y = as.numeric(counts[i, md$sample_id]), 
                           total = ntot[md$sample_id], md)
    fit_tmp <- glmer(formula, weights = total, family = binomial, 
                     data = data_tmp)
    ## Fit contrasts one by one
    out <- apply(K, 1, function(k){
      contr_tmp <- glht(fit_tmp, linfct = matrix(k, 1))
      summ_tmp <- summary(contr_tmp)
      out <- c(summ_tmp$test$coefficients, summ_tmp$test$pvalues)
      names(out) <- c("coeff", "pval")
      return(out)
    })
    return(t(out))
  })
  ### Extract fitted contrast coefficients
  coeffs <- lapply(fit_binomial, function(x){
    x[, "coeff"]
  })
  coeffs <- do.call(rbind, coeffs)
  colnames(coeffs) <- paste0("coeff_", contrast_names)
  rownames(coeffs) <- rownames(counts)
  
  ### Extract p-values
  pvals <- lapply(fit_binomial, function(x){
    x[, "pval"]
  })
  pvals <- do.call(rbind, pvals)
  colnames(pvals) <- paste0("pval_", contrast_names)
  rownames(pvals) <- rownames(counts)
  
  ### Adjust the p-values
  adjp <- apply(pvals, 2, p.adjust, method = "BH")
  colnames(adjp) <- paste0("adjp_", contrast_names)
  
  return(list(coeffs = coeffs, pvals = pvals, adjp = adjp))
}
