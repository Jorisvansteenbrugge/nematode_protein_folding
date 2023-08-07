pacman::p_load(tidyverse,ggplot2, ComplexHeatmap, grattantheme)

#load("processed_data//2023-08-01_folding_fatcat_output_all_pvalueadjusted_inc_clusters.RData")





# Function to calculate the overlap between two vectors
calculate_overlap <- function(vec1, vec2) {
  intersect_length <- length(intersect(vec1, vec2))
  return(intersect_length)
}

calculate_overlap_percentage <- function(vec1, vec2) {
  intersect_length <- length(intersect(vec1, vec2))
  union_length <- length(union(vec1, vec2))
  overlap_percentage <- (intersect_length / union_length) * 100
  return(overlap_percentage)
}



Cluster_overlap_heatmap <- function(comparisons.filtered, k_split = 5, seed = 10) {

  u_clusters <- unique(comparisons.filtered$Cluster)

  cluster_GO_terms <- lapply(u_clusters, function(cluster){
    cluster_pool <- comparisons.filtered |>
      filter(Cluster == cluster)|>
      select(Query_GO_ID, Subject_GO_ID) |>
      unlist() |>
      stringr::str_split(pattern=";") |>
      unlist() |>
      gsub(pattern = "\\s+", replacement = "") |>
      unique()
    return(cluster_pool)
  }
  )

  names(cluster_GO_terms) <- u_clusters


                                        # Calculate the overlap between all list elements
  overlap_matrix <- matrix(0, nrow = length(cluster_GO_terms), ncol = length(cluster_GO_terms))
  for (i in 1:length(cluster_GO_terms)) {
    for (j in 1:length(cluster_GO_terms)) {
      overlap_matrix[i, j] <- calculate_overlap_percentage(cluster_GO_terms[[i]], cluster_GO_terms[[j]])
    }
  }

                                        # Set row and column names for the heatmap
  row_names <- names(cluster_GO_terms)
  col_names <- names(cluster_GO_terms)
  colnames(overlap_matrix) <- col_names
  rownames(overlap_matrix) <- row_names


  col_fun <- circlize::colorRamp2(c(0,80, 100), c("white", grattan_yellow, grattan_red))
  cluster_overlap_heatmap <- Heatmap(overlap_matrix, col = col_fun, border = T,
                                     row_names_gp = gpar(fontsize = 8),
                                     column_names_gp = gpar(fontsize = 8),
                                     heatmap_legend_param = list(
                                       title = "Percent overlap",
                                       at = c(0,50,80,100),
                                       labels = c("0%", "50%","80%", "100%")
                                     ),
                                     row_km = 5,
                                     column_km = 5)
  cluster_overlap_heatmap

  # Filter the overlap matrix to keep only values above 80% threshold
  overlap_threshold <- 80
  overlap_subset <- overlap_matrix > overlap_threshold


                                     ))

  return(list(full = cluster_overlap_heatmap,
              subset = cluster_subset_heatmap))



  }

get_subcluster_names <- function(overlap_matrix, row_clusters) {
  overlap_matrix[row_clusters,] |>rownames()
}

Plot_cluster_GOchart <- function(cluster_GO_terms, clusters, abundance_threshold = 1){
  data <- reshape2::melt(cluster_GO_terms) |>
    rename(GO = value, Cluster = L1) |>
    filter(Cluster %in% clusters)

  keep <- names(which(table(data$GO) > abundance_threshold))
  keep <- keep[nzchar(keep)]

  data <- data |>filter(GO %in% keep)

  plt <- data |>
    ggplot() +
    geom_bar(aes(x = GO)) + theme_grattan() +
    theme(axis.text.x = element_text(angle=90))

  return(plt)
}

set.seed(seed)
cluster_overlap_heatmap <- Heatmap(overlap_matrix, col = col_fun, border = T,
                                     row_names_gp = gpar(fontsize = 8),
                                     column_names_gp = gpar(fontsize = 8),
                                     heatmap_legend_param = list(
                                       title = "Percent overlap",
                                       at = c(0,50,80,100),
                                       labels = c("0%", "50%","80%", "100%")
                                     ),
                                     row_km = 5,
                                   column_km = 5)
cluster_overlap_heatmap

set.seed(10)
row_clusters = row_order(cluster_overlap_heatmap)

set.seed(10)
col_clusters = column_order(cluster_overlap_heatmap)



plt5<-Plot_cluster_GOchart(cluster_GO_terms, get_subcluster_names(overlap_matrix, row_clusters[['5']]))
plt1<-Plot_cluster_GOchart(cluster_GO_terms, get_subcluster_names(overlap_matrix, row_clusters[['1']]),
                           abundance_threshold = 0)
