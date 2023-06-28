# Load required libraries
library(dplyr)

# Read the edge table
edges <- comparisons.filtered

# Read the cluster assignment table
clusters <- read.csv("cytoscape_clusters_new_new.csv")


# Merge the edge and cluster tables based on source and target IDs

merged_data <- merge(edges, clusters, by.x = "Query", by.y = "name", all.x = TRUE)
merged_data <- merge(merged_data, clusters, by.x = "Subject", by.y = "name", all.x = TRUE)



# Calculate the average Identity score per cluster
avg_identity <- merged_data %>%
  group_by(X__fastGreedyCluster.x) %>%
  summarise(avg_Identity = mean(global_identity, na.rm = TRUE))

# Print the result
avg_identity$avg_Identity |> mean()

merged_data$global_identity |>mean()
