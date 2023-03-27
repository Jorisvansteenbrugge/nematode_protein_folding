
library(dplyr)
library(ggplot2)
data <- read.csv2(
  "/home/joris/scratch/folding_comparision/knowable_genome/knowable_cluster_table_annotated_filtered.csv",
  sep = '\t')

gene_c <- "Gpal_D383_g06418.t1"
annotations_c <- data |> filter(Unannotated_gene == gene_c)

freq_table <- annotations_c$Annotation[1:10] |> strsplit(split=';') |>unlist() |>table() |>as.data.frame()
colnames(freq_table) <- c("word", "freq")


freq_table |>
  ggplot(aes(x=word, y = freq)) +
  geom_point() +
  geom_segment( aes(x = word, xend=word , y = 0, yend =freq)) +
  theme(axis.text.x = element_text(angle=90)) +
  coord_flip()
