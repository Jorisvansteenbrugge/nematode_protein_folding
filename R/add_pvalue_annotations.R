pacman::p_load(tidyverse, magrittr, kableExtra, magick, webshot)



setwd("/home/joris/scratch/folding_comparison/knowable_genome/")

annotation_table <- read_tsv(
  "knowable_cluster_table_annotated_filtered.csv"
)

annotation_table_unfiltered <- read_tsv(
  "knowable_cluster_table_annotated.csv"
)

network_table <- read_csv2(
  "../folding_fatcat_output_all_pvalueadjusted.csv"
)

network_table$Query   %<>% stringr::str_replace_all(pattern = '.pdb', replacement = '.t1')
network_table$Subject %<>% stringr::str_replace_all(pattern = '.pdb', replacement = '.t1')



annotation_table_pvalue_gpal <- annotation_table |>
  filter(Species=="Gpal") |>
  left_join(
    network_table,
    by = join_by(
      Unannotated_gene == Subject,
      Target ==Query ) ) |>
  select(c("Cluster","Species", "Unannotated_gene", "Target", "P-value", "Annotation", "Pfam Domains"))

annotation_table_pvalue_mchit <- annotation_table |>
  filter(Species=="Mchit") |>
  left_join(
    network_table,
    by = join_by(
      Unannotated_gene == Query,
      Target == Subject ) ) |>
  select(c("Cluster","Species", "Unannotated_gene", "Target", "P-value", "Annotation", "Pfam Domains"))

annotation_table_pvalue |>
  arrange(Cluster, Species, Unannotated_gene, `P-value`) |>
  kable(format='html', booktabs = T) |>
  cat(file="~/org-roam/alphafold_data/knowable_table.html")

annotation_table_pvalue_gpal$Unannotated_gene |> unique() |> length()
annotation_table_pvalue_mchit$Unannotated_gene |> unique() |> length()


                                        # Number of annotatable genes
annotation_table_pvalue_mchit |>
  group_by(Unannotated_gene) |>
  summarise(mean=mean(`P-value`,na.rm = T)) |>
  filter(!is.nan(mean)) |>
  nrow()

annotation_table_pvalue |>
  filter(Species=='Gpal') |>
  group_by(Unannotated_gene) |>
  summarise(mean=mean(`P-value`,na.rm = T)) |>
  filter(!is.nan(mean)) |>
  nrow()


                                        # Gene ids that remain unannotated
annotation_table_pvalue_mchit |>
  group_by(Unannotated_gene) |>
  summarise(mean=mean(`P-value`,na.rm = T)) |>
  filter(is.nan(mean)) |>
  select(Unannotated_gene) |>
  write.csv2(file = 'mchit_unannotated_gene_ids.csv')

annotation_table_pvalue_gpal |>
  group_by(Unannotated_gene) |>
  summarise(mean=mean(`P-value`,na.rm = T)) |>
  filter(is.nan(mean)) |>
  select(Unannotated_gene) |>
  write.csv2(file = 'gpal_unannotated_gene_ids.csv')
