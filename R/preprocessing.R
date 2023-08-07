library(tidyverse)
                                        # Data Import ------------------------
setwd("~/scratch/folding_comparison")



comparisons <- read_tsv(
  "raw_data/folding_fatcat_output_inc_rb1_globalidentities.tsv"
) |> rbind(
       read_tsv(
         "raw_data/folding_fatcat_output_rbp_vs_all_globalidentities.tsv"
       )
     )

###



annotations_pal <- read.csv2(
  "raw_data/pallida_annotations.csv",
  sep = "\t"
)

rownames(annotations_pal) <- annotations_pal$SeqName

gpal_ids_clean <- str_replace_all(comparisons$Subject, pattern = '.pdb', replacement = "")

comparisons$Subject_Description <- annotations_pal[gpal_ids_clean, "Description"]
comparisons$Subject_GO_ID <- annotations_pal[gpal_ids_clean, "GO.IDs"]
comparisons$Subject_Enzyme_Names <- annotations_pal[gpal_ids_clean, "Enzyme.Names"]
comparisons$Subject_InterPro_GO_Names <- annotations_pal[gpal_ids_clean, "InterPro.GO.Names"]
comparisons$Subject_InterPro_IDs <- annotations_pal[gpal_ids_clean, "InterPro.IDs"]


# Mchit functional annotations ----
annotations_mchit <- read.csv2(
  "raw_data/chitwoodi_annotations.csv",
  sep = '\t'
  )

rownames(annotations_mchit) <- annotations_mchit$SeqName

mc_ids_clean <- str_replace_all(comparisons$Query, pattern = '.pdb', replacement = "")

comparisons$Query_Description <- annotations_mchit[mc_ids_clean, "Description"]
comparisons$Query_GO_ID <- annotations_mchit[mc_ids_clean, "GO.IDs"]
comparisons$Query_Enzyme_Names <- annotations_mchit[mc_ids_clean, "Enzyme.Names"]
comparisons$Query_InterPro_GO_Names <- annotations_mchit[mc_ids_clean, "InterPro.GO.Names"]
comparisons$Query_InterPro_IDs <- annotations_mchit[mc_ids_clean, "InterPro.IDs"]


# Pvalue adjustment ----------------


log_p_zero_replacement <- 0.0000000000001

adj_factor <- comparisons$Subject %>%
  unique() %>% length

comparisons$P.adjusted <- as.numeric(comparisons$`P-value`) * adj_factor
comparisons$P.adjusted[comparisons$P.adjusted > 1] <- 1
comparisons$P.adjustedLog <- -log10(comparisons$P.adjusted)
comparisons$P.adjustedLog[is.infinite(comparisons$P.adjustedLog)] <- log_p_zero_replacement

comparisons.filtered <- comparisons %>%
  filter(P.adjusted <=0.05)
dim(comparisons.filtered)

####Test FDR
comparisons.fdr <- comparisons
comparisons.fdr$P.adjusted <- as.numeric(comparisons$`P-value`) |> p.adjust(method="bonferroni")

comparisons.fdr.filtered <- comparisons.fdr |>
  filter(P.adjusted <= 0.05)
dim(comparisons.fdr.filtered)



#comparisons$P.adjusted <- p.adjust.custom(
#  comparisons,
#  method = "bonferroni"
                                        #)



node_table <- read_csv("raw_data/cytoscape_node_table.csv") |>
  filter(species=="Mchit") |>
  select(c("__fastGreedyCluster",'id')) |>
  setNames(c("Cluster", "id"))


#load("processed_data/2023-08-01_folding_fatcat_output_all_pvalueadjusted.RData")


comparisons.filtered <- dplyr::left_join(comparisons.filtered, node_table, join_by(Query == id) )


                                        # Transmembrane filter --------------------

cluster_prune <- comparisons.filtered |>
  filter( grepl("membrane", Subject_Description) | grepl("membrane", Query_Description)) |>
  select(Cluster) |>
  unlist() |>
  unique()

comparisons.filtered %<>% filter(!Cluster %in% cluster_prune)

save(comparisons, comparisons.filtered,
    file = paste0("processed_data/",
                  paste(
                    Sys.Date(),
                    "folding_fatcat_output_all_pvalueadjusted_inc_clusters_membranePrune.RData",
                    sep = '_' )
                  )
    )

paste0(
  "processed_data/",
  paste(Sys.Date(), "folding_fatcat_output_all_pvalueadjusted.RData",sep="_")
) %>%
  save(comparisons, comparisons.filtered, file = .)
