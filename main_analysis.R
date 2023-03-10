pacman::p_load(tidyverse, RCy3, magrittr,usethis, testthat, here)
Sys.setenv(RETICULATE_PYTHON = "/home/joris/miniconda3/bin/python")

library(foldingRproject)
                                        # Function definitions ----




# Data Import ------------------------
setwd("~/scratch/folding_comparision/")

comparisons <- read_tsv(
  "folding_fatcat_output_inc_rb1_globalidentities.tsv"
) |> rbind(
       read_tsv(
         "folding_fatcat_output_rbp_vs_all_globalidentities.tsv"
       )
     )




# Pvalue adjustment ----------------



adj_factor <- comparisons$Subject %>%
  unique() %>% length
comparisons$P.adjusted <- as.numeric(comparisons$`P-value`) * adj_factor

comparisons.filtered <- comparisons %>%
  filter(P.adjusted <=0.05)

#comparisons$P.adjusted <- p.adjust.custom(
#  comparisons,
#  method = "bonferroni"
#)

write.csv2(comparisons, "folding_fatcat_output_all_pvalueadjusted.csv")

                                        # Overview 3m comparisons
counts_per_query <- comparisons.filtered |>
  group_by(Query) |>
  summarise(count = n())

png(filename="histogram.png")
 counts_per_query |>
   ggplot() +
   geom_histogram(
     aes( x = count),
     position = position_dodge(width = 0.1),
     binwidth = 1) +
   scale_x_continuous(breaks=c(seq(0,30,5),40,50,60,70,80,90,100))+
   xlab("Number of hits per query") +
   ylab("Count") +
   ggtitle("Number of hits per query") +
   theme_minimal()
dev.off()


comparisons_matrix <- comparisons |>
  select('Query','Subject','P-value') |>
  pivot_wider(names_from = 'Subject', values_from = `P-value`) |>
  as.matrix()

rownames(comparisons_matrix) <- comparisons_matrix[,1]

comparisons_matrix <- comparisons_matrix[,-1] |>
  apply(MARGIN=2, as.numeric)

reorder_matrix <- function(matrix, cutoff = 0.05) {
  matrix <- matrix[,(matrix < cutoff )|>
                    colSums() |>
                    order(decreasing = T)]
  matrix <- matrix[(matrix < cutoff )|>
                    rowSums() |>
                    order(decreasing = T),]
  return(matrix)
}

png(filename="~/org-roam/alphafold_data/heatmap.png")
  comparisons_matrix |>
    reorder_matrix(cutoff=0.05) |>
    apply(MARGIN=2, function(x) as.numeric(x < 0.05)) |>
    gplots::heatmap.2(trace = "none",
                    labRow = FALSE,
                    labCol = FALSE,
                    dendrogram = 'none',
                    col = c("#000000","#00FF00"),
                    Rowv=FALSE,
                    Colv=FALSE,
                    key=FALSE,
                    xlab = 'G. pallida proteins',
                    ylab = 'M. chitwoodi proteins')
dev.off()

png(filename="~/org-roam/alphafold_data/supp_pairwise_identity.png")
comparisons.filtered$
  global_identity |>
  density()|>
  plot(main="Supplemental - Pairwise Sequence Identities")|>
  abline(
    v = mean(comparisons.filtered$global_identity)
  ) |>
  abline(
    v=
      mean(comparisons.filtered$global_identity) +
      0.5*sd(comparisons.filtered$global_identity)
  ) |>
  abline(
    v=
      mean(comparisons.filtered$global_identity) -
      0.5*sd(comparisons.filtered$global_identity)
  )
  text(x = mean(comparisons.filtered$global_identity) , labels = "Mean")
dev.off()



                                        # Vulcano Plot -----

plot.df <- cbind(
  as.numeric(comparisons$global_identity),
  -1 * log10(as.numeric(comparisons$P.adjusted + 0.00000000001))
)

plot.df%<>% as.data.frame()
colnames(plot.df) <- c("Identity","Pval")

pvalue_log <- -1 * log10(0.05)
plot.df %>%
  ggplot +
  geom_point(
    aes(
      x = Identity,
      y = Pval,
      colour = Pval >= pvalue_log
      )
    ) +
  scale_colour_manual(
    name = "Significant",
    values = setNames(c("blue", "black"),
                      c(TRUE, FALSE)
                      )
    ) +
  xlab("Percent Identity") +
  ylab("-Log10(pvalue)") +
  ggtitle("Identity vs P-value") +
  theme_bw()

############################################
## ggsave(                                ##
##   filename = 'identity_vs_pvalue.pdf', ##
##   device = 'pdf')                      ##
## ggsave(                                ##
##   filename = 'identity_vs_pvalue.png', ##
##   device = 'png')                      ##
############################################



                                        # How many N is in each group -----
 ident_threshold <- 30
low_ident_noSig <- comparisons %>%
  filter(P.adjusted >= 0.05 & global_identity < ident_threshold)
low_ident_Sig <- comparisons %>%
  filter(P.adjusted < 0.05 & global_identity < ident_threshold)
h_ident_noSig <- comparisons %>%
  filter(`P-value` >= 0.05 & global_identity >= ident_threshold)
h_ident_Sig <- comparisons %>%
  filter(P.adjusted < 0.05 & global_identity >= ident_threshold)

low_ident_noSig %>% nrow
low_ident_Sig %>% nrow
h_ident_noSig%>% nrow
h_ident_Sig%>% nrow



                                        # Gpal functional annotations ----

library(stringr)
annotations_pal <- read.csv2(
  "/mnt/nemahomes/steen176/Genomes/Annotations/G_pallida_annotations_onlyt1pruned.csv",
  sep = "\t")

rownames(annotations_pal) <- annotations_pal$SeqName

gpal_nodes <- cbind(
  comparisons.filtered$Subject,
  rep(
    "Gpal", times = length(comparisons.filtered$Subject)
    )
  )
colnames(gpal_nodes) <- c('id','species')

gpal_ids_clean <- str_replace_all(gpal_nodes[,'id'], pattern = '.pdb', replacement = "")

gpal_nodes %<>% cbind(.,
                      annotations_pal[
                        gpal_ids_clean,
                        c(
                          "Description",
                          "GO.IDs",
                          "Enzyme.Names",
                          "InterPro.GO.Names",
                          "InterPro.IDs"
                          )
                        ]
                      )

# Mchit functional annotations ----
annotations_mchit <- read.csv2(
  "/mnt/dataHDD/scratch/chitwoodi_genomes/annotation/omicsbox_table_onlyt1.txt",
  sep = '\t'
  )

rownames(annotations_mchit) <- annotations_mchit$SeqName

mchit_nodes <- cbind(
  comparisons.filtered$Query,
  rep(
    "Mchit", times = length(comparisons.filtered$Query)
    )
)

colnames(mchit_nodes) <- c('id','species')

mchit_ids_clean <- str_replace_all(mchit_nodes[,'id'], pattern = '.pdb', replacement = "")

mchit_nodes %<>% cbind(.,
                      annotations_mchit[
                        mchit_ids_clean,
                        c(
                          "Description",
                          "GO.IDs",
                          "Enzyme.Names",
                          "InterPro.GO.Names",
                          "InterPro.IDs"
                          )
                        ]
                      )
                                        # Format Node and Edge data frames ----
nodes <- rbind(
  gpal_nodes,
  mchit_nodes
  ) %>%
  data.frame

edges <- comparisons.filtered[,c('Query',"Subject", "Identity",'P.adjusted')] %>% data.frame
colnames(edges) <- c("source","target", "Identity",'P.adjusted')

write.csv2(edges, file = 'cytoscape_edges_r.csv')

                                        # Export Cytoscape
createNetworkFromDataFrames(nodes,edges)

style.name = "myStyle"
defaults <- list(NODE_SHAPE="circle",
                 NODE_SIZE=30,
                 EDGE_TRANSPARENCY=30,
                 NODE_LABEL_POSITION="W,E,c,0.00,0.00")

nodeLabels <- mapVisualProperty('node label','id','p')
nodeFills <- mapVisualProperty('node fill color','species','d',c("Gpal","Mchit"), c("#FF9900","#66AAAA"))
edgeWidth <- mapVisualProperty('edge width','Identity','p')


createVisualStyle(style.name, defaults, list(nodeFills,edgeWidth))
setVisualStyle(style.name)


                                        # Cluster Degrees ------
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
cluster_sizes = cluster_table$X__fastGreedyCluster %>% table

cluster_sizes[which(cluster_sizes != 2)] %>% mean


## 1-1 cluster names
cl_names <- cluster_sizes[which(cluster_sizes == 2)] %>%
  names %>%
  sapply(function(x)
    paste0('Cl', x)) %>%
  as.character

cluster_descriptions %>%
  filter(Cluster %in% cl_names) %>%
  View


