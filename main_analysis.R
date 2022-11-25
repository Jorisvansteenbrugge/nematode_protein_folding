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

comparisons$P.adjusted <- p.adjust.custom(comparisons)

comparisons.filtered <- comparisons %>%
  filter(P.adjusted <= 0.05)


                                        # Plot identisy all vs significant only

all_identity <- cbind(
  as.numeric(comparisons$global_identity),
  rep("All", times = nrow(comparisons)))

sig_identity <- cbind(
  comparisons.filtered$global_identity,
  rep("Significant only", times = nrow(comparisons.filtered))
)

identities <- rbind(all_identity, sig_identity) %>% data.frame
colnames(identities) <- c("Identity", "Set")
Identities$Identity %<>% as.numeric

all_rmsd <- cbind(
  as.numeric(comparisons$opt.rmsd),
  rep("All", times = nrow(comparisons))
)
sig_rmsd <- cbind(
  comparisons.filtered$opt.rmsd,
  rep("Significant only", times = nrow(comparisons.filtered)))

RMSDs <- rbind(all_rmsd, sig_rmsd) %>% data.frame(stringsAsFactors = FALSE)
colnames(RMSDs) <- c("RMSD", "Set")
RMSDs$RMSD %<>% as.numeric

Identities %>%
  ggplot +
  geom_density(
    aes(
      x = Identity,
      fill = Set
      ),
    alpha = 0.5
    ) +
  scale_x_continuous(
    name = "Percent Sequence Identity",
    breaks = c(
      seq(1, 30, 2),
      seq(30, 50, 5), 80)
  ) +
  ggtitle("Sequence Identity")

                                        # Vulcano Plot -----

plot.df <- cbind(
  as.numeric(comparison),
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
  filter(P.value >= 0.05 & global_identity >= ident_threshold)
h_ident_Sig <- comparisons %>%
  filter(P.adjusted < 0.05 & global_identity >= ident_threshold)

low_ident_noSig %>% nrow
low_ident_Sig %>% nrow
h_ident_noSig%>% nrow
h_ident_Sig%>% nrow


                                        # Playground ----


Lookup_effector_family <- function(fam_name) {
 effector_ids <- effector_list |>
    filter(Family == fam_name)


  return(comparisons |>
    filter(Subject %in% effector_ids$ID)
    )
}

vap_rmsd <- c(0.526, 0.313, 0.344, 0.816, 0.673, 0.644, 0.765)
vap_rmsd|>min()
vap_rmsd|>max()
vap_rmsd |>mean()

sprysec_b_genes <- c("Gpal_D383_g12450.pdb", "Gpal_D383_g12451.pdb","Gpal_D383_g12443.pdb", "Gpal_D383_g12448.pdb", "Gpal_D383_g12403.pdb",
                   #  "Gpal_D383_g12854.t4",
                     "Gpal_D383_g16728.pdb", "Gpal_D383_g15175.pdb", "Gpal_D383_g16477.pdb")
sprysec_rmsd <- c(3.017, 1.579, 0.966, 1.732, 1.249, 1.883, 1.033, 1.367)
sprysec_rmsd |> min()
sprysec_rmsd |> max()
sprysec_rmsd |> mean()



effector_list <- read.csv(
  "All_pallida_named_effectos.txt",
  sep = '\t',
  header = FALSE,
  col.names = c("ID", "Family")
)

Lookup_effector_family("VAP_CDS")


plot(
  subset$Identity,
  subset$opt.rmsd,
  main="Identity vs RMSD",
  xlab="Percent Identity",
  ylab="Optimised RMSD")


subset$opt.rmsd %>% density %>%  plot(main="Optomised RMSD", xlab='RMSD')

sprysec_example = "Gpal_D383_g01418.pdb"
sprysec_example2 = "Gpal_D383_g12494.pdb"
subset_sprysec <- comparisons.filtered %>%
  filter(Subject==sprysec_example2)


subset_sprysec$opt.rmsd %>% density %>%  plot(main="Optomised RMSD", xlab='RMSD')

                                        # Heatmap ----

library(gplots)
pvalue_matrix <- comparisons %>%
  .[,c("Query","Subject","P.adjusted")] %>%
  reshape(
    direction = 'wide',
    idvar="Query",
    timevar = "Subject"
    )

rownames(pvalue_matrix) <- pvalue_matrix[,1]
pvalue_matrix[,1] <- NULL

d <- pvalue_matrix %>%
  as.matrix %>%
  add(0.0000001) %>%
  log10 %>%
  multiply_by(-1)

heatmap.2(
  d,
  trace = "none"
  )


                                        # Gpal functional annotations ----
reticulate::source_python("misc.py")
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

gpal_ids_clean <- strip_pdb_ids(gpal_nodes[,'id'])

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

mchit_ids_clean <- strip_pdb_ids(mchit_nodes[,'id'])

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
edges$P.adjustedRev <- 1-edges$P.adjusted
edges$P.adjustedLog <- -1*log(edges$P.adjusted+.Machine$double.xmin)
write.csv2(edges, file = 'cytoscape_edges.csv')


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


                                        # merge cluster number and Identity value
reticulate::source_python('parse_cluster_info.py')
cluster_identities <- get_cluster_identities()

cluster_identities$Identity %<>%
  sapply(function(x) sub(pattern = ",", replacement = ".",x = x)) %>%   as.numeric()
cluster_identities$P.adjusted %<>%
  sapply(function(x) sub(pattern = ",", replacement = ".",x = x)) %>%   as.numeric()



cluster_sizes <- cluster_identities$Cluster %>% table
cluster_sizes_pruned <- cluster_sizes %>%
  magrittr::is_less_than(2) %>%
  `!` %>%
  which %>%
  names

cluster_identities_pruned <- cluster_identities %>%
  filter(Cluster %in% cluster_sizes_pruned)

#rownames(cluster_identities_pruned) <- cluster_identities_pruned$Cluster


                                        # Plot Identity vs pvalue -> per cluster
cluster_identities_pruned %>%
  ggplot +
  geom_point(aes(x=Identity, y = -1*log10(P.adjusted+0.000001)))  +
  facet_wrap(~Cluster)

clusterNames <- unique(cluster_identities$Cluster)
cluster_table <- read.csv2('cytoscape_clusters_new.csv', sep=',')


                                        # Cluster Average Identity ----
                                        # to do

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


cluster_descriptions <- sapply(clusterNames, function(cluster){
  cl_clean <- sub('Cl','',cluster)
  avg_cluster_identity <- cluster_identities %>%
    filter(Cluster == cluster) %>%
    .$Identity %>%
    mean

  description_gpal <- get_most_occuring('Description', 'Gpal',cl_clean)
  GOnames_gpal <- get_most_occuring('InterPro.GO.Names', 'Gpal',cl_clean)

  description_mchit <- get_most_occuring('Description', 'Mchit',cl_clean)

  GOnames_mchit <- get_most_occuring('InterPro.GO.Names', 'Mchit',cl_clean)

  return(
    c(
      cluster,
      avg_cluster_identity,
      description_gpal,
      description_mchit,
      GOnames_gpal,
      GOnames_mchit
      )


  )
})  %>% do.call(rbind.data.frame, .)

colnames(cluster_descriptions) <- c(
  "Cluster",
  "AvgClusterIdentity",
  "gpalDescription",
  "mchitDescription",
  "gpalGOnames",
  "mchitGOnames"
  )
cluster_descriptions$AvgClusterIdentity %<>% as.numeric
View(cluster_descriptions)


exp_table <- read.csv2("./blasting/within_cluster_stats2.0.csv", sep=';')

get_perc_translated <- function(spec_a, spec_b){

  average_perc <- sapply(1:125, function (i){

    cl_id <- paste0("Cl",i)
    current_spec_a <- exp_table |>
      filter(Cluster == cl_id &
             Exp_species == spec_a) |>
      select(Exp_value) |>
      as.numeric()
    current_spec_b <- exp_table |>
      filter(Cluster == cl_id &
             Exp_species == spec_b) |>
      select(Exp_value) |>
      as.numeric()

    translate_percent <- current_spec_b/current_spec_a*100
    if (translate_percent > 100) translate_percent <- 100
      return(translate_percent)
  })

  #print(average_perc)

  return(mean(average_perc))
}


get_perc_translated("G. pallida", "G. rostochiensis")
get_perc_translated("G. pallida", "H. schachtii")


get_perc_translated("M. chitwoodi", "M. hapla")
get_perc_translated("M. chitwoodi", "M. incognita")
get_perc_translated("M. chitwoodi", "M. graminicola")


get_perc_translated("G. pallida", "M. chitwoodi")

exp_table |>
  ggplot() +
  geom_bar(
    aes(
      x=Species,
      y=log2(Exp_value),
      fill=Exp_species),
    stat='identity'
  )  +
  facet_wrap(~Cluster)
