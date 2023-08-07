pacman::p_load(tidyverse, magrittr,usethis, testthat,
               here, gridExtra,gridGraphics, gplots,
               colorRamp2, RColorBrewer, reshape2, ComplexUpset,
               figpatch, ComplexHeatmap, patchwork)


#require(RCy3)
library(grattantheme)

setwd("~/scratch/folding_comparison/")

source("R/histogram_hits_per_query.R")
source("R/cluster_overlaps.R")

load("processed_data//2023-08-01_folding_fatcat_output_all_pvalueadjusted_inc_clusters.RData")



                                        # Histogram number of hits per query

histograms <- Histogram_hit_per_query()






structure_plot.df <- comparisons |>
  select('Query','Subject','P.adjustedLog') |>
  pivot_wider(names_from = 'Subject', values_from = `P.adjustedLog`) |>
  dplyr::select(-"Query") |>
  as.matrix() |>
  apply(MARGIN=2, as.numeric)
colscale <- color.scale(colors.use = brewer.pal(11, "Reds"), y.val = as.numeric(structure_plot.df))
structure_map <- grid.grabExpr(draw(ComplexHeatmap::Heatmap(structure_plot.df,
                        name = "log10P",
                        col = colscale,
                        row_names_side = NULL,
                        column_names_side = NULL,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        column_title = "G. pallida structures",
                        row_title = "M. chitwoodi structures",
                        heatmap_legend_param = list(
                          title = "-log10(P)",
                          at = c(0,1.3, 5, 10,15)
                        )
                        )
                        )
                        )


## png('Figures/subfigures/heatmap_structures.png', 2500, 2800, pointsize = 100)
## #svg('subfigures/heatmap_structures.svg')
## colscale <- color.scale(colors.use = brewer.pal(11, "Reds"), y.val = as.numeric(structure_plot.df))
## heatmap(structure_plot.df, scale = 'none', col=colscale, margins=c(10,8))
## dev.off()
## png("subfigures/heatmap_structures_scale.png")
## scale.legend(input=as.numeric(structure_plot.df),start.zero=F,col.scale=colscale,line.col="black")
## dev.off()



sequence_plot.df <- comparisons |>
  select('Query','Subject','global_identity') |>
  pivot_wider(names_from = 'Subject', values_from = 'global_identity') |>
  dplyr::select(-'Query') |>
  as.matrix() |>
  apply(MARGIN=2, as.numeric)
colscale <- color.scale(colors.use = brewer.pal(11, "Reds"), y.val = as.numeric(sequence_plot.df))
sequence_map <- grid.grabExpr(draw(ComplexHeatmap::Heatmap(sequence_plot.df,
                        name = "% identity",
                        col = colscale,
                        row_names_side = NULL,
                        column_names_side = NULL,
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        column_title = "G. pallida sequences",
                        row_title = "M. chitwoodi sequences",
                        heatmap_legend_param = list(
                          title = "% identity",
                          at = c(0,15, 30, 45, 60, 75))
                        )
                      )
                    )




                                        # Vulcano Plot -----

plot.df <- cbind(
  as.numeric(comparisons$global_identity),
  comparisons$P.adjustedLog
)

plot.df%<>% as.data.frame()
colnames(plot.df) <- c("Identity","Pval")

pvalue_log <- -1 * log10(0.05)
volcano_plt <- plot.df %>%
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
    values = setNames(c(grattan_orange, grattan_darkred),
                      c(TRUE, FALSE)
                      )
    ) +
  xlab("Percent Identity") +
  ylab("-log10(pvalue)") +
  geom_hline(yintercept=pvalue_log, linetype='dashed', col = 'red')+
  scale_x_continuous(limits = c(0,100))+
  ggtitle("") +
  theme_grattan(chart_type = 'scatter')




design <- "
  12233
  4455#
  4455#
"
volcano_plt + histograms$gpal_hist + histograms$mchit_hist +
  wrap_elements(sequence_map) + wrap_elements(structure_map) +
  plot_layout(design = design) & plot_annotation(tag_levels = 'A')

ggsave(file="Figures/Panel_comparisons.png", width = 5000, height = 3000, unit='px')
ggsave(file="Figures/Panel_comparisons.pdf", width = 5000, height = 3000, unit='px')




                                        ##### Cluster statistics
data <- comparisons.filtered |>
  group_by(Cluster) |>
  summarise(
    Mchit = n_distinct(Query),
    Gpal = n_distinct(Subject)
  ) |>
  gather(key = "variable", value = "value", -Cluster)

stacked_proportion_plt <- data |>
  group_by(Cluster) |>
  mutate(proportion = value / sum(value)) |>
  ggplot() +
  geom_bar(aes(x = Cluster, y = proportion, fill = variable), stat = 'identity', position = 'stack') +
  theme_grattan(legend = 'bottom') +
  labs(title = "Relative Cluster Contribution per Species", y = "Relative Contribution")


total_plt <- data |>
  group_by(Cluster) |>
  summarise(Total = sum(value)) |>
  ggplot() +
  geom_bar(aes(x = Cluster, y = Total), fill = grattan_blue, stat = 'identity') +
  theme_grattan() +
  #theme(axis.title.y = element_text(angle=90)) +
  labs(title = "Number of Proteins per Cluster", y = "No. Proteins")



mchit_cl_data <- data |>
  group_by(Cluster) |>
  mutate(proportion = value / sum(value)) |>
  filter(variable == "Mchit")

Gpal_cl_data <- data |>
  group_by(Cluster) |>
  mutate(proportion = value / sum(value)) |>
  filter(variable == "Gpal")

n <- length(mchit_cl_data$proportion)
sum_recip_mchit <- sum(1/mchit_cl_data$proportion)
sum_recip_gpal <- sum(1/Gpal_cl_data$proportion)

harmonic_mchit <- n / sum_recip_mchit
harmonic_gpal <- n / sum_recip_gpal
harmonic_mchit
harmonic_gpal


cluster_heatmaps <- Cluster_overlap_heatmap(comparisons.filtered)

design <- "
  13
  23
"
total_plt + stacked_proportion_plt + wrap_elements(grid.grabExpr(draw(cluster_overlap_heatmap))) +
                                      plot_layout(design = design) + plot_annotation(tag_levels = "A")


ggsave(file="Figures/Cluster_stats.png")
ggsave(file="Figures/Cluster_stats.pdf")


######################

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


#
