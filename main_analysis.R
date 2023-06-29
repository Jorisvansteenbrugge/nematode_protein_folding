pacman::p_load(tidyverse, magrittr,usethis, testthat, here, gridExtra,gridGraphics, gplots, colorRamp2, RColorBrewer)
require(RCy3)



# Data Import ------------------------
setwd("~/scratch/folding_comparision/")

comparisons <- read_tsv(
  "folding_fatcat_output_inc_rb1_globalidentities.tsv"
) |> rbind(
       read_tsv(
         "folding_fatcat_output_rbp_vs_all_globalidentities.tsv"
       )
     )

###


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

#comparisons$P.adjusted <- p.adjust.custom(
#  comparisons,
#  method = "bonferroni"
#)

#write.csv2(comparisons, "folding_fatcat_output_all_pvalueadjusted.csv")


                                        # Histogram number of hits per query
png("subfigures/gpal_hist.png")
comparisons.filtered |>
  group_by(Query) |>
  summarise(count = n()) |>
   ggplot() +
   geom_histogram(
     aes( x = count),
     position = position_dodge(width = 0.1),
     binwidth = 1) +
  scale_x_continuous(breaks=c(seq(0,30,10),40,50,60,70,80,90,100))+
  coord_cartesian(ylim = c(0,150),
                  xlim = c(0,100)) +
   xlab("Number of hits per query") +
   ylab("Count") +
   ggtitle("") +
    theme_minimal()
dev.off()
png("subfigures/mchit_hist.png")
comparisons.filtered |>
  group_by(Subject) |>
  summarise(count = n()) |>
   ggplot() +
   geom_histogram(
     aes( x = count),
     position = position_dodge(width = 0.1),
     binwidth = 1) +
  scale_x_continuous(breaks=c(seq(0,30,10),40,50,60,70,80,90,100))+
  coord_cartesian(ylim = c(0,150),
                  xlim = c(0,100))+
   xlab("Number of hits per query") +
   ylab("Count") +
   ggtitle("") +
   theme_minimal()
dev.off()



reorder_matrix <- function(matrix, cutoff = 0.05, decreasing = T) {
  matrix <- matrix[,(matrix < cutoff )|>
                    colSums() |>
                    order(decreasing = decreasing)]
  matrix <- matrix[(matrix < cutoff )|>
                    rowSums() |>
                    order(decreasing = decreasing),]
  return(matrix)
}







structure_plot.df <- comparisons |>
  select('Query','Subject','P.adjustedLog') |>
  pivot_wider(names_from = 'Subject', values_from = `P.adjustedLog`) |>
  dplyr::select(-"Query") |>
  as.matrix() |>
  apply(MARGIN=2, as.numeric)

## # Define the breakpoints (-5, 1.3, 11)
## breakpoints <- c(0, 1.3, max(structure_plot.df, na.rm=T))

## # Create a custom color palette
## my_palette <- colorRampPalette(c("blue", "yellow"))(n = 10)

## # Replace colors for values below 1.3 with shades of blue

## my_palette[1:(breakpoints[2] - breakpoints[1]) / (breakpoints[3] - breakpoints[1]) * 10] <- colorRampPalette(c("lightblue", "blue"))(n = (breakpoints[2] - breakpoints[1]) / (breakpoints[3] - breakpoints[1]) * 10 )


## gplots::heatmap.2(
##           structure_plot.df,
##           trace = "none",
##           labRow = FALSE,
##           sepcolor = NA,
##           labCol = FALSE,
##           dendrogram = 'both',
##           col = my_palette,
##           #Rowv=FALSE,
##           #Colv=FALSE,
##           key=TRUE,
##           xlab = 'G. pallida Structures',
##           ylab = 'M. chitwoodi Structures')

color.scale <- function(colors.use,y.val){
                        if(class(y.val) == "matrix"){
                            scale.use <-  c(as.numeric(rownames(y.val))[nrow(y.val)],as.numeric(rownames(y.val))[1],as.numeric(rownames(y.val))[1]-as.numeric(rownames(y.val))[2])
                        }
                        if(class(y.val) == "numeric"){
                            scale.use <- c(min(y.val,na.rm=T),max(y.val,na.rm=T),diff(c(min(y.val,na.rm=T),max(y.val,na.rm=T)))/100)
                        }
                        seq1 <- seq(scale.use[1],scale.use[2],scale.use[3])    ###adjust for optimal effect
                        ramp <- colorRamp(colors.use, bias=1)
                        colors.out <- rgb(ramp(seq(0, 1, length = length(seq1))), max = 255)
                        return(colors.out)
                       }

scale.legend <- function(input,bg.col,col.scale,x.lab,start.zero,lwd.nu,cex.nu,xlab.line,ylab.line,line.col){
                          if(missing(input)){                                   stop("no input defined")}
                          if(missing(bg.col)){                                  bg.col <- "white"}
                          if(missing(col.scale)){                               col.scale <- c("#8E0152","#C51B7D","#DE77AE","#F1B6DA","#FDE0EF","#F7F7F7","#E6F5D0","#B8E186","#7FBC41","#4D9221","#276419")}
                          if(missing(x.lab)){                                   x.lab <- ""}
                          if(missing(start.zero)){                              start.zero <- TRUE}
                          if(missing(lwd.nu)){                                  lwd.nu <- 2}
                          if(missing(cex.nu)){                                  cex.nu <- 1}
                          if(missing(xlab.line)){                               xlab.line <- 2.5}
                          if(missing(ylab.line)){                               ylab.line <- 2.5}
                          if(missing(line.col)){                                line.col <- "white"}

                          ###Remove NA

                          ###Make the color scale
                          colors.use <- color.scale(col.scale,input)

                          ###first plotting layer
                          if(start.zero==TRUE){x.range <- c(0,max(input,na.rm=T))}
                          if(start.zero==FALSE){x.range <- c(min(input,na.rm=T),max(input,na.rm=T))}

                          x.val <- axis.output(x.range[1],x.range[2])[[1]]
                              x.val <- x.val[c(1:which(x.val > max(x.range))[1])]
                          y.val <- hist(input,breaks=x.val,plot=FALSE)$counts

                          y.range <- c(0,max(y.val))
                          x.range <- c(min(x.val),max(x.val))

                          plot(x=x.val,y=c(y.val,y.val[length(y.val)]),type="s",xaxs="i",axes=F,xlab="",ylab="",xlim=x.range,ylim=y.range)
                          col.grad <- seq(x.range[1],x.range[2],length.out=length(colors.use)+1)
                          for(i in 2:length(col.grad)){
                              rect(xleft=col.grad[i-1],xright=col.grad[i]+10000,ybottom=-40^6,ytop=40^6,border=NA,col=colors.use[i-1])
                          }
                          points(x=x.val,y=c(y.val,y.val[length(y.val)]),type="s",col=line.col,lwd=lwd.nu,lend=1)

                          ###close it up, maxe axes etc
                          box()
                          fixed.axis <- axis.output(x.range[1],x.range[2])
                          axis(1,fixed.axis[[1]],labels=F,tcl=-0.2)
                          if(start.zero){axis(1,c(0,fixed.axis[[2]]),c(0,fixed.axis[[2]]),mgp=c(3,0.7,0),cex.axis=cex.nu,tcl=-0.5,las=1)}
                          if(start.zero==FALSE){axis(1,fixed.axis[[2]],fixed.axis[[2]],mgp=c(3,0.7,0),cex.axis=cex.nu,tcl=-0.5,las=1)}
                          mtext(side=1,text=x.lab,cex=cex.nu+0.5,line=xlab.line)
                          fixed.axis <- axis.output(y.range[1],y.range[2])
                          axis(2,fixed.axis[[1]],labels=F,tcl=-0.2)
                          axis(2,fixed.axis[[2]],fixed.axis[[2]],mgp=c(3,0.7,0),cex.axis=cex.nu,tcl=-0.5,las=1)
                          mtext(side=2,"Counts",cex=cex.nu+0.5,line=ylab.line)
                         }


axis.output <- function(min.nu,max.nu){

                        output <- NULL; output <- as.list(output)
                        for(i in sort(c(10^c(seq(-2,8,by=1)),2*10^c(seq(-2,8,by=1)),4*10^c(seq(-2,8,by=1)),6*10^c(seq(-2,8,by=1)),8*10^c(seq(-2,8,by=1))),decreasing=T)){
                            ###determine the scale of the differences
                            if(max.nu-min.nu < i*2){
                                small.axis <- seq(floor(min.nu/i)*i,ceiling(max.nu/i)*i*2,by=i/20)
                                large.axis <- seq(floor(min.nu/i)*i,ceiling(max.nu/i)*i*2,by=i/4)[-1]
                            }
                        }
                        output[[1]] <- small.axis
                        output[[2]] <- large.axis
                        return(output)
}

png('subfigures/heatmap_structures.png', 2500, 2800, pointsize = 100)
#svg('subfigures/heatmap_structures.svg')
colscale <- color.scale(colors.use = brewer.pal(11, "Reds"), y.val = as.numeric(structure_plot.df))
heatmap(structure_plot.df, scale = 'none', col=colscale, margins=c(10,8))
dev.off()
png("subfigures/heatmap_structures_scale.png")
scale.legend(input=as.numeric(structure_plot.df),start.zero=F,col.scale=colscale,line.col="black")
dev.off()


png("subfigures/heatmap_sequences.png", 2500, 2800, pointsize = 100)
sequence_plot.df <- comparisons |>
  select('Query','Subject','global_identity') |>
  pivot_wider(names_from = 'Subject', values_from = 'global_identity') |>
  dplyr::select(-'Query') |>
  as.matrix() |>
  apply(MARGIN=2, as.numeric)
colscale <- color.scale(colors.use = brewer.pal(11, "Reds"), y.val = as.numeric(sequence_plot.df))
heatmap(sequence_plot.df, scale = 'none', col=colscale, margins=c(10,8))
dev.off()
png("subfigures/heatmap_sequences_scale.png")
#par(fig=c(0,0.3,0.5,1), new =T)
scale.legend(input=as.numeric(sequence_plot.df),start.zero=F,col.scale=colscale,line.col="black")
dev.off()

  ## reorder_matrix(cutoff=25, decreasing = F) |>
  ## #apply(MARGIN=2, function(x) as.numeric(x > 25)) |>
  ##   gplots::heatmap.2(trace = "none",
  ##                   labRow = F,
  ##                   labCol = F,
  ##                   dendrogram = 'column',
  ##                   col = colorRampPalette(c("#faf3f2", "black")),
  ##                   #Rowv=FALSE,
  ##                   #Colv=FALSE,
  ##                   key=TRUE,
  ##                   xlab = 'G. pallida Sequences',
  ##                   ylab = 'M. chitwoodi Sequences')
#dev.off()

a <- comparisons |>
  select('Query','Subject','global_identity') |>
  pivot_wider(names_from = 'Subject', values_from = 'global_identity') |>
  dplyr::select(-'Query') |>
  as.matrix() |>
  apply(MARGIN=2, as.numeric) |>
  reorder_matrix(cutoff=25, decreasing = F) |>
  apply(MARGIN=2, function(x) as.numeric(x > 25))



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
  comparisons$P.adjustedLog
)

plot.df%<>% as.data.frame()
colnames(plot.df) <- c("Identity","Pval")

pvalue_log <- -1 * log10(0.05)

png("subfigures/vulcano.png",600,600)
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
  ylab("-log10(pvalue)") +
  geom_hline(yintercept=pvalue_log, linetype='dashed', col = 'red')+
  scale_x_continuous(limits = c(0,100))+
  ggtitle("") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text  = element_text(face = "bold", size = 18),
    axis.title = element_text(size = 24),
    )
dev.off()

sequence_map <- png::readPNG('subfigures/heatmap_sequences_e.png')
sequence_map <- rasterGrob(sequence_map, interpolate = T)

structure_map <- png::readPNG('subfigures/heatmap_structures_e.png')
structure_map <- rasterGrob(structure_map, interpolate = T)

vulcano_fig <- png::readPNG('subfigures/vulcano.png')|>
  rasterGrob(interpolate = T)

png(filename = "~/org-roam/alphafold_data/Panel_comparisons.png",width=2800, height=2800)
gridExtra::grid.arrange(
             sequence_map,
             structure_map,
             vulcano_fig,
             gpal_hist,
             mchit_hist,
             ncol=3,nrow=2
             #layout_matrix = cbind(c(1,1,4), c(2,2,5), c(3,3,NA))
           )
dev.off()


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


#
