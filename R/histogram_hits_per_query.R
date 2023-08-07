
Histogram_hit_per_query <- function(gpal_hist_png = "Figures/subfigures/gpal_hist.png", mchit_hist_png = "Figures/subfigures/mchit_hist.png") {

  #png(gpal_hist_png)
 gpal_hist <- comparisons.filtered |>
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
   theme_grattan() +
   theme(axis.title.y = element_text(angle=90)) +
   labs( x = "Number of hits per G. pallida query", y = "Count") #+
    #dev.off()

#png(mchit_hist_png)
mchit_hist <-comparisons.filtered |>
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
  theme_grattan() +
  theme(axis.title.y = element_text(angle=90)) +
   labs( x = "Number of hits per M. chitwoodi query", y = "Count") #+
#dev.off()
  return(list("mchit_hist" = mchit_hist, "gpal_hist" =  gpal_hist))
}



## reorder_matrix <- function(matrix, cutoff = 0.05, decreasing = T) {
##   matrix <- matrix[,(matrix < cutoff )|>
##                     colSums() |>
##                     order(decreasing = decreasing)]
##   matrix <- matrix[(matrix < cutoff )|>
##                     rowSums() |>
##                     order(decreasing = decreasing),]
##   return(matrix)
## }
