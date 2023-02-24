

#' Adjusted pvalue based on the number of tests per query sequence
#'
#' @param comparison_df dataframe containing the pvalues and querys to correct multiple testing for
#' @param df_N_col column name containg the non-unique query names that will be used to determine 'n' comparisons that were done
#' @param df_P_col column name containing the pvalues
#' @param method method to for multiple testing correction, can be any implemented in `p.adjust`
#' @author Joris van Steenbrugge
#' @export
#' @examples
#' data <- data.frame( "Query" = rep(name, n), "P.value" = sample( seq(0,1, by=0.01), n))
#' data$p.adjusted <- p.adjust.custom(data, df_N_col = "Query", df_P_col = "P.value")
#'
p.adjust.custom <- function(comparison_df, df_N_col = "Query" , df_P_col = "P-value", method = "fdr" ) {
  unique_querries <- comparison_df[[df_N_col]] |> unique()
  adjusted_pvals <- c()
  for ( query in unique_querries ){
    comparisons_c <- comparison_df |>
      filter( !!sym(df_N_col) == query )
    comparisons_c$P.adjusted <- p.adjust(
      comparisons_c[[df_P_col]],
      method = method
    )
    adjusted_pvals <- c(adjusted_pvals, comparisons_c$P.adjusted)
  }
  return(adjusted_pvals)
}


#' Plot a density of the pairwise sequences identities of all proteins vs only the structurally significant proteins
#'
#' @param df_all dataframe containing all comparisons with a column called `global_identity`
#' @param df_sig dataframe containing only the significant comparisons with a column called `global_identity`
#' @author Joris van Steenbrugge
#' @export
#' @examples
#' plot_identities(comparisons, comparisons.filtered)
plot_identities <- function(df_all, df_sig) {
  all_identity <- cbind(
    as.numeric(df_all$global_identity),
    rep("All", times = nrow(df_all)))

  sig_identity <- cbind(
    df_sig$global_identity,
    rep("Significant only", times = nrow(df_sig))
  )

  identities <- rbind(all_identity, sig_identity) %>% data.frame
  colnames(identities) <- c("Identity", "Set")
  identities$Identity %<>% as.numeric

  plt <- identities %>%
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

  return(plt)
}

#' Plot a vulcano style plot of the structural alignment pvalue versus the protein sequence identity for the
#' pairwise comparisons
#' @param df dataframe containing all comparisons with a column called `global_identity` and `P.adjusted`
#' @author Joris van Steenbrugge
#' @export
#' @examples
#' plot_vulcano(comparisons)
plot_vulcano <- function(df) {
  plot.df <- cbind(
    as.numeric(df$global_identity),
    -1 * log10(as.numeric(df$P.adjusted + 0.00000000001))
  )

  plot.df%<>% as.data.frame()
  colnames(plot.df) <- c("Identity","Pval")

  pvalue_log <- -1 * log10(0.05)
  plt <- plot.df %>%
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
    geom_hline(yintercept = pvalue_log) +
    geom_vline(xintercept = 30) +
    theme_bw()

  return(plt)
}
