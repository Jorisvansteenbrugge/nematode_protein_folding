

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
