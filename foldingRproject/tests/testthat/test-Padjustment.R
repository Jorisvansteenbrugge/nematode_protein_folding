

generate_test_frame <- function(name, n) {
  data <- data.frame(
    "Query" = rep(name, n),
    "P.value" = sample( seq(0,1, by=0.01), n)
  )
  return(data)
}

data_a <- generate_test_frame("A", 31)
data_b <- generate_test_frame("B", 41)
data_merged <- rbind(data_a, data_b)

  test_that("P adjustment single FDR", {
    adjusted <- p.adjust(
      data_a$P.value,
      method = "fdr"
    )
    adjusted_custom <- p.adjust.custom(
      data_a,
      df_N_col = "Query",
      df_P_col = "P.value",
      method = "fdr"
    )
    expect_equal(adjusted, adjusted_custom)
  })

test_that("P adjustment single bonferroni", {
    adjusted <- p.adjust(
      data_a$P.value,
      method = "bonferroni"
    )
    adjusted_custom <- p.adjust.custom(
      data_a,
      df_N_col = "Query",
      df_P_col = "P.value",
      method = "bonferroni"
    )
    expect_equal(adjusted, adjusted_custom)
  })


test_that("P adjustment double FDR", {
  adjusted_a <- p.adjust(
    data_a$P.value,
    method = "fdr"
  )
  adjusted_b <- p.adjust(
    data_b$P.value,
    method = "fdr"
  )
  adjusted <- c(adjusted_a, adjusted_b)
  adjusted_custom <- p.adjust.custom(
    data_merged,
     df_N_col = "Query",
      df_P_col = "P.value",
      method = "fdr"
  )
  expect_equal(adjusted, adjusted_custom)
})

test_that("P adjustment double bonferroni", {
  adjusted_a <- p.adjust(
    data_a$P.value,
    method = "bonferroni"
  )
  adjusted_b <- p.adjust(
    data_b$P.value,
    method = "bonferroni"
  )
  adjusted <- c(adjusted_a, adjusted_b)
  adjusted_custom <- p.adjust.custom(
    data_merged,
     df_N_col = "Query",
      df_P_col = "P.value",
      method = "bonferroni"
  )
  expect_equal(adjusted, adjusted_custom)
})
