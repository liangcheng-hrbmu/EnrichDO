library(EnrichDO)

test_that("doEnrich", {
  demo.data <- c(1636, 351, 102, 2932, 3077, 348, 4137, 54209, 5663, 5328, 23621, 3416, 3553)
  res <- sapply(1:5, function(m) {
    aa <- doEnrich(interestGenes = demo.data, m = m)
    result <- aa@enrich
    bb <- dplyr::filter(result, level < m)$p
    return(all(as.numeric(bb) == 1))
  })
  expect_true(all(res) == TRUE)
})

test_that("viewDetailResult and viewSummaryResult work correctly", {
  demo.data <- c(1636, 351, 102, 2932, 3077)
  test_res <- doEnrich(interestGenes = demo.data, maxGsize = 100, minGsize = 10)

  detail_df <- viewDetailResult(test_res)
  expect_true(is.data.frame(detail_df))
  expect_true(nrow(detail_df) > 0)

  summary_df <- viewSummaryResult(test_res)
  expect_true(is.data.frame(summary_df))
  expect_true("geneRatio" %in% colnames(summary_df))
  expect_true("bgRatio" %in% colnames(summary_df))
})
