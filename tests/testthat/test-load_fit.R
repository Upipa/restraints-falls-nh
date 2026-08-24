test_that("load_fit returns the shipped posterior draws", {
  draws <- load_fit()

  expect_s3_class(draws, "draws_df")
  expect_equal(nrow(draws), 4000L)
  expect_true(all(c("delta_fall", "rr_fall", "delta_major") %in% names(draws)))
  expect_false(any(grepl("^x_rep\\[", names(draws))))
})
