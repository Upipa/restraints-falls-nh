test_that("prepare_model_data reproduces the analysis data structure", {
  prep <- prepare_model_data()
  sd <- prep$stan_data

  expect_identical(sd$K, 47L)
  expect_identical(sd$N_x, 2766L)
  expect_identical(sd$N_multi, 2174L)
  expect_identical(sd$N_bin, 350L)
  expect_identical(sd$N_mm, 349L)
  expect_identical(sd$N_bm, 70L)
})

test_that("multinomial cells are non-negative and sum to n_y", {
  prep <- prepare_model_data()
  sd <- prep$stan_data

  expect_gte(min(sd$c_multi), 0)
  expect_gte(min(sd$c_mm), 0)
  expect_equal(unname(rowSums(sd$c_multi)), sd$n_y_multi)
  expect_equal(unname(rowSums(sd$c_mm)), sd$n_y_mm)
})

test_that("position maps partition the observation groups", {
  sd <- prepare_model_data()$stan_data

  expect_equal(sort(c(sd$pos_multi, sd$pos_bin)), seq_len(sd$N_both))
  expect_equal(sort(c(sd$pos_mm, sd$pos_bm)), seq_len(sd$N_y_only))
})

test_that("the quadratic term is standardized on the pooled observations", {
  sd <- prepare_model_data()$stan_data
  pooled <- c(sd$t2_x, sd$t2_mm, sd$t2_bm)

  expect_equal(sd(pooled), 1, tolerance = 0.01)
  expect_lt(max(abs(pooled)), 5)
})

test_that("missing columns produce an informative error", {
  expect_snapshot(
    prepare_model_data(data.frame(anno = 2020)),
    error = TRUE
  )
})
