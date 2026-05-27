test_that("db_connect has expected arguments", {
  args <- formals(db_connect)
  expect_named(args, c("uid", "pwd"))
})

test_that("db_connect fails gracefully without credentials", {
  expect_error(db_connect(uid = "", pwd = ""))
})
