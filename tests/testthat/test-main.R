test_that("multiplication works", {
  readings <- matrix(rnorm(1000,-10, 10), ncol = 10)
  calc_adj_frag(readings, 20, 10, 0.1) |> expect_no_error()
  calc_adj_frag(readings, 20, 10) |> expect_no_error()
})
