test_that("multiplication works", {
  mat <- matrix(rnorm(1000,-10, 10), ncol = 10)
  calc_adj_frag(mat, 20, 10, 0.1) |> expect_no_error()
  calc_adj_frag(mat, 20, 10) |> expect_no_error()
  ridge(mat[1:4, ], mat[1:4 + 1, ], 0.1, FALSE, NULL) |>
    fragilityRow() |>
    expect_no_error()
})
