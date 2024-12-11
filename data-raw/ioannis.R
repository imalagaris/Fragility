{
  ncols <- 8
  readings <- matrix(rnorm(1000), ncol = ncols)
  colnames(readings) <- seq_len(ncols) |> (\(x) paste0(LETTERS[x], x))()
  t_window <- 4
  t_step <- 1
  lambda <- 0.1


  (n_tps <- nrow(readings))
  (n_elec <- ncol(readings))

  electrode_list <- colnames(readings)

  # Number of steps
  (n_steps <- floor((n_tps - t_window) / t_step) + 1)

  (scaling <- 10^floor(log10(max(readings))))
  readings <- readings / scaling
}



res <- lapply(seq_len(n_steps), function(iw) {
  ## Sample indices for the selected window
  si <- seq_len(t_window - 1) + (iw - 1) * t_step
  ## measurements at time point t
  xt <- readings[si, ]
  ## measurements at time point t plus 1
  xtp1 <- readings[si + 1, ]

  ## Coefficient matrix A (adjacency matrix)
  ## each column is coefficients from a linear regression
  ## formula: xtp1 = xt*A + E
  if (is.null(lambda)) {
    Ai <- ridgesearchlambdadichomotomy(xt, xtp1, intercept = FALSE, iw = iw)
  } else {
    Ai <- ridge(xt, xtp1, intercept = FALSE, lambda = lambda, iw = iw)
  }

  R2 <- ridgeR2(xt, xtp1, Ai)

  list(si = si, xt = xt, xtp1 = xtp1, Ai = Ai, R2 = R2)
})

res[[1]]
(predOfXtp1 <- with(res[[1]], xt %*% Ai))
res[[1]]$xtp1
E <- res[[1]]$xtp1 - predOfXtp1

e <- unlist(lapply(res, `[[`, "Ai"))

dim(A) <- c(n_elec, n_elec, n_steps)
dimnames(A) <- list(
  Electrode1 = electrode_list,
  Electrode2 = electrode_list,
  Step = seq_len(n_steps)
)


R2 <- unlist(lapply(res, `[[`, "R2"))
dim(R2) <- c(n_elec, n_steps)
dimnames(R2) <- list(
  Electrode = electrode_list,
  Step = seq_len(n_steps)
)

lambdas <- rep(lambda, length(res))

x1 <- res[[1]]$xtp1
pred <- with(res[[1]], xt %*% Ai)
sst <- x1 |> apply(2, \(x) sum((x - mean(x))^2))
sse <- (x1 - pred) |> apply(2, \(x) sum(x^2))
rsq <- 1 - sse / sst



ridgeRSQ <- function(xt, xtp1, A) {
  sse <- (xtp1 - xt %*% A) |> apply(2, \(x) sum(x^2))
  sst <- apply(xt, 2, \(x) sum((x - mean(x))^2))
  1 - sse / sst
}



