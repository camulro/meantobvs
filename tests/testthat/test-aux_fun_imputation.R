#Tests for the internal imputation machinery: mice.imputation,
#buildimputation and extimputation (user supplied
#imputations) and the rank deficient design matrix model.rankdef.matrix.

test_that("model.rankdef.matrix returns an intercept for a response only model", {
  mf <- model.frame(y ~ 1, lm.test.data(), na.action = NULL)
  out <- meantobvs:::model.rankdef.matrix(mf)

  expect_equal(dim(out), c(24, 1))
  expect_identical(colnames(out), "(Intercept)")
  expect_identical(as.vector(out), rep(1, 24))
})

test_that("model.rankdef.matrix builds the full design matrix for numeric covariates", {
  mf <- model.frame(y ~ x1 + x2, lm.test.data(), na.action = NULL)
  out <- meantobvs:::model.rankdef.matrix(mf)

  expected <- model.matrix(y ~ x1 + x2, mf)
  expect_equal(dim(out), c(24, 3))
  expect_identical(colnames(out), c("(Intercept)", "x1", "x2"))
  expect_equal(as.matrix(out), as.matrix(expected))
})

test_that("model.rankdef.matrix keeps all dummies of a factor", {
  data <- lm.test.data()
  data$x1 <- factor(rep(c("L1", "L2", "L3"), 8))
  mf <- model.frame(y ~ x1 + x2, data, na.action = NULL)

  out <- meantobvs:::model.rankdef.matrix(mf)

  expect_equal(dim(out), c(24, 5))
  expect_identical(colnames(out), c("(Intercept)", "x1L1", "x1L2", "x1L3", "x2"))
  expect_identical(as.vector(out[, "(Intercept)"]), rep(1, 24))
  expect_identical(as.vector(out[, "x1L1"]), rep(c(1, 0, 0), 8))
})

test_that("mice.imputation returns an imputation array", {
  data <- lm.test.data(missing = TRUE)
  fulldataframe <- model.frame(y ~ x1 + x2, data, na.action = NULL)

  imp <- mice.imputation(fulldataframe, n.imp = 2, seed = 1, parallel = FALSE)

  expect_named(imp, c("imputation.array", "logEvents"))
  expect_s3_class(imp$imputation.array, "meantobvs.imputation")
  expect_equal(dim(imp$imputation.array), c(24, 3, 2))
  expect_identical(dimnames(imp$imputation.array)[[2]],
                   c("(Intercept)", "x1", "x2"))
  expect_identical(dimnames(imp$imputation.array)[[3]], c("1", "2"))

  #the missing values on x1 (rows 4 and 11) are imputed in every slice
  for (s in 1:2) {
    expect_false(is.na(imp$imputation.array[4, "x1", s]))
    expect_false(is.na(imp$imputation.array[11, "x1", s]))
  }
})

test_that("mice.imputation is reproducible for a given seed", {
  data <- lm.test.data(missing = TRUE)
  fulldataframe <- model.frame(y ~ x1 + x2, data, na.action = NULL)

  imp1 <- mice.imputation(fulldataframe, n.imp = 2, seed = 1, parallel = FALSE)
  imp2 <- mice.imputation(fulldataframe, n.imp = 2, seed = 1, parallel = FALSE)
  expect_identical(imp1, imp2)
})

test_that("mice.imputation with the mean method imputes the observed means", {
  data <- lm.test.data(missing = TRUE)
  fulldataframe <- model.frame(y ~ x1 + x2, data, na.action = NULL)

  imp <- mice.imputation(fulldataframe, n.imp = 2, imp.mice.method = "mean",
                         seed = 1, parallel = FALSE)

  truth <- mean(data$x1, na.rm = TRUE)
  for (s in 1:2) {
    expect_equal(unname(imp$imputation.array[4, "x1", s]), truth, tolerance = 1e-6)
    expect_equal(unname(imp$imputation.array[11, "x1", s]), truth, tolerance = 1e-6)
    #the fully observed covariate x2 is left unchanged
    expect_equal(unname(imp$imputation.array[, "x2", s]), data$x2, tolerance = 1e-12)
  }
})

test_that("mice.imputation imputes any number of datasets", {
  data <- lm.test.data(missing = TRUE)
  fulldataframe <- model.frame(y ~ x1 + x2, data, na.action = NULL)

  imp <- mice.imputation(fulldataframe, n.imp = 3, seed = 1, parallel = FALSE)
  expect_equal(dim(imp$imputation.array)[3], 3)
})

test_that("mice.imputation returns the elapsed time with time.test = TRUE", {
  data <- lm.test.data(missing = TRUE)
  fulldataframe <- model.frame(y ~ x1 + x2, data, na.action = NULL)

  expect_s3_class(
    mice.imputation(fulldataframe, n.imp = 1, seed = 1, parallel = FALSE,
                    time.test = TRUE),
    "difftime"
  )
})

test_that("mice.imputation accepts a visit sequence", {
  data <- lm.test.data(missing = TRUE)
  fulldataframe <- model.frame(y ~ x1 + x2, data, na.action = NULL)

  imp <- mice.imputation(fulldataframe, n.imp = 1, visit.seq = c("x1", "x2"),
                         seed = 1, parallel = FALSE)
  expect_equal(dim(imp$imputation.array)[3], 1)
})

test_that("buildimputation imputes and records the imputation info", {
  data <- lm.test.data(missing = TRUE)

  output <- capture.output({
    imputation <- meantobvs:::buildimputation(
      c(x1 = TRUE, x2 = FALSE), y ~ x1 + x2, data, imp.predict.mat = NULL,
      n.imp = 2, maxit = 5, n = 24, q = 2, p0 = 1, imp.mice.method = "pmm",
      imp.seed = 1, parallelmice = FALSE, n.core = NULL,
      obsnotNA = 1:24, ordvars = c("(Intercept)", "x1", "x2")
    )
  })

  expect_gt(length(output), 0)

  imputation.array <- imputation$imputation.array
  expect_equal(dim(imputation.array), c(24, 3, 2))
  expect_identical(dimnames(imputation.array)[[2]],
                   c("(Intercept)", "x1", "x2"))
  expect_identical(unique(c(imputation.array[, "(Intercept)", ])), 1)
  expect_false(is.na(imputation.array[4, "x1", ]))
  expect_false(is.na(imputation.array[11, "x1", ]))

  info <- imputation$imp.info
  expect_identical(info$n.imp, 2)
  expect_identical(info$imp.seed, 1)
  expect_identical(info$imp.mice.method, "pmm")
  expect_false(info$parallelmice)
  expect_null(info$imp.predict.mat)
  expect_identical(unname(info$NAvars), 1L)
  expect_true(is.data.frame(info$loggedEvents) || is.null(info$loggedEvents))
})

test_that("buildimputation is reproducible for a given seed", {
  data <- lm.test.data(missing = TRUE)

  run <- function() {
    imputation <- NULL
    capture.output({
      imputation <- meantobvs:::buildimputation(
        c(x1 = TRUE, x2 = FALSE), y ~ x1 + x2, data, imp.predict.mat = NULL,
        n.imp = 2, maxit = 5, n = 24, q = 2, p0 = 1, imp.mice.method = "pmm",
        imp.seed = 1, parallelmice = FALSE, n.core = NULL,
        obsnotNA = 1:24, ordvars = c("(Intercept)", "x1", "x2")
      )
    })
    imputation$imputation.array
  }

  expect_identical(run(), run())
})

test_that("buildimputation with the mean method imputes observed means", {
  data <- lm.test.data(missing = TRUE)

  capture.output({
    imputation <- meantobvs:::buildimputation(
      c(x1 = TRUE, x2 = FALSE), y ~ x1 + x2, data, imp.predict.mat = NULL,
      n.imp = 2, maxit = 5, n = 24, q = 2, p0 = 1, imp.mice.method = "mean",
      imp.seed = 1, parallelmice = FALSE, n.core = NULL,
      obsnotNA = 1:24, ordvars = c("(Intercept)", "x1", "x2")
    )
  })

  truth <- mean(data$x1, na.rm = TRUE)
  expect_equal(unname(imputation$imputation.array[4, "x1", ]),
               rep(truth, 2), tolerance = 1e-6)
  expect_equal(unname(imputation$imputation.array[11, "x1", ]),
               rep(truth, 2), tolerance = 1e-6)
})

test_that("buildimputation defaults to the sequential mice when sizes are small", {
  data <- lm.test.data(missing = TRUE)

  capture.output({
    imputation <- meantobvs:::buildimputation(
      c(x1 = TRUE, x2 = FALSE), y ~ x1 + x2, data, imp.predict.mat = NULL,
      n.imp = 2, maxit = 5, n = 24, q = 2, p0 = 1, imp.mice.method = "pmm",
      imp.seed = 1, parallelmice = NULL, n.core = NULL,
      obsnotNA = 1:24, ordvars = c("(Intercept)", "x1", "x2")
    )
  })
  expect_false(imputation$imp.info$parallelmice)
})

test_that("buildimputation accepts a prediction matrix on complete data", {
  data <- lm.test.data() #no missings, so every which(NAvars) check is skipped
  pmat <- matrix(c(0, 1, 1, 0), 2, 2,
                 dimnames = list(c("x1", "x2"), c("x1", "x2")))

  capture.output({
    imputation <- meantobvs:::buildimputation(
      c(x1 = FALSE, x2 = FALSE), y ~ x1 + x2, data, imp.predict.mat = pmat,
      n.imp = 2, maxit = 5, n = 24, q = 2, p0 = 1, imp.mice.method = "pmm",
      imp.seed = 1, parallelmice = FALSE, n.core = NULL,
      obsnotNA = 1:24, ordvars = c("(Intercept)", "x1", "x2")
    )
  })

  expect_equal(dim(imputation$imputation.array), c(24, 3, 2))
  expect_identical(imputation$imp.info$imp.predict.mat, pmat)
})

test_that("buildimputation rejects unknown row variables in imp.predict.mat", {
  data <- lm.test.data()
  pmat <- matrix(0, 1, 2, dimnames = list("zzz", c("x1", "x2")))

  expect_test_error(
    meantobvs:::buildimputation(
      c(x1 = FALSE, x2 = FALSE), y ~ x1 + x2, data, imp.predict.mat = pmat,
      n.imp = 2, maxit = 5, n = 24, q = 2, p0 = 1, imp.mice.method = "pmm",
      imp.seed = 1, parallelmice = FALSE, n.core = NULL,
      obsnotNA = 1:24, ordvars = c("(Intercept)", "x1", "x2")
    ),
    "Row variables in imp.predict.mat not found in data"
  )
})

test_that("buildimputation rejects unknown column variables in imp.predict.mat", {
  data <- lm.test.data()
  pmat <- matrix(0, 2, 1, dimnames = list(c("x1", "x2"), "zzz"))

  expect_test_error(
    meantobvs:::buildimputation(
      c(x1 = FALSE, x2 = FALSE), y ~ x1 + x2, data, imp.predict.mat = pmat,
      n.imp = 2, maxit = 5, n = 24, q = 2, p0 = 1, imp.mice.method = "pmm",
      imp.seed = 1, parallelmice = FALSE, n.core = NULL,
      obsnotNA = 1:24, ordvars = c("(Intercept)", "x1", "x2")
    ),
    "Column variables in imp.predict.mat not found in data"
  )
})

test_that("extimputation builds the array from supplied imputations", {
  data <- lm.test.data(missing = TRUE)
  imputations <- lm_test_imputations(data)

  out <- meantobvs:::extimputation(
    y ~ x1 + x2, imputations, n0 = nrow(data),
    framefull = model.frame(y ~ x1 + x2, data, na.action = NULL),
    ordvars = c("(Intercept)", "x1", "x2"), obsnotNA = 1:nrow(data),
    p0 = 1, NAvars = c(x1 = TRUE, x2 = FALSE)
  )

  expect_equal(dim(out$imputation.array), c(24, 3, 1))
  expect_identical(dimnames(out$imputation.array)[[2]],
                   c("(Intercept)", "x1", "x2"))
  expect_equal(out$imp.info$n.imp, 1)
  expect_equal(out$n.imp, 1)
  expect_identical(unname(out$imp.info$NAvars), 1L)

  m <- out$imputation.array
  expect_identical(unique(c(m[, "(Intercept)", ])), 1)
  #the missing cells were imputed with the observed mean by the fixture
  truth <- mean(data$x1, na.rm = TRUE)
  expect_equal(unname(m[4, "x1", 1]), truth, tolerance = 1e-6)
  expect_equal(unname(m[11, "x1", 1]), truth, tolerance = 1e-6)
  expect_equal(unname(m[, "x2", 1]), data$x2, tolerance = 1e-12)
})

test_that("extimputation accepts a list of imputed datasets", {
  data <- lm.test.data(missing = TRUE)
  imputations <- lm_test_imputations(data)
  imps_list <- replicate(2, data.frame(imputations[, , 1]), simplify = FALSE)

  out <- meantobvs:::extimputation(
    y ~ x1 + x2, imps_list, n0 = nrow(data),
    framefull = model.frame(y ~ x1 + x2, data, na.action = NULL),
    ordvars = c("(Intercept)", "x1", "x2"), obsnotNA = 1:nrow(data),
    p0 = 1, NAvars = c(x1 = TRUE, x2 = FALSE)
  )

  expect_equal(dim(out$imputation.array), c(24, 3, 2))
  expect_equal(out$n.imp, 2)
  expect_identical(dimnames(out$imputation.array)[[2]],
                   c("(Intercept)", "x1", "x2"))
})

test_that("extimputation drops observations with a missing response", {
  data <- gd25_response_na_data() #y[3] is NA
  imputations <- lm_test_imputations(data)

  out <- meantobvs:::extimputation(
    y ~ x1 + x2, imputations, n0 = nrow(data),
    framefull = model.frame(y ~ x1 + x2, data, na.action = NULL),
    ordvars = c("(Intercept)", "x1", "x2"),
    obsnotNA = setdiff(1:nrow(data), 3),
    p0 = 1, NAvars = c(x1 = TRUE, x2 = FALSE)
  )

  expect_equal(dim(out$imputation.array), c(23, 3, 1))
})

test_that("extimputation requires an array or a list", {
  data <- lm.test.data()

  expect_test_error(
    meantobvs:::extimputation(
      y ~ x1 + x2, c(x1 = 1, x2 = 2), n0 = nrow(data),
      framefull = model.frame(y ~ x1 + x2, data, na.action = NULL),
      ordvars = c("(Intercept)", "x1", "x2"), obsnotNA = 1:nrow(data),
      p0 = 1, NAvars = c(x1 = TRUE, x2 = FALSE)
    ),
    "Imputations should be given as an array or list"
  )
})

test_that("extimputation validates array imputations", {
  data <- lm.test.data()
  framefull <- model.frame(y ~ x1 + x2, data, na.action = NULL)

  wrong.nrow <- array(0, dim = c(23, 3, 1),
                      dimnames = list(NULL, c("y", "x1", "x2"), "1"))
  expect_test_error(
    meantobvs:::extimputation(
      y ~ x1 + x2, wrong.nrow, n0 = nrow(data), framefull = framefull,
      ordvars = c("(Intercept)", "x1", "x2"), obsnotNA = 1:nrow(data),
      p0 = 1, NAvars = c(x1 = TRUE, x2 = FALSE)
    ),
    "24xn.varsxn.imp"
  )

  no.names <- array(0, dim = c(24, 3, 1))
  expect_test_error(
    meantobvs:::extimputation(
      y ~ x1 + x2, no.names, n0 = nrow(data), framefull = framefull,
      ordvars = c("(Intercept)", "x1", "x2"), obsnotNA = 1:nrow(data),
      p0 = 1, NAvars = c(x1 = TRUE, x2 = FALSE)
    ),
    "imp.datasets should have dimension names"
  )

  missing.var <- array(0, dim = c(24, 2, 1),
                       dimnames = list(NULL, c("y", "x1"), "1"))
  expect_test_error(
    meantobvs:::extimputation(
      y ~ x1 + x2, missing.var, n0 = nrow(data), framefull = framefull,
      ordvars = c("(Intercept)", "x1", "x2"), obsnotNA = 1:nrow(data),
      p0 = 1, NAvars = c(x1 = TRUE, x2 = FALSE)
    ),
    "Variables: x2; not given by imp.datasets"
  )
})

test_that("extimputation validates list imputations", {
  data <- lm.test.data()
  framefull <- model.frame(y ~ x1 + x2, data, na.action = NULL)

  wrong.nrow <- list(array(0, dim = c(23, 3, 1),
                           dimnames = list(NULL, c("y", "x1", "x2"), "1")))
  expect_test_error(
    meantobvs:::extimputation(
      y ~ x1 + x2, wrong.nrow, n0 = nrow(data), framefull = framefull,
      ordvars = c("(Intercept)", "x1", "x2"), obsnotNA = 1:nrow(data),
      p0 = 1, NAvars = c(x1 = TRUE, x2 = FALSE)
    ),
    "24xn.vars\\) matrices"
  )

  no.names <- list(array(0, dim = c(24, 3, 1)))
  expect_test_error(
    meantobvs:::extimputation(
      y ~ x1 + x2, no.names, n0 = nrow(data), framefull = framefull,
      ordvars = c("(Intercept)", "x1", "x2"), obsnotNA = 1:nrow(data),
      p0 = 1, NAvars = c(x1 = TRUE, x2 = FALSE)
    ),
    "imp.datasets should have dimension names"
  )

  missing.var <- list(array(0, dim = c(24, 2, 1),
                            dimnames = list(NULL, c("y", "x1"), "1")))
  expect_test_error(
    meantobvs:::extimputation(
      y ~ x1 + x2, missing.var, n0 = nrow(data), framefull = framefull,
      ordvars = c("(Intercept)", "x1", "x2"), obsnotNA = 1:nrow(data),
      p0 = 1, NAvars = c(x1 = TRUE, x2 = FALSE)
    ),
    "Variables: x2; not given by imp.datasets"
  )
})
