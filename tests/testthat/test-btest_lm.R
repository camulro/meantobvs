#Tests for the btest.lm function

#Check that the result is of the desired form and basic outputs without factors
test_that("btest.lm computes posterior probabilities for named models", {
  #build test deterministic data
  data <- lm.test.data()

  models <- list(null = y ~ 1, x1 = y ~ x1, full = y ~ x1 + x2)

  #capture all text written by the function on console
  output <- capture.output({
    result <- btest.lm(
      data = data,
      models = models,
      BF.method = "BIC",
      prior.betas = "gZellner" #does not apply with BIC
    )
  })

  expect_gt(length(output), 0) #output is not empty (greater than 0)

  expect_s3_class(result, "meantobtest") #checks class(result)

  expect_identical(result$nullmodel, "null")
  expect_identical(names(result$lBFi0), c("null.to.null", "x1.to.null", "full.to.null"))
  expect_identical(result$lBFi0[["null.to.null"]], 0)
  expect_identical(result$models, models)

  expect_equal(sum(result$PostProbi), 1) #check probs
  expect_equal(sum(result$priorprobs), 1)
  expect_equal(result$BF.method, "BIC")
  expected.PostProbi <- exp(result$lBFi0) * result$priorprobs /
    sum(exp(result$lBFi0) * result$priorprobs)
  #check computed inclprobs
  expect_equal(unname(result$PostProbi), unname(expected.PostProbi))
  expect_equal(output[1], "")

  expect_null(result$prior.betas)
  expect_null(result$imp.info)

  expect_length(result$modelspool, 3)
  expect_length(result$PostProbi, 3)
  expect_length(result$priorprobs, 3)
  expect_length(result$prior.models, 1)

  expect_true(all(vapply(result$modelspool, inherits, logical(1), "lm"))) #checks class
  expect_true(all(is.finite(result$lBFi0)))
  expect_false(anyNA(result$lBFi0)) #check there are no NAs
  expect_false(anyNA(result))
})

#Check other options
test_that("btest.lm correctly performs with other prior methods", {
  #build test deterministic data
  data <- lm.test.data()

  models <- list(null = y ~ 1, x1 = y ~ x1, full = y ~ x1 + x2)

  output <- capture.output({
    result <- btest.lm(
      data = data,
      models = models,
      BF.method = "gprior",
      prior.models = "User",
      priorprobs = c(2, 1, 1)
    )
  })

  expect_gt(length(output), 0) #output is not empty (greater than 0)

  expect_identical(result$prior.betas, "gZellner")
  expect_identical(result$prior.models, "User")

  expect_equal(sum(result$PostProbi), 1) #check probs
  expect_equal(sum(result$priorprobs), 1)
  expected.PostProbi <- exp(result$lBFi0) * result$priorprobs /
    sum(exp(result$lBFi0) * result$priorprobs)
  #check computed inclprobs
  expect_equal(unname(result$PostProbi), unname(expected.PostProbi))
  expect_equal(output[1], "")

  expect_null(result$imp.info)
})

#Check how external imputed data is treated
test_that("btest.lm uses supplied imputations", {
  #Build test data with missings and impute with the mean
  data <- lm.test.data(missing = TRUE)
  imputations <- lm.test.imputations(data)

  models <- list(null = y ~ 1, x1 = y ~ x1, full = y ~ x1 + x2)

  output <- capture.output({
    result <- btest.lm(
      data = data,
      models = models,
      BF.method = "BIC",
      imp.datasets = imputations,
      n.imp = 3, parallelmice = TRUE #to check imp.info
    )
  })

  #recover imputed data used
  imp.array <- unserialize(
    memDecompress(
      result$compress.imp.array, type = "xz")
  )

  expect_true("compress.imp.array" %in% names(result))

  expect_equal(unname(imp.array[,c("x1", "x2"),]),
               unname(imputations[,c("x1", "x2"),]))

  expect_equal(result$imp.info$n.imp, 1)
  expect_equal(sum(result$PostProbi), 1)
  expect_equal(sum(result$priorprobs), 1)

  expect_null(result$imp.info$parallelmice)
})

#Check that probabilities are computed correctly
test_that("btest.lm's posterior is correct", {
  data <- lm.test.data()
  #Full enumeration of models to compare with bvs output
  models <- list(null = y ~ 1, x1 = y ~ x1, x2 = y ~ x2, full = y ~ x1 + x2)
  #Check all BF available
  pb <- lm.prior.betas()

  for (j in pb) {

    result.bvs <- bvs.lm(
      y ~ x1 + x2,
      data = data,
      BF.method = "gprior",
      prior.betas = j,
      prior.models = "Constant"
    )

    result.Btest <- btest.lm(
      data = data,
      models = models,
      BF.method = "gprior",
      prior.betas = j
    )

    expect_identical(result.bvs$prior.betas, result.Btest$prior.betas)
    expect_identical(result.bvs$prior.models, result.Btest$prior.models)

    expect_equal(unname(result.bvs$modelsprob[c(4,2,3,1), "Post"]),
                 unname(result.Btest$PostProbi))
    expect_equal(unname(result.bvs$priorprobs[c(1,3)]),
                 unname(result.Btest$priorprobs[c(1,4)]))
  }
})

#Check that errors are given correctly and other printed messages
test_that("btest.lm reports errors and warnings", {
  data <- lm.test.data()

  expect.test.error(
    btest.lm(data = data, models = c("M0", "M1")),
    "Argument models should be a list."
  )

  expect.test.error(
    btest.lm(data = data, models = list(y ~ 1, y ~ x1), null.model = 1),
    "Please provide a name for the competing models."
  )

  expect.test.error(
    btest.lm(data = data, models = list(y ~ 1)),
    "At least two models have to be provided."
  )

  models <- list(null = y ~ 1, x1 = y ~ x1)

  expect.test.error(
    btest.lm(data = data, models = models, null.model = "unknown"),
    "null model provided is not in the list"
  )

  expect.test.error(
    btest.lm(data = data, models = models, null.model = "x1"),
    "The given null model does not coincide with the one with the"
  )

  models2 <- list(x1 = y ~ x1, x2 = y ~ x2, full = y ~ x1 + x2)

  expect.test.error(
    btest.lm(data = data, models = models2, null.model = "x2"),
    "The simplest \\(null\\) model may not be nested in all the others\\."
  )

  expect.test.error(
    btest.lm(data = data, models = models, prior.models = "SB"),
    "Only priors 'ScottBerger', 'Constant' and 'User' supported."
  )

  expect.test.error(
    btest.lm(data = data, models = models,
                    BF.method = "gprior", prior.betas = "gprior"),
    "prior.betas must be one of 'gZellner', 'Robust',"
  )

  expect.test.error(
    btest.lm(data = data, models = list(null = y ~ -1, full = y ~ x1 + x2),
             null.model = "null"),
    "Null model shall contain the intercept term"
  )

  expect.test.error(
    btest.lm(data = data, models = list(null = y ~ -1, full = y ~ x1 + x2)),
    "Null model shall contain the intercept term"
  )

  output <- capture.output({
    result <- btest.lm(
      data = data,
      models = models2
    )

    print(result)
  })

  expect_equal(output[14],
               "log Bayes factors (expressed in relation to x2)")

})
