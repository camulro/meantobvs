#Tests for the bvs.lm function

#Check that the result is of the desired form and basic outputs without factors
test_that("bvs.lm returns the full posterior for complete data", {
  #build test deterministic data
  data <- lm.test.data()

  #capture all text written by the function on console
  output <- capture.output({
    result <- bvs.lm(
      y ~ x1 + x2,
      data = data,
      BF.method = "BIC",
      prior.betas = "gZellner", #does not apply with BIC
      n.keep = "all"
    )
  })

  expect_gt(length(output), 0) #output is not empty (greater than 0)
  expect_s3_class(result, "meantobvs") #checks class(result)
  expect_s3_class(result$lmnull, "lm") #checks class(result$lmnull)
  expect_s3_class(result$lmfull, "lm") #checks class(result$lmfull)

  expect_identical(result$method, "Full") #Check names
  expect_identical(result$variables, c("x1", "x2"))
  expect_identical(names(result$inclprob), c("x1", "x2"))
  expect_identical(names(result$postprobdim), as.character(1:3))
  expect_identical(colnames(result$modelsprob), c("x1", "x2", "Post"))
  expect_identical(rownames(result$modelsprob), as.character(1:4))

  expect_equal(length(result$prior.models), 1)
  expect_equal(result$p, 2)
  expect_equal(result$n, nrow(data))
  expect_equal(result$k, 1) #Only intercept fixed
  expect_equal(sum(result$modelsprob[, "Post"]), 1) #check probs
  expect_equal(sum(result$priorprobs), 1)
  expect_equal(length(result$priorprobs), 2 + 1)
  expect_equal(sum(result$postprobdim), 1)
  expect_equal(length(result$postprobdim), 2 + 1)
  expect_equal(length(result$HPMbin), 2) #check summaries
  expect_equal(length(result$MPMbin), 2)
  expected.inclprob <- colSums(result$modelsprob[, c("x1", "x2")] *
    result$modelsprob[, "Post"])
  #check computed inclprobs
  expect_equal(unname(result$inclprob), unname(expected.inclprob))
  expect_equal(output[5],
               "  Competing variables: x1 x2 ")
  expect_equal(output[6],
               "The problem has a total of 4 competing models.")
  expect_equal(output[7],
               "Of these, the  4 most probable (a posteriori) are kept.")

  expect_null(result$prior.betas)
  expect_null(result$imp.info)

  expect_false(anyNA(result$modelsprob)) #check there are no NAs
  expect_false(anyNA(result))
  expect_true(all(is.finite(result$modelsprob)))
  expect_true(is.finite(result$C))
})

#Check other options
test_that("bvs.lm correctly performs with other prior methods", {
  #build test deterministic data
  data <- lm.test.data()

  #capture all text written by the function on console
  output <- capture.output({
    result <- bvs.lm(
      y ~ x1 + x2,
      data = data,
      BF.method = "gprior",
      n.keep = 1,
      prior.models = "User",
      priorprobs = c(2, 1, 1)
    )
  })

  expect_gt(length(output), 0) #output is not empty (greater than 0)

  expect_identical(rownames(result$modelsprob), as.character(1))
  expect_identical(result$prior.betas, "gZellner")
  expect_identical(result$prior.models, "User")

  expect_equal(result$priorprobs, c(`1` = 1/2, `2` = 1/4, `3` = 1/4))
  expect_equal(sum(result$priorprobs), 1) #check probs
  expect_equal(length(result$priorprobs), 2 + 1)
  expect_equal(sum(result$postprobdim), 1)
  expect_equal(length(result$postprobdim), 2 + 1)
  # expected.inclprob <- result$modelsprob[, c("x1", "x2")] *
  #                      result$modelsprob[, "Post"]
  # #check computed inclprobs
  # expect_equal(unname(result$inclprob), unname(expected.inclprob))

  expect_null(result$imp.info)
})

#Check how external imputed data is treated
test_that("bvs.lm uses supplied imputations", {
  #Build test data with missings and impute with the mean
  data <- lm.test.data(missing = TRUE)
  imputations <- lm.test.imputations(data)

  output <- capture.output({
    result <- suppressWarnings(bvs.lm(
      y ~ x1 + x2,
      data = data,
      imp.datasets = imputations,
      n.imp = 3, parallelmice = TRUE #to check imp.info
    ))
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
  expect_equal(nrow(result$modelsprob), 4)
  expect_equal(sum(result$modelsprob[, "Post"]), 1)

  expect_null(result$imp.info$parallelmice)
})

#Check that probabilities computed are the same as other packages
test_that("bvs.lm's posterior is correct", {
  data <- lm.test.data()
  #Check all BF available
  pb <- lm.prior.betas()[-3] #Liangetal fails in BayesVarSel

  #just intercept fixed
  for (j in pb) {

    result.bvs <- bvs.lm(
      y ~ x1 + x2,
      data = data,
      BF.method = "gprior",
      prior.betas = j
    )

    result.other <- BayesVarSel::Bvs(
      y ~ x1 + x2,
      data = data,
      prior.betas = j,
      n.keep = 4
    )

    expect_equal(nrow(result.bvs$modelsprob), 4)
    expect_equal(sum(result.bvs$modelsprob[, "Post"]), 1)
    expect_equal(unname(result.bvs$modelsprob[, "Post"]),
                 result.other$modelsprob[, "prob"])
    expect_equal(result.bvs$inclprob,
                 result.other$inclprob)
    expect_equal(result.bvs$postprobdim,
                 result.other$postprobdim)
    expect_equal(result.bvs$HPMbin,
                 result.other$HPMbin)
  }

  #Build a noise variable
  withr::with_seed(1, #set seed for reproducibility
                   x3 <- rnorm(nrow(data))
  )
  data$x3 <- x3

  #x1 fixed
  for (j in pb) {

    result.bvs <- bvs.lm(
      y ~ x1 + x2 + x3,
      null.model = y ~ x1,
      data = data,
      BF.method = "gprior",
      prior.betas = j
    )

    result.other <- BayesVarSel::Bvs(
      y ~ x1 + x2 + x3,
      null.model = y ~ x1,
      data = data,
      prior.betas = j,
      n.keep = 4
    )

    expect_equal(nrow(result.bvs$modelsprob), 4)
    expect_equal(sum(result.bvs$modelsprob[, "Post"]), 1)
    #robust prior does not provide the exact probabilities:
    expect_equal(unname(result.bvs$modelsprob[, "Post"]),
                 result.other$modelsprob[, "prob"], tolerance = 10^6)
    expect_equal(result.bvs$inclprob,
                 result.other$inclprob)
    expect_equal(result.bvs$postprobdim,
                 result.other$postprobdim, tolerance = 10^6)
    expect_equal(result.bvs$HPMbin,
                 result.other$HPMbin)
  }
})

#Check that errors are given correctly and other printed messages
test_that("bvs.lm reports errors and warnings", {
  data <- lm.test.data()
  #Build a noise variable
  withr::with_seed(1, #set seed for reproducibility
    x3 <- rnorm(nrow(data))
  )
  data$x3 <- x3

  expect.test.error(
    bvs.lm(y ~ x1, data = data, null.model = z ~ 1),
    "response in the full and null model does not coincide"
  )

  expect.test.error(
    bvs.lm(y ~ x1 + x2, data = data, null.model = y ~ x3),
    "Error in var: x3; null model not nested in full"
  )

  expect.test.error(
    bvs.lm(y ~ x1, data = data, null.model = y ~ x1),
    "number of fixed variables is equal to the number of"
  )

  expect.test.error(
    bvs.lm(y ~ x1, data = data),
    "The number of variables has to be greater than 1"
  )

  expect.test.error(
    bvs.lm(y ~ x1 + x2, data = data, n.keep = "two"),
    "Only n.keep='all' or type the exact number"
  )

  expect.test.error(
    bvs.lm(y ~ x1 + x2, data = data, n.keep = "all", prior.models = "SB"),
    "Only priors 'ScottBerger', 'Constant' and 'User' supported"
  )

  expect.test.error(
    bvs.lm(y ~ x1 + x2, data = data, n.keep = "all",
                  BF.method = "gprior", prior.betas = "gprior"),
    "prior.betas must be one of 'gZellner', 'Robust', "
  )

  expect.test.error(
    bvs.lm(y ~ x1 + x2, data = data, null.model = y ~ -1),
    "Null model shall contain the intercept term"
  )

  data <- lm.test.data(missing = TRUE)
  output <- capture.output({
    result <- bvs.lm(
      y ~ x1 + x2,
      data = data,
      n.keep = 5,
      n.imp = 39,
      imp.seed = 9011975
    )

    print(result)
  })

  expect_equal(output[1],
               "The number of models to keep (5) is larger than the total number of models (4) and it has been set to 4.")
  expect_equal(output[23],
               "The 4 most probable models and their probabilities are:")
  expect_equal(trimws(output[16]),
               "our dearest *missing* <3")
})
