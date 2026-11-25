#Tests for the gibbsbvs.lm function

#Check that the result is of the desired form and basic outputs
test_that("gibbsbvs.lm returns Gibbs summaries", {
  #build test deterministic data
  data <- lm.test.data()

  #capture all text written by the function on console
  output <- capture.output({
    result <- suppressWarnings(
      gibbsbvs.lm(
        y ~ x1 + x2,
        data = data,
        BF.method = "BIC",
        prior.betas = "gZellner", #does not apply with BIC
        init.model = "Null",
        n.iter = 12,
        n.burnin = 3,
        n.thin = 2,
        Gibbs.seed = 123
      )
    )
  })

  expect_gt(length(output), 0) #output is not empty (greater than 0)
  expect_s3_class(result, "meantobvs") #checks class(result)
  expect_s3_class(result$lmnull, "lm") #checks class(result$lmnull)
  expect_s3_class(result$lmfull, "lm") #checks class(result$lmfull)

  expect_identical(result$method, "Gibbs") #Check names
  expect_identical(result$variables, c("x1", "x2"))
  expect_identical(names(result$inclprob), c("x1", "x2"))
  expect_identical(names(result$inclprobRB), c("x1", "x2"))
  expect_identical(names(result$postprobdim), as.character(1:3))
  expect_identical(colnames(result$modelslogBF), c("x1", "x2", "logBF"))

  expect_equal(nrow(result$modelslogBF), 6)
  expect_equal(ncol(result$modelslogBF), 3)
  expect_equal(length(result$prior.models), 1)
  expect_equal(result$p, 2)
  expect_equal(result$n, nrow(data))
  expect_equal(result$k, 1) #Only intercept fixed
  expect_equal(sum(result$priorprobs), 1) #check probs
  expect_equal(length(result$priorprobs), 2 + 1)
  expect_equal(sum(result$postprobdim), 1)
  expect_equal(length(result$postprobdim), 2 + 1)
  expect_equal(length(result$HPMbin), 2) #check summaries
  expect_equal(length(result$MPMbin), 2)
  expect_equal(output[5],
               "  Competing variables: x1 x2 ")
  expect_equal(output[6],
               "The problem has a total of 4 competing models.")
  expect_equal(output[7],
               "Of these, 15 are sampled with replacement.")
  expect_equal(output[8],
               "Then, 6 are kept and used to construct the summaries.")


  expect_null(result$prior.betas)
  expect_null(result$imp.info)

  expect_false(anyNA(result$modelslogBF)) #check there are no NAs
  expect_false(anyNA(result))
  expect_true(all(is.finite(result$modelslogBF)))
  expect_true(is.finite(result$C))
})

#Check other options
test_that("gibbsbvs.lm correctly performs other prior methods", {
  #build test deterministic data
  data <- lm.test.data()

  #capture all text written by the function on console
  output <- capture.output({
    result <- suppressWarnings(
      gibbsbvs.lm(
        y ~ x1 + x2,
        data = data,
        BF.method = "gprior",
        prior.betas = "gZellner",
        prior.models = "User",
        priorprobs = c(2, 1, 1),
        init.model = c(0, 1),
        n.iter = 1,
        n.burnin = 0,
        Gibbs.seed = 123
      ))
  })

  expect_gt(length(output), 0) #output is not empty (greater than 0)

  expect_identical(rownames(result$modelslogBF), as.character(1))
  expect_identical(result$prior.betas, "gZellner")
  expect_identical(result$prior.models, "User")

  expect_equal(result$priorprobs, c(`1` = 1/2, `2` = 1/4, `3` = 1/4))
  expect_equal(sum(result$priorprobs), 1) #check probs
  expect_equal(length(result$priorprobs), 2 + 1)
  expect_equal(sum(result$postprobdim), 1)
  expect_equal(length(result$postprobdim), 2 + 1)
  expect_equal(length(result$inclprob), 2)
  expect_equal(length(result$inclprobRB), 2)
  expect_equal(result$postprobs, NA_real_) #no posterior computed with 1 iter
  expect_equal(result$C, NA_real_)

  expect_null(result$imp.info)
})

#Check how external imputed data is treated
test_that("gibbsbvs.lm uses supplied imputations", {
  #Build test data with missings and impute with the mean
  data <- lm.test.data(missing = TRUE)
  imputations <- lm.test.imputations(data)

  output <- capture.output({
    result <- suppressWarnings(
      gibbsbvs.lm(
        y ~ x1 + x2,
        data = data,
        imp.datasets = imputations,
        n.imp = 3,
        parallelmice = TRUE, #to check imp.info
        n.iter = 12,
        n.burnin = 0,
        n.thin = 1,
        Gibbs.seed = 123
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
  expect_equal(nrow(result$modelslogBF), 12)
  expect_equal(sum(result$postprobdim), 1)

  expect_null(result$imp.info$parallelmice)
})


#Check that errors are given correctly and other printed messages
test_that("gibbsbvs.lm reports errors and warnings", {
  data <- lm.test.data()

  expect.test.error(
    suppressWarnings(gibbsbvs.lm(
      y ~ x1 + x2,
      data = data,
      init.model = c(1, 0, 1),
    )),
    "Initial model with incorrect length"
  )

  #Build a noise variable
  withr::with_seed(1, #set seed for reproducibility
                   x3 <- rnorm(nrow(data))
  )
  data$x3 <- x3

  expect.test.error(suppressWarnings(
    gibbsbvs.lm(y ~ x1, data = data, null.model = z ~ 1)),
    "response in the full and null model does not coincide"
  )

  expect.test.error(suppressWarnings(
    gibbsbvs.lm(y ~ x1 + x2, data = data, null.model = y ~ x3)),
    "Error in var: x3; null model not nested in full"
  )

  expect.test.error(suppressWarnings(
    gibbsbvs.lm(y ~ x1, data = data, null.model = y ~ x1)),
    "number of fixed variables is equal to the number of"
  )

  expect.test.error(suppressWarnings(
    gibbsbvs.lm(y ~ x1, data = data)),
    "The number of variables has to be greater than 1"
  )

  expect.test.error(suppressWarnings(
    gibbsbvs.lm(y ~ x1 + x2, data = data, null.model = y ~ -1)),
    "Null model shall contain the intercept term"
  )

  output <- capture.output({
    result <- suppressWarnings(
      gibbsbvs.lm(
        y ~ x1 + x2,
        data = data,
        n.iter = 10,
        Gibbs.seed = 123
      ))

    print(result)
  })

  expect_equal(output[8],
               "Then, 10 are kept and used to construct the summaries.")
  expect_equal(output[16],
               "The  1  most probable models among the visited ones are:")

  expect_warning(
    capture.output(gibbsbvs.lm(y ~ x1 + x2, data = data, n.iter = 2)),
    "The number of variables is small enough"
  )
})

