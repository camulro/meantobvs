#Tests for the gibbsbvs.lm function

#Check that the result is of the desired form and basic outputs without factors
test_that("gibbsbvsGD25 returns Gibbs summaries", {
  data <- lm.test.data()

  output <- capture.output({
    result <- suppressWarnings(gibbsbvsGD25(
      y ~ x1 + x2,
      data = data,
      init.model = "Null",
      n.iter = 12,
      n.burnin = 3,
      n.thin = 2,
      n.imp = 2,
      Gibbs.seed = 123,
      imp.seed = 1
    ))
  })

  expect_gt(length(output), 0) #output is not empty (greater than 0)
  expect_s3_class(result, "meantobvs") #checks class(result)
  expect_s3_class(result$lmnull, "lm") #checks class(result$lmnull)
  expect_s3_class(result$lmfull, "mipo") #checks class(result$lmfull)

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
  expect_equal(output[6], "x1, x2")
  expect_equal(output[7],
               "The problem has a total of 4 competing models.")
  expect_equal(output[8],
               "Of these, 15 are sampled with replacement.")
  expect_equal(output[9],
               "Then, 6 are kept and used to construct the summaries.")


  expect_equal(result$imp.info$n.imp, 2)
  expect_equal(result$imp.info$imp.seed, 1)

  expect_false(anyNA(result$modelslogBF)) #check there are no NAs
  expect_false(anyNA(result))
  expect_true(all(is.finite(result$modelslogBF)))
  expect_true(is.finite(result$C))
  expect_error(capture.output(print(result)), NA)
  expect_error(capture.output(summary(result)), NA)
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

#Check that errors are given correctly and other printed messages
test_that("gibbsbvsGD25 reports errors and warnings", {
  data <- lm.test.data()

  expect.test.error(
    suppressWarnings(gibbsbvsGD25(
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
  data$x3 <- as.factor(x3)

  expect.test.error(
    gibbsbvsGD25(y ~ x1 + x3, data = data),
    "This method is only for continuous covariates"
  )

  data <- data[,-4] #remove x3
  expect.test.error(suppressWarnings(
    gibbsbvsGD25(y ~ x1, data = data)),
    "The number of variables has to be greater than 1"
  )

  output <- capture.output({
    result <- suppressWarnings(
      gibbsbvsGD25(
        y ~ x1 + x2,
        data = data,
        n.iter = 10,
        Gibbs.seed = 123
      ))

    print(result)
  })

  expect_equal(output[9],
               "Then, 10 are kept and used to construct the summaries.")
  expect_equal(output[17],
               "The  1  most probable models among the visited ones are:")

  expect_warning(
    capture.output(gibbsbvsGD25(y ~ x1 + x2, data = data, n.iter = 2)),
    "The number of variables is small enough"
  )
})

#just one imputation
test_that("missingGibbsGD25 accepts a single imputation", {
  data <- lm.test.data()

  capture.output({
    result <- suppressWarnings(gibbsbvsGD25(
      y ~ x1 + x2,
      data = data,
      init.model = "Null",
      n.iter = 4, n.burnin = 0, n.thin = 1,
      n.imp = 1,
      Gibbs.seed = 1, imp.seed = 1
    ))
  })
  expect_equal(dim(result$modelslogBF), c(4, 3))
  expect_true(all(is.finite(result$modelslogBF)))

})
