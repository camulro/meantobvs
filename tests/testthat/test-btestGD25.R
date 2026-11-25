#Tests for the btestGD25 function

#Check that the result is of the desired form and basic outputs
test_that("btestGD25 computes posterior probabilities for named models", {

  #build test deterministic data
  data <- lm.test.data()

  models <- list(null = y ~ 1, x1 = y ~ x1, full = y ~ x1 + x2)

  #capture all text written by the function on console
  output <- capture.output({
    result <- btestGD25(
      data = data,
      models = models,
      n.imp = 2,
      imp.seed = 1
    )
  })

  expect_gt(length(output), 0) #output is not empty (greater than 0)

  expect_s3_class(result, "meantobtest") #checks class(result)

  expect_identical(result$nullmodel, "null")
  expect_identical(names(result$lBFi0), c("null.to.null", "x1.to.null", "full.to.null"))
  expect_identical(result$lBFi0[["null.to.null"]], 0)
  expect_identical(result$models, models)

  #the Constant prior is the default: 1/N per model
  expect_identical(result$prior.models, "Constant")
  expect_equal(unname(result$priorprobs), rep(1 / 3, 3))

  expect_equal(sum(result$PostProbi), 1) #check probs
  expect_equal(sum(result$priorprobs), 1)
  expected.PostProbi <- exp(result$lBFi0) * result$priorprobs /
    sum(exp(result$lBFi0) * result$priorprobs)
  #check computed inclprobs
  expect_equal(unname(result$PostProbi), unname(expected.PostProbi))

  expect_equal(result$imp.info$n.imp, 2)
  expect_equal(result$imp.info$imp.seed, 1)

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
test_that("btestGD25 correctly performs with other prior methods", {
  #build test deterministic data
  data <- lm.test.data()

  models <- list(x1 = y ~ x1, full = y ~ x1 + x2)

  output <- capture.output({
    result <- btestGD25(
      data = data,
      models = models,
      prior.models = "User",
      priorprobs = c(2, 1, 1),
      n.imp = 2,
      imp.seed = 1
    )
  })

  #check that null model is added to the list of competing models
  expect_true(any(grepl("added to the list of competing models", output)))
  expect_identical(result$nullmodel, "(y ~ 1)")
  expect_identical(names(result$PostProbi), c("x1", "full", "(y ~ 1)"))
  expect_identical(names(result$lBFi0),
                   c("x1.to.((y ~ 1))", "full.to.((y ~ 1))", "(y ~ 1).to.((y ~ 1))"))
  expect_identical(result$lBFi0[["(y ~ 1).to.((y ~ 1))"]], 0)
  expect_identical(names(result$modelspool), c("x1", "full", "(y ~ 1)"))

  expect_gt(length(output), 0) #output is not empty (greater than 0)

  expect_identical(result$prior.models, "User")

  expect_equal(sum(result$PostProbi), 1) #check probs
  expect_equal(sum(result$priorprobs), 1)
  expected.PostProbi <- exp(result$lBFi0) * result$priorprobs /
    sum(exp(result$lBFi0) * result$priorprobs)
  #check computed inclprobs
  expect_equal(unname(result$PostProbi), unname(expected.PostProbi))

  expect_equal(result$imp.info$n.imp, 2)
  expect_equal(result$imp.info$imp.seed, 1)
})

#Check that errors are given correctly and other printed messages
test_that("btestGD25 reports errors and warnings", {
  data <- lm.test.data()

  expect.test.error(
    btestGD25(data = data, models = c("M0", "M1")),
    "Argument models should be a list."
  )

  expect.test.error(
    btestGD25(data = data, models = list(y ~ 1)),
    "At least two models have to be provided."
  )

  models <- list(null = y ~ 1, full = y ~ x1 + x2)

  expect.test.error(
    btestGD25(data = data, models = models, prior.models = "SB"),
    "Only priors 'ScottBerger', 'Constant' and 'User' supported."
  )

  output <- capture.output({
    result <- btestGD25(
      data = data,
      models = models
    )

    print(result)
  })

  expect_equal(output[14],
               "log Bayes factors (expressed in relation to null)")

  #user model prior
  models <- list(null = y ~ 1, x1 = y ~ x1, full = y ~ x1 + x2)

  expect.test.error(
    btestGD25(data = data, models = models, prior.models = "User",
              n.imp = 1, imp.seed = 1),
    "no prior probabilities provided"
  )
  expect.test.error(
    btestGD25(data = data, models = models, prior.models = "User",
              priorprobs = c(0.5, 0.5), n.imp = 1, imp.seed = 1),
    "incorrect length"
  )
  expect.test.error(
    btestGD25(data = data, models = models, prior.models = "User",
              priorprobs = c(1, -1, 1), n.imp = 1, imp.seed = 1),
    "must be positive"
  )
  expect.test.error(
    btestGD25(data = data, models = models, prior.models = "User",
              priorprobs = c(0, 0, 0), n.imp = 1, imp.seed = 1),
    "must be positive"
  )
  expect.test.error(
    btestGD25(data = data, models = models, prior.models = "uniform",
              n.imp = 1, imp.seed = 1),
    "Only priors 'ScottBerger', 'Constant' and 'User' supported"
  )

  #Build a noise variable
  withr::with_seed(1, #set seed for reproducibility
                   x3 <- rnorm(nrow(data))
  )
  data$x3 <- as.factor(x3)

  expect.test.error(
    btestGD25(list(y ~ 1, full = y ~ x1 + x3), data = data),
    "This method is only for continuous covariates"
  )

})
