#Tests for the bvsGD25 function

#Check that the result is of the desired form and basic outputs
test_that("bvsGD25 returns the full posterior for complete data", {
  #build test deterministic data
  data <- lm.test.data()

  output <- capture.output({
    result <- bvsGD25(y ~ x1 + x2,
                      data = data,
                      n.keep = "all",
                      n.imp = 2,
                      imp.seed = 1)
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
  expect_equal(output[6], "x1, x2")
  expect_equal(output[7],
               "The problem has a total of 4 competing models.")
  expect_equal(output[8],
               "Of these, the  4 most probable (a posteriori) are kept.")
  expect_equal(result$imp.info$n.imp, 2)
  expect_equal(result$imp.info$imp.seed, 1)

  expect_false(anyNA(result$modelsprob)) #check there are no NAs
  expect_false(anyNA(result))
  expect_true(all(is.finite(result$modelsprob)))
  expect_true(is.finite(result$C))
})


#Check other options
test_that("bvsGD25 correctly performs with other prior methods", {
  #build test deterministic data
  data <- lm.test.data()

  #capture all text written by the function on console
  output <- capture.output({
    result <- bvsGD25(y ~ x1 + x2,
                      data = data,
                      prior.models = "User",
                      priorprobs = c(2, 1, 1),
                      n.keep = 1,
                      n.imp = 2,
                      imp.seed = 1)
  })

  expect_gt(length(output), 0) #output is not empty (greater than 0)

  expect_identical(rownames(result$modelsprob), as.character(1))
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

  expect_equal(result$imp.info$n.imp, 2)
  expect_equal(result$imp.info$imp.seed, 1)
})

#Check that errors are given correctly and other printed messages
test_that("bvsGD25 reports errors and warnings", {
  data <- lm.test.data()
  #Build a noise variable
  withr::with_seed(1, #set seed for reproducibility
                   x3 <- rnorm(nrow(data))
  )
  data$x3 <- as.factor(x3)

  expect.test.error(
    bvsGD25(y ~ x1 + x3, data = data),
    "This method is only for continuous covariates"
  )

  data <- data[,-4] #remove x3
  expect.test.error(
    bvsGD25(y ~ x1, data = data),
    "The number of variables has to be greater than 1"
  )

  expect.test.error(
    bvsGD25(y ~ x1 + x2, data = data, n.keep = "two"),
    "Only n.keep='all' or type the exact number"
  )

  expect.test.error(
    bvsGD25(y ~ x1 + x2, data = data, n.keep = "all", prior.models = "SB"),
    "Only priors 'ScottBerger', 'Constant' and 'User' supported"
  )

  expect.test.error(
    bvsGD25(y ~ x1 + x2, data = data, prior.models = "User",
                n.imp = 1, imp.seed = 1),
    "no prior probabilities provided"
  )

  expect.test.error(
    bvsGD25(y ~ x1 + x2, data = data, prior.models = "User",
                priorprobs = c(1, 2), n.imp = 1, imp.seed = 1),
    "length of prior"
  )

  #accepts single imputations?
  capture.output({
    result1 <- bvsGD25(y ~ x1 + x2, data = data,
                       n.keep = "all", n.imp = 1, imp.seed = 1)
  })
  expect_equal(sum(result1$modelsprob[, "Post"]), 1)
  expect_true(is.finite(result1$C))

  output <- capture.output({
    result <- bvsGD25(
      y ~ x1 + x2,
      data = data,
      n.keep = 5
    )

    print(result)
  })

  expect_equal(output[1],
               "The number of models to keep (5) is larger than the total number of models (4) and it has been set to 4.")
  expect_equal(output[16],
               "The 4 most probable models and their probabilities are:")
})

