#Tests for factor related arguments of meantobvs: marginal.factors argument,
#factor detection, ordered levels, interaction terms, model.rankdef.matrix,
#no representation of saturated or oversaturated models.

#Test bvs.lm factor performance
test_that("bvs.lm marginalizes factors by default", {
  data <- lm.testfactor.data()

  output <- capture.output({
    result <- bvs.lm(y ~ x1 + x2 + f1 + f2, data = data, n.keep = "all")
  })

  expect_gt(length(output), 0)
  expect_s3_class(result, "meantobvs") #checks class(result)
  expect_s3_class(result$lmnull, "lm") #checks class(result$lmnull)
  expect_s3_class(result$lmfull, "lm") #checks class(result$lmfull)

  expect_identical(result$variables, c("x1", "x2", "f1", "f2"))
  expect_identical(dim(result$positions), c(2L, 7L))
  expect_identical(rownames(result$positions), c("f1", "f2"))
  expect_identical(colnames(result$positions),
                   c("x1", "x2", "f1L1", "f1L2", "f1L3", "f2L1", "f2L2"))
  expect_identical(colnames(result$modelsprob), c("x1", "x2", "f1", "f2", "Post"))
  expect_identical(unname(result$positionsx), c(TRUE, TRUE, FALSE, FALSE))
  #the prior for the dummies is added to the one for the covariates/factors
  expect_identical(result$prior.models, c("ScottBerger", "ScottBerger"))
  expect_identical(output[6], "  Factors: f1 f2 ")

  expect_equal(result$p, 4) #covariates and factors
  expect_equal(result$k, 1)
  expect_equal(dim(result$modelsprob), c(16, 5))
  expect_equal(dim(result$modelsrankdefprob), c(2^2*(2^3-3)*(2^2-2), 8))
  expect_equal(sum(result$modelsprob[, "Post"]), 1) #check probs
  expect_equal(sum(result$priorprobs), 1)
  expect_equal(length(result$priorprobs), 4 + 1)
  expect_equal(sum(result$postprobdim), 1)
  expect_equal(length(result$postprobdim), 4 + 1)
  expect_equal(length(result$HPMbin), 4) #check summaries
  expect_equal(length(result$MPMbin), 4)
  expected.inclprob <- colSums(result$modelsprob[, c("x1", "x2", "f1", "f2")] *
                                 result$modelsprob[, "Post"])
  #check computed inclprobs
  expect_equal(unname(result$inclprob), unname(expected.inclprob))
  expect_all_equal(result$modelsrankdefprob[,"f2L1"], 0)
  mrd <- result$modelsrankdefprob[,c("f1L1", "f1L2", "f1L3")]
  for (i in 1:19) {
    #check representant used
    expect_false(identical(mrd[i,], c(1,1,1)) ||
                   identical(mrd[i,], c(1,0,1)) ||
                   identical(mrd[i,], c(1,1,0)))
  }

  expect_equal(output[5],
               "  Numerical covariates: x1 x2 ")
  expect_equal(output[7],
               "The problem has a total of 16 competing models.")
  expect_equal(output[8],
               "Of these, the  16 most probable (a posteriori) are kept.")

  #outputs only reported when there are factor(s) to marginalize
  expect_true(is.matrix(result$positions))
  expect_true(is.matrix(result$modelsrankdefprob))
  expect_true(result$marginal.factors)

  expect_false(anyNA(result$modelsprob)) #check there are no NAs
  expect_false(anyNA(result))
  expect_true(all(is.finite(result$modelsprob)))
  expect_true(is.finite(result$C))

})

#Test bvs.lm factor performance when no factor marginalization is asked
test_that("bvs.lm treats the factor dummies as covariates without marginalization", {
  data <- lm.testfactor.data()

  output <- capture.output({
    result <- bvs.lm(y ~ x1 + x2 + f1, data = data, n.keep = "all",
                            marginal.factors = FALSE)
  })

  expect_false(result$marginal.factors)

  expect_identical(result$variables, c("x1", "x2", "f1L2", "f1L3"))
  expect_identical(result$prior.models, "ScottBerger")
  expect_identical(colnames(result$modelsprob), c("x1", "x2", "f1L2", "f1L3", "Post"))

  expect_equal(result$p, 4) #two dummies of the factor enter as competing variables
  expect_equal(result$k, 1)
  expect_equal(dim(result$modelsprob), c(16, 5))
  expect_equal(sum(result$modelsprob[, "Post"]), 1)

  #no factor-specific outputs are reported?
  expect_null(result$positions)
  expect_null(result$positionsx)
  expect_null(result$modelsrankdefprob)

  expect_false(anyNA(result$modelsprob)) #check there are no NAs
  expect_false(anyNA(result))
  expect_true(all(is.finite(result$modelsprob)))
  expect_true(is.finite(result$C))

})

#Test factors prior probabilities computed
test_that("bvs.lm computes the prior over factor model space correctly", {
  data <- lm.testfactor.data()

  output <- capture.output({
    result <- bvs.lm(y ~ x1 + f1 + f2, data = data, n.keep = "all")
  })

  C <- result$C
  post <- result$modelsprob[,"Post"]
  prior <- result$priorprobs

  #just intercept, x1 or f2
  expect_equal(unname(post[6]), unname(result$modelsrankdefprob[20, "Post"]))
  expect_equal(unname(post[3]), unname(result$modelsrankdefprob[1, "Post"]))
  expect_equal(unname(post[8]), unname(result$modelsrankdefprob[10, "Post"]))
  #x1 and f2
  expect_equal(unname(post[4]), unname(result$modelsrankdefprob[11, "Post"]))

  BF <- post * C / prior[rowSums(result$modelsprob[,1:3]) + 1]
  #compute BFs:
  BF[5] #just f1

  mf <- model.frame(y ~ x1 + f1 + f2, data)
  lmnull <- lm(y ~ 1, data, x = T, y = T)
  lBF <- function(model) BF.BIC.lm(y = lmnull$y,
                                   X = model.rankdef.matrix(mf)[,c(1,model + 1)],
                                   SS0 = as.numeric(crossprod(lmnull$residuals)))

  submodels <- result$modelsrankdefprob[, 1:6]

  compBF <-
    exp(
      lBF(which(submodels[2,] == 1)) +
        logScottBerger.d(c(1,0), c(TRUE, FALSE), c(3,2))) +
      exp(lBF(which(submodels[4,] == 1)) +
        logScottBerger.d(c(1,0), c(TRUE, FALSE), c(3,2))) +
      exp(lBF(which(submodels[6,] == 1)) +
        logScottBerger.d(c(1,0), c(TRUE, FALSE), c(3,2))) +
      exp(lBF(which(submodels[8,] == 1)) +
        logScottBerger.d(c(2,0), c(TRUE, FALSE), c(3,2))
    )
  expect_equal(compBF, BF[5])






  #fijando un factor

  # expect_false(result$marginal.factors)
  #
  # expect_identical(result$variables, c("x1", "x2", "f1L2", "f1L3"))
  # expect_identical(result$prior.models, "ScottBerger")
  # expect_identical(colnames(result$modelsprob), c("x1", "x2", "f1L2", "f1L3", "Post"))
  #
  # expect_equal(result$p, 4) #two dummies of the factor enter as competing variables
  # expect_equal(result$k, 1)
  # expect_equal(dim(result$modelsprob), c(16, 5))
  # expect_equal(sum(result$modelsprob[, "Post"]), 1)
  #
  # #no factor-specific outputs are reported?
  # expect_null(result$positions)
  # expect_null(result$positionsx)
  # expect_null(result$modelsrankdefprob)
  #
  # expect_false(anyNA(result$modelsprob)) #check there are no NAs
  # expect_false(anyNA(result))
  # expect_true(all(is.finite(result$modelsprob)))
  # expect_true(is.finite(result$C))

})

#Test bvs.lm factor performance when interactions
test_that("bvs.lm marginalizes factors when factor interaction", {
  data <- lm.testfactor.data()

  # output <- capture.output({
  #   result <- bvs.lm(y ~ x1 + x2 + f1 * f2, data = data, n.keep = "all")
  # })

  expect_gt(length(output), 0)
  expect_s3_class(result, "meantobvs") #checks class(result)
  expect_s3_class(result$lmnull, "lm") #checks class(result$lmnull)
  expect_s3_class(result$lmfull, "lm") #checks class(result$lmfull)

  expect_identical(result$variables, c("x1", "x2", "f1", "f2"))
  expect_identical(dim(result$positions), c(2L, 7L))
  expect_identical(rownames(result$positions), c("f1", "f2"))
  expect_identical(colnames(result$positions),
                   c("x1", "x2", "f1L1", "f1L2", "f1L3", "f2L1", "f2L2"))
  expect_identical(colnames(result$modelsprob), c("x1", "x2", "f1", "f2", "Post"))
  expect_identical(unname(result$positionsx), c(TRUE, TRUE, FALSE, FALSE))
  #the prior for the dummies is added to the one for the covariates/factors
  expect_identical(result$prior.models, c("ScottBerger", "ScottBerger"))
  expect_identical(output[6], "  Factors: f1 f2 ")

  expect_equal(result$p, 4) #covariates and factors
  expect_equal(result$k, 1)
  expect_equal(sum(result$modelsprob[, "Post"]), 1) #check probs
  expect_equal(sum(result$priorprobs), 1)
  expect_equal(length(result$priorprobs), 4 + 1)
  expect_equal(sum(result$postprobdim), 1)
  expect_equal(length(result$postprobdim), 4 + 1)
  expect_equal(length(result$HPMbin), 4) #check summaries
  expect_equal(length(result$MPMbin), 4)
  expected.inclprob <- colSums(result$modelsprob[, c("x1", "x2", "f1", "f2")] *
                                 result$modelsprob[, "Post"])
  #check computed inclprobs
  expect_equal(unname(result$inclprob), unname(expected.inclprob))
  expect_equal(output[5],
               "  Numerical covariates: x1 x2 ")
  expect_equal(output[7],
               "The problem has a total of 16 competing models.")
  expect_equal(output[8],
               "Of these, the  16 most probable (a posteriori) are kept.")

  #outputs only reported when there are factor(s) to marginalize
  expect_true(is.matrix(result$positions))
  expect_true(is.matrix(result$modelsrankdefprob))
  expect_true(result$marginal.factors)

  expect_false(anyNA(result$modelsprob)) #check there are no NAs
  expect_false(anyNA(result))
  expect_true(all(is.finite(result$modelsprob)))
  expect_true(is.finite(result$C))

})

#Test errors related with factors
test_that("bvs.lm reports errors related with factors", {
  data <- lm.testfactor.data()

  ############################### bvs.lm(y ~ f1, data = data, n.keep = "all")

  expect.test.error(
    bvs.lm(y ~ x1 + x2 + f1, data = data, n.keep = "all",
                  prior.models.dummies = "User"),
    "Only priors 'ScottBerger' and 'Constant' supported"
  )

  output <- capture.output({
    result <- bvs.lm(y ~ x1 + x2 + f1, data = data, n.keep = "all",
                     prior.models = "User",
                     priorprobs = c(1, 2, 3, 4))
  })
  #priorprobs are renormalized to sum one in bvs.lm
  expect_equal(unname(result$priorprobs), c(0.1, 0.2, 0.3, 0.4))

  expect.test.error(
    bvs.lm(y ~ x1 + x2 + f1, data = data, n.keep = "all",
           prior.models = "User", priorprobs = c(0.2, 0.3)),
    "User prior selected but"
  )

  output <- capture.output({
    result <- bvs.lm(y ~ x1 + x2 + f1 + f2, data = data, null.model = y ~ f1,
                     n.keep = "all")
  })

  expect_identical(result$variables, c("x1", "x2", "f2"))
  expect_identical(result$prior.models, c("ScottBerger", "ScottBerger"))
  expect_equal(result$p, 3)
  expect_equal(result$k, 2) #intercept and fixed factor
  expect_equal(sum(result$modelsprob[, "Post"]), 1)

  expect_true(is.matrix(result$positions))

  expect_false(anyNA(result$modelsprob)) #check there are no NAs
  expect_false(anyNA(result))
  expect_true(all(is.finite(result$modelsprob)))
  expect_true(is.finite(result$C))
})

#Factor imputation performed
test_that("bvs.lm imputes missing factors", {
  data <- lm.testfactor.data(missing = TRUE)

  output <- capture.output({
    result <- bvs.lm(y ~ x1 + x2 + f1, data = data, n.imp = 2,
                            imp.seed = 1, n.keep = "all")
  })

  expect_identical(result$variables, c("x1", "x2", "f1"))
  expect_true(inherits(result$lmfull, "mipo"))
  expect_true(is.matrix(result$positions))

  expect_equal(result$imp.info$n.imp, 2)
  expect_equal(unname(result$imp.info$NAvars), 3:5)
  expect_equal(as.vector(result$lmfull$pooled$term),
               c("(Intercept)", "x1", "x2", "f1L2", "f1L3"))
  expect_equal(sum(result$modelsprob[, "Post"]), 1)

  expect_true("compress.imp.array" %in% names(result))
  #recover imputed data used
  imp.array <- unserialize(
    memDecompress(
      result$compress.imp.array, type = "xz")
  )

  expect_equal(dimnames(imp.array)[[2]],
               c("(Intercept)", "x1", "x2", "f1L1", "f1L2", "f1L3"))
  expect_all_true(rowSums(imp.array[,c("f1L1", "f1L2", "f1L3"),1]) == 1)

})

#Test model.rankdef.matrix
test_that("model.rankdef.matrix returns rank defficient matrices when factors", {
  data <- lm.testfactor.data(missing = TRUE)

  formula <- y ~ x1 + x2 + f1
  mf <- model.frame(as.formula(formula), data, na.action = NULL)
  mrf <- meantobvs:::model.rankdef.matrix(mf)

  expect_equal(colnames(mrf),
               c("(Intercept)", "x1", "x2", "f1L1", "f1L2", "f1L3"))
  expect_all_true(as.vector(is.na(mrf[c(4,11), c("f1L1", "f1L2","f1L3")])))

  })
