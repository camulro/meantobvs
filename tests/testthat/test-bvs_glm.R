#Tests for the bvs.glm function

#Check that the result is of the desired form and basic outputs without factors
test_that("bvs.glm returns the full posterior for complete data", {
  #go through all families available
  families <- glm.families() #with default links
  for (i in names(families)) {
    family.obj <- families[[i]]
    #build test deterministic data for family i
    data <- glm.family.data(i)

    #capture all text written by the function on console
    output <- capture.output({
      result <- bvs.glm(
        y ~ x1 + x2,
        data = data,
        family = family.obj,
        BF.method = "BIC",
        prior.betas = "gZellner", #does not apply with BIC
        n.keep = "all"
      )
    })

    expect_gt(length(output), 0) #output is not empty (greater than 0)
    expect_s3_class(result, "meantobvs") #checks class(result)
    expect_s3_class(result$glmnull, "glm") #checks class(result$lmnull)
    expect_s3_class(result$glmfull, "glm") #checks class(result$lmfull)

    expect_identical(result$method, "Full") #Check names
    expect_identical(result$variables, c("x1", "x2"))
    expect_identical(names(result$inclprob), c("x1", "x2"))
    expect_identical(names(result$postprobdim), as.character(1:3))
    expect_identical(colnames(result$modelsprob), c("x1", "x2", "Post"))
    expect_identical(rownames(result$modelsprob), as.character(1:4))
    expect_identical(result$glmnull$family$family, family.obj$family)

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
    expect_true(result$glmfull$converged)
  }
})

#Check other options
test_that("bvs.glm correctly performs with other prior methods", {
  #go through all families available
  families <- glm.families() #with default links
  for (i in names(families)) {
    family.obj <- families[[i]]
    #build test deterministic data for family i
    data <- glm.family.data(i)

    output <- capture.output({
      result <- bvs.glm(
        y ~ x1 + x2,
        data = data,
        family = family.obj,
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
    #   result$modelsprob[, "Post"]
    # #check computed inclprobs
    # expect_equal(unname(result$inclprob), unname(expected.inclprob))

    expect_true(result$glmfull$converged)

    expect_null(result$imp.info)
  }
})

#Check how external imputed data is treated
test_that("bvs.glm uses supplied imputations", {
  #go through all families available
  families <- glm.families() #with default links
  for (i in names(families)) {
    family.obj <- families[[i]]
    #build test deterministic data for family i
    data <- glm.family.data(i, missing = TRUE)
    imputations <- glm.test.imputations(data)

    output <- capture.output({
      result <- suppressWarnings(bvs.glm(
        y ~ x1 + x2,
        data = data,
        family = family.obj,
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
  }
})

#Check that probabilities computed are the same as other packages
#just for families available on both
test_that("bvs.glm's posterior is correct", {

  #Check all BF available
  pb <- glm.prior.betas()[-4] #ZellnerSiow is wrong?

  #go through all families available on both packages
  families <- glm.bas.families() #with default links

  #Build a noise variable
  withr::with_seed(1, #set seed for reproducibility
                   x3 <- rnorm(30)
  )

  for (i in names(families)) {
    family <- family.obj<- families[[i]]
    #build test deterministic data for family i
    data <- glm.family.data(i)
    data$x3 <- x3

    print(i)

    #just intercept fixed
    for (j in pb) {

      result.bvs <- bvs.glm(
        y ~ x1 + x2,
        data = data,
        family = family.obj,
        BF.method = "gprior",
        prior.betas = j,
        n.keep = 4
      )

      result.other <- BAS::bas.glm(
        y ~ x1 + x2,
        data = data,
        family = family,
        betaprior = prior.betas.gprior.glm(j, n = nrow(data), p = 2),
        method = "deterministic",
        n.models = 4
      )

      expect_equal(nrow(result.bvs$modelsprob), 4)
      expect_equal(sum(result.bvs$modelsprob[, "Post"]), 1)

      #compute postprobs in the same order to compare
      other.models <- result.bvs$modelsprob[,1:2]
      other.postprobs <- numeric(4)
      for (r in 1:4) { #r is each row on other.models
        aux <- which(other.models[r,] == 1)

        for (l in 1:4) {
          if (identical(result.other$which[[l]][-1],
                        as.integer(unname(aux)))) {
            other.postprobs[r] <- result.other$postprobs[l]
          }
        }
      }

      expect_equal(unname(result.bvs$modelsprob[, "Post"]),
                   other.postprobs, tolerance = 10^(-5))
      expect_equal(unname(result.bvs$inclprob),
                   result.other$probne0[-1], tolerance = 10^(-5))

      modeldim <- rowSums(other.models)
      other.postprobdim <- sapply(0:2, function (k)
        sum(other.postprobs[which(modeldim == k)]))
      expect_equal(unname(result.bvs$postprobdim),
                   other.postprobdim, tolerance = 10^(-5))

      other.HPM <- other.models[which.max(other.postprobs),]
      expect_equal(result.bvs$HPMbin, other.HPM)
      other.MPM <- which(result.other$probne0[-1] > 0.5)
      expect_identical(which(unname(result.bvs$inclprob) > 0.5),
                       other.MPM)

      expect_true(result.bvs$glmfull$converged)
    }

    #x1 fixed
    for (j in pb) {

      result.bvs <- bvs.glm(
        y ~ x1 + x2 + x3,
        null.model = y ~ x1,
        data = data,
        family = family.obj,
        BF.method = "gprior",
        prior.betas = j,
        n.keep = 4
      )

      ##To check particular BF:
      # #Build design matrix
      # X <- lm.test.complete.design(data)[, c("(Intercept)", "x1")]
      # X2 <- lm.test.complete.design(data)[, c("(Intercept)", "x1", "x2")]
      #
      # #Fit models
      # null.fit <- glm(y ~ 1, data = data, family = family.obj, y = T, x = T)
      # model.fit <- glm.fit(x = X, y = data$y, family = family.obj)
      # model.fit2 <- glm.fit(x = X2, y = data$y, family = family.obj)
      #
      # #Get arguments to compute BF approx
      # n <- nrow(data)
      # k <- 1L #number of competing variables
      #
      # args <- prior.betas.gprior.glm.fit(j, n = n, p = 3)
      # if (args$type == "fixed") {
      #   result2.1 <- meantobvs:::BF.gprior.glm.fit(
      #     y = data$y,
      #     X = X2,
      #     glm0 = null.fit,
      #     g = args$g, n = n, k = k + 1L
      #   )
      #
      #   result2.2 <- meantobvs:::BF.gprior.glm.fit(
      #     y = data$y,
      #     X = X,
      #     glm0 = null.fit,
      #     g = args$g, n = n, k = k
      #   )
      #
      #   result2 <- result2.1 - result2.2
      #
      # } else {
      #   result2.1 <- meantobvs:::BF.hypergprior.glm.fit(
      #     y = data$y,
      #     X = X2,
      #     glm0 = null.fit,
      #     prior.betas.args = args,
      #     n = n, k = k + 1L
      #   )
      #
      #   result2.2 <- meantobvs:::BF.hypergprior.glm.fit(
      #     y = data$y,
      #     X = X,
      #     glm0 = null.fit,
      #     prior.betas.args = args,
      #     n = n, k = k
      #   )
      #
      #   result2 <- result2.1 - result2.2
      # }
      # result2
      # result.bvs$modelsprob[,"Post"] * result.bvs$C /(1/6)

      result.other <- BAS::bas.glm(
        y ~ x1 + x2 + x3,
        include.always = ~ x1,
        data = data,
        family = family,
        betaprior = prior.betas.gprior.glm(j, n = nrow(data), p = 3),
        method = "deterministic",
        n.models = 4
      )

      expect_equal(nrow(result.bvs$modelsprob), 4)
      expect_equal(sum(result.bvs$modelsprob[, "Post"]), 1)

      #compute postprobs in the same order to compare
      other.models <- result.bvs$modelsprob[,1:2]
      other.postprobs <- numeric(4)
      for (r in 1:4) { #r is each row on other.models
        aux <- which(other.models[r,] == 1)

        for (l in 1:4) {
          if (identical(result.other$which[[l]][c(-1,-2)],
                        as.integer(unname(aux) + 1))) {
            other.postprobs[r] <- result.other$postprobs[l]
          }
        }
      }

      expect_equal(unname(result.bvs$modelsprob[, "Post"]),
                   other.postprobs, tolerance = 10^(-4))
      expect_equal(unname(result.bvs$inclprob),
                   result.other$probne0[c(-1,-2)], tolerance = 10^(-5))

      modeldim <- rowSums(other.models)
      other.postprobdim <- sapply(0:2, function (k)
        sum(other.postprobs[which(modeldim == k)]))
      expect_equal(unname(result.bvs$postprobdim),
                   other.postprobdim, tolerance = 10^(-5))

      other.HPM <- other.models[which.max(other.postprobs),]
      expect_equal(result.bvs$HPMbin, other.HPM)
      other.MPM <- which(result.other$probne0[c(-1,-2)] > 0.5)
      expect_identical(which(unname(result.bvs$inclprob) > 0.5),
                       other.MPM)

      expect_true(result.bvs$glmfull$converged)
    }
  }
})

#Check that errors are given correctly and other printed messages
test_that("bvs.glm reports errors and warnings", {
  data <- glm.test.data()
  #Build a noise variable
  withr::with_seed(1, #set seed for reproducibility
                   x3 <- rnorm(nrow(data))
  )
  data$x3 <- x3

  expect.test.error(
    bvs.glm(y ~ x1, data = data, family = binomial(), null.model = z ~ 1),
    "response in the full and null model does not coincide"
  )

  expect.test.error(
    bvs.glm(y ~ x1 + x2, data = data, family = binomial(), null.model = y ~ x3),
    "Error in var: x3; null model not nested in full"
  )

  expect.test.error(
    bvs.glm(y ~ x1, data = data, family = binomial(), null.model = y ~ x1),
    "number of fixed variables is equal to the number of"
  )

  expect.test.error(
    bvs.glm(y ~ x1, data = data, family = binomial()),
    "The number of variables has to be greater than 1"
  )

  expect.test.error(
    bvs.glm(y ~ x1 + x2, data = data, family = binomial(), n.keep = "two"),
    "Only n.keep='all' or type the exact number"
  )

  expect.test.error(
    bvs.glm(y ~ x1 + x2, data = data, family = binomial(),
                   n.keep = "all", prior.models = "SB"),
    "Only priors 'ScottBerger', 'Constant' and 'User' supported"
  )

  expect.test.error(
    bvs.glm(y ~ x1 + x2, data = data, family = binomial(), n.keep = "all",
                  BF.method = "gprior", prior.betas = "gprior"),
    "Prior.betas must be one of 'gZellner', 'Robust', "
  )

  expect.test.error(
    bvs.glm(y ~ x1 + x2, data = data, null.model = y ~ -1),
    "Null model shall contain the intercept term"
  )

  output <- capture.output({
    result <- bvs.glm(
      y ~ x1 + x2,
      data = data,
      family = binomial(),
      n.keep = 5
    )

    print(result)
  })

  expect_equal(output[1],
               "The number of models to keep (5) is larger than the total number of models (4) and it has been set to 4.")
  expect_equal(output[16],
               "The 4 most probable models and their probabilities are:")

  #family-link warnings:
  expect.test.error(
    bvs.glm(y ~ x1 + x2, data = data, family = quasi()),
    "family quasi non-recognized"
  )
  expect.test.error(
    bvs.glm(y ~ x1 + x2, data = data, family = binomial("inverse")),
    "link inverse non-recognized"
  )
})

