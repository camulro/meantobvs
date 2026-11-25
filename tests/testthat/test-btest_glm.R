#Tests for the btest.glm function

#Check that the result is of the desired form and basic outputs without factors
test_that("btest.glm computes posterior probabilities for named models", {

  #list of competing models
  models <- list(null = y ~ 1, x1 = y ~ x1, full = y ~ x1 + x2)

  #go through all families available
  families <- glm.families() #with default links
  for (i in names(families)) {
    family.obj <- families[[i]]
    #build test deterministic data for family i
    data <- glm.family.data(i)

    #capture all text written by the function on console
    output <- capture.output({
      result <- btest.glm(
        data = data,
        models = models,
        family = family.obj,
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

    expect_true(all(vapply(result$modelspool, inherits, logical(1), "glm"))) #checks class
    expect_true(all(is.finite(result$lBFi0)))
    expect_false(anyNA(result$lBFi0)) #check there are no NAs
    expect_false(anyNA(result))
  }
})


#Check other options
test_that("btest.glm correctly performs with other prior methods", {
  #list of competing models
  models <- list(null = y ~ 1, x1 = y ~ x1, full = y ~ x1 + x2)

  #go through all families available
  families <- glm.families() #with default links
  for (i in names(families)) {
    family.obj <- families[[i]]
    #build test deterministic data for family i
    data <- glm.family.data(i)

    #capture all text written by the function on console
    output <- capture.output({
      result <- btest.glm(
        data = data,
        models = models,
        family = family.obj,
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
  }
})

#Check how external imputed data is treated
test_that("btest.glm uses supplied imputations", {

  #list of competing models
  models <- list(null = y ~ 1, x1 = y ~ x1, full = y ~ x1 + x2)

  #go through all families available
  families <- glm.families() #with default links
  for (i in names(families)) {

    family.obj <- families[[i]]
    #build test deterministic data for family i
    data <- glm.family.data(i, missing = TRUE)
    imputations <- glm.test.imputations(data)

    #capture all text written by the function on console
    output <- capture.output({
      result <- btest.glm(
        data = data,
        models = models,
        family = family.obj,
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
  }
})

#Check that probabilities are computed correctly
test_that("btest.glm's posterior is correct", {

  #Full enumeration of models to compare with bvs output
  models <- list(null = y ~ 1, x1 = y ~ x1, x2 = y ~ x2, full = y ~ x1 + x2)
  #Check all BF available
  pb <- glm.prior.betas()

  #go through all families available
  families <- glm.families() #with default links
  for (i in names(families)) {
    family.obj <- families[[i]]
    #build test deterministic data for family i
    data <- glm.family.data(i)

    for (j in pb) {

      result.bvs <- bvs.glm(
        y ~ x1 + x2,
        data = data,
        family = family.obj,
        BF.method = "gprior",
        prior.betas = j,
        prior.models = "Constant"
      )

      result.btest <- btest.glm(
        data = data,
        models = models,
        family = family.obj,
        BF.method = "gprior",
        prior.betas = j
      )

      expect_identical(result.bvs$prior.betas, result.btest$prior.betas)
      expect_identical(result.bvs$prior.models, result.btest$prior.models)

      expect_equal(unname(result.bvs$modelsprob[c(4,2,3,1), "Post"]),
                   unname(result.btest$PostProbi))
      expect_equal(unname(result.bvs$priorprobs[c(1,3)]),
                   unname(result.btest$priorprobs[c(1,4)]))
    }
  }
})

#Check that errors are given correctly and other printed messages
test_that("btest.glm reports errors and warnings", {
  data <- glm.test.data()

  expect.test.error(
    btest.glm(data = data, family = binomial(), models = c("M0", "M1")),
    "Argument models should be a list."
  )

  expect.test.error(
    btest.glm(data = data, family = binomial(),
              models = list(y ~ 1, y ~ x1), null.model = 1),
    "Please provide a name for the competing models."
  )

  expect.test.error(
    btest.glm(data = data, family = binomial(), models = list(y ~ 1)),
    "At least two models have to be provided."
  )

  models <- list(null = y ~ 1, x1 = y ~ x1)

  expect.test.error(
    btest.glm(data = data, family = binomial(), models = models, null.model = "unknown"),
    "null model provided is not in the list"
  )

  expect.test.error(
    btest.glm(data = data, family = binomial(), models = models, null.model = "x1"),
    "The given null model does not coincide with the one"
  )

  models2 <- list(x1 = y ~ x1, x2 = y ~ x2, full = y ~ x1 + x2)

  expect.test.error(
    btest.glm(data = data, models = models2, family = binomial(),  null.model = "x2"),
    "The simplest \\(null\\) model may not be nested in all the others\\."
  )

  expect.test.error(
    btest.glm(data = data, family = binomial(), models = models, prior.models = "SB"),
    "Only priors 'ScottBerger', 'Constant' and 'User' supported."
  )

  expect.test.error(
    btest.glm(data = data, models = models, family = binomial(),
             BF.method = "gprior", prior.betas = "gprior"),
    "Prior.betas must be one of 'gZellner', 'Robust',"
  )

  expect.test.error(
    btest.glm(data = data, models = list(null = y ~ -1, full = y ~ x1 + x2),
             null.model = "null"),
    "Null model shall contain the intercept term"
  )

  expect.test.error(
    btest.glm(data = data, models = list(null = y ~ -1, full = y ~ x1 + x2)),
    "Null model shall contain the intercept term"
  )

  output <- capture.output({
    result <- btest.glm(
      data = data,
      family = binomial(),
      models = models2
    )

    print(result)
  })

  expect_equal(output[14],
               "log Bayes factors (expressed in relation to x2)")

})







# for (i in seq_along(glm_families())) {
#   fam_string <- names(glm_families())[i]
#   family_obj <- glm_families()[[i]]
#
#   for (BF.method in c("BIC", "TBF")) {
#     test_that(sprintf(
#       "btest.glm computes posterior probabilities for %s with %s",
#       encode_family(family_obj), BF.method
#     ), {
#       data <- glm_family_data(fam_string)
#       models <- list(
#         null = y ~ 1,
#         x1 = y ~ x1,
#         full = y ~ x1 + x2
#       )
#
#       capture.output({
#         result <- suppressWarnings(btest.glm(
#           data = data,
#           models = models,
#           family = family_obj,
#           BF.method = BF.method
#         ))
#       })
#
#       expect_s3_class(result, "meantobtest")
#       expect_identical(result$nullmodel, "null")
#       expect_identical(names(result$lBFi0), c("null.to.null", "x1.to.null", "full.to.null"))
#       expect_equal(result$lBFi0[["null.to.null"]], 0, tolerance = 1e-6)
#       expect_true(all(is.finite(result$lBFi0)))
#       expect_equal(sum(result$PostProbi), 1, tolerance = 1e-6)
#       expect_length(result$modelspool, 3)
#       expect_true(all(vapply(result$modelspool, inherits, logical(1), "glm")))
#       expect_identical(result$family$family, family_obj$family)
#       expect_identical(result$BF.method, BF.method)
#       if (BF.method == "TBF") {
#         expect_identical(result$prior.betas, "gZellner")
#       }
#     })
#   }
# }
#
# for (fam_string in names(glm_bas_families())) {
#   family_obj <- glm_bas_family(fam_string)
#
#   test_that(sprintf(
#     "btest.glm supports the gprior method for %s",
#     encode_family(family_obj)
#   ), {
#     data <- glm_family_data(fam_string)
#     models <- list(
#       null = y ~ 1,
#       x1 = y ~ x1,
#       full = y ~ x1 + x2
#     )
#
#     capture.output({
#       result <- suppressWarnings(btest.glm(
#         data = data,
#         models = models,
#         family = family_obj,
#         BF.method = "gprior"
#       ))
#     })
#
#     expect_s3_class(result, "meantobtest")
#     expect_identical(result$nullmodel, "null")
#     expect_equal(result$lBFi0[["null.to.null"]], 0, tolerance = 1e-6)
#     expect_true(all(is.finite(result$lBFi0)))
#     expect_equal(sum(result$PostProbi), 1, tolerance = 1e-6)
#     expect_identical(result$family$family, family_obj$family)
#     expect_identical(result$BF.method, "gprior")
#     expect_identical(result$prior.betas, "gZellner")
#   })
# }
#
# for (i in seq_along(glm_families())) {
#   fam_string <- names(glm_families())[i]
#   family_obj <- glm_families()[[i]]
#
#   test_that(sprintf(
#     "btest.glm accepts external imputations for %s with BIC",
#     encode_family(family_obj)
#   ), {
#     data <- glm_family_data(fam_string, missing = TRUE)
#     imputations <- glm_test_imputations(data)
#     models <- list(
#       null = y ~ 1,
#       x1 = y ~ x1,
#       full = y ~ x1 + x2
#     )
#
#     capture.output({
#       result <- suppressWarnings(btest.glm(
#         data = data,
#         models = models,
#         family = family_obj,
#         BF.method = "BIC",
#         imp.datasets = imputations
#       ))
#     })
#
#     expect_s3_class(result, "meantobtest")
#     expect_true("compress.imp.array" %in% names(result))
#     expect_equal(result$imp.info$n.imp, 1)
#     expect_equal(sum(result$PostProbi), 1, tolerance = 1e-6)
#     expect_length(result$modelspool, 3)
#     expect_identical(result$family$family, family_obj$family)
#   })
# }
#

