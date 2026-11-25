#Tests for the model priors performance

#Test the result of Constant prior
test_that("logConstant assigns 1/2^p to every model", {
  for (p in 1:5) {
    expect_equal(logConstant(p), -p * log(2))
    expect_equal(exp(logConstant(p)), 1 / 2^p)
  }
})

#Test the result of ScottBerger prior
test_that("logScottBerger is uniform over the model size", {
  #Build set of competing models
  models <- as.matrix(expand.grid(rep(list(0:1), 3)))

  #it matches its defining expression
  for (m in list(c(0, 0, 0), c(1, 0, 0), c(1, 1, 0), c(1, 1, 1))) {
    expect_equal(logScottBerger(p = 3, model = m),
                 -log(3 + 1) - lchoose(3, sum(m)))
  }

  #check that the total prior mass over the whole model space is 1
  mass <- sum(exp(apply(models, 1,
                        function(model) logScottBerger(p = 3, model = model))))
  expect_equal(mass, 1)

  #check that probability is uniform among the p+1 model sizes
  bysize <- vapply(0:3, function(k) {
    ind <- (rowSums(models) == k)
    sum(exp(apply(models[ind, , drop = FALSE], 1,
                  function(model) logScottBerger(p = 3, model = model))))
  }, numeric(1))
  expect_equal(bysize, rep(1 / 4, 4))
})

#Test how selects the prior
test_that("checkforprior.models selects the requested prior", {
  constant <- meantobvs:::checkforprior.models("Constant", NULL, p = 2)
  expect_equal(constant(c(1, 0)), logConstant(p = 2))

  scottberger <- meantobvs:::checkforprior.models("ScottBerger", NULL, p = 2)
  expect_equal(scottberger(c(1, 0)), logScottBerger(p = 2, model = c(1, 0)))

  user <- meantobvs:::checkforprior.models("User", c(0.2, 0.3, 0.5), p = 2)
  expect_equal(user(c(1, 0)), logUser(p = 2, model = c(1, 0),
                                      priorprobs = c(0.2, 0.3, 0.5)))
})

#Test the result of User prior
test_that("logUser follows the user given model size probabilities", {
  pp <- c(0.1, 0.2, 0.7)
  models <- as.matrix(expand.grid(rep(list(0:1), 2)))

  expect_equal(logUser(p = 2, model = c(1, 0), priorprobs = pp),
               log(pp[2]) - log(sum(pp)) - lchoose(2, 1))

  #the total mass over the model space is 1, whether or not priorprobs sum 1
  expect_equal(sum(exp(apply(models, 1,
                             function(model) logUser(p = 2, model = model,
                                                     priorprobs = pp)))), 1)
  expect_equal(sum(exp(apply(models, 1,
                             function(model) logUser(p = 2, model = model,
                                                     priorprobs = rep(1, 3))))), 1)
})

#Test how selects User prior
test_that("checkforprior.models validates the User prior", {
  p <- 2

  expect_test_error(
    meantobvs:::checkforprior.models("User", NULL, p = p),
    "no prior probabilities provided"
  )
  expect_test_error(
    meantobvs:::checkforprior.models("User", "x", p = p),
    "no numeric probabilities provided"
  )
  expect_test_error(
    meantobvs:::checkforprior.models("User", c(1, NA, 1), p = p),
    "some prior probabilities not provided"
  )
  expect_test_error(
    meantobvs:::checkforprior.models("User", c(1, 1), p = p),
    "length of prior"
  )
  expect_test_error(
    meantobvs:::checkforprior.models("User", c(1, -1, 1), p = p),
    "must be positive"
  )
  expect_test_error(
    meantobvs:::checkforprior.models("User", c(0, 0, 0), p = p),
    "must be positive"
  )
  expect_test_error(
    meantobvs:::checkforprior.models("uniform", NULL, p = p),
    "Only priors 'ScottBerger', 'Constant' and 'User' supported"
  )
})

#prior over dummy submodels:
test_that("model prior over the dummy submodel space", {
  #no active factor: 0
  expect_equal(logScottBerger.d(delta = c(0, 0), tau = c(FALSE, FALSE),
                                l = c(3, 4)), 0)
  expect_equal(logConstant.d(tau = c(FALSE, FALSE), l = c(3, 4)), 0)

  #single active factor:
  ##-lchoose(3,1) - log(3 - 1) = -log(3) - log(2)
  expect_equal(logScottBerger.d(delta = c(1, 0), tau = c(TRUE, FALSE),
                                l = c(3, 4)), -log(3) - log(2))
  ##-log(2^3 - 1 - 3) = -log(4)
  expect_equal(logConstant.d(tau = c(TRUE, FALSE), l = c(3, 4)), -log(4))

  #two active factors
  expect_equal(logScottBerger.d(delta = c(1, 2), tau = c(TRUE, TRUE),
                                l = c(3, 4)),
               -log(3) - lchoose(4, 2) - log(2) - log(3))
  expect_equal(logConstant.d(tau = c(TRUE, TRUE), l = c(3, 4)),
               -log(4) - log(2^4 - 1 - 4))
})

#Test how selects the prior over submodel space
test_that("checkforprior.models.dummies selects the requested prior", {

  l <- c(2, 3) #levels of 2 factors
  constant.d <- meantobvs:::checkforprior.models.dummies("Constant", l)
  expect_equal(constant.d(tau = c(T,F)), logConstant.d(c(T,F), l))
  expect_equal(constant.d(tau = c(F,T)), logConstant.d(c(F,T), l))
  expect_equal(constant.d(tau = c(T,T)), logConstant.d(c(T,T), l))

  scottberger.d <- meantobvs:::checkforprior.models.dummies("ScottBerger", l)
  expect_equal(scottberger.d(c(1,0), c(T,F)), logScottBerger.d(c(1,0), c(T,F), l))
  expect_equal(scottberger.d(c(0,2), c(F,T)), logScottBerger.d(c(0,2), c(F,T), l))
  expect_equal(scottberger.d(c(1,2), c(T,T)), logScottBerger.d(c(1,2), c(T,T), l))

})

