#Tests for Bayes factors functions

## lm functions:

#Check BIC in lm
test_that("BF.BIC.lm matches the BIC expression", {
  #build test deterministic data
  data <- lm.test.data()
  #Build design matrix
  X <- lm.test.complete.design(data)[, c("(Intercept)", "x1")]
  #Fit models
  null.fit <- lm(y ~ 1, data = data)
  model.fit <- lm.fit(x = X, y = data$y)

  #Get arguments to compute BF approx
  SS0 <- crossprod(null.fit$residuals)
  SSE <- crossprod(model.fit$residuals)
  n <- nrow(data)
  k <- 1

  result <- meantobvs:::BF.BIC.lm(
    y = data$y,
    X = X,
    SS0 = SS0,
    n = n,
    k = k,
    p0 = 1
  )

  expect_true(is.finite(result))
  expect_true(is.numeric(result))
  expect_true(is.double(result))
  expect_false(is.na(result))

  expected <- n/2 * log(as.numeric(SS0/SSE)) - k/2 * log(n)

  expect_equal(as.numeric(result), expected)
})

#Check TBF in lm: just gZellner approximation for now. It should be improved
#to cover the other BF approximations
test_that("BF.TBF.lm - gZellner matches the TBF expression", {
  #build test deterministic data
  data <- lm.test.data()
  #Build design matrix
  X <- lm.test.complete.design(data)[, c("(Intercept)", "x1")]
  #Fit models
  null.fit <- lm(y ~ 1, data = data)
  model.fit <- lm.fit(x = X, y = data$y)

  #Get arguments to compute BF approx
  SS0 <- crossprod(null.fit$residuals)
  SSE <- crossprod(model.fit$residuals)

  #Get TBF function
  lTBF.method <- lTBF.function("gZellner", nrow(data), 2, devnull = 0)

  result <- meantobvs:::BF.TBF.lm(
    y = data$y,
    X = X,
    SS0 = SS0,
    lTBF.method = lTBF.method,
    n = nrow(data),
    k = 1,
    p0 = 1
  )

  expect_true(is.finite(result))
  expect_true(is.numeric(result))
  expect_true(is.double(result))
  expect_false(is.na(result))

  expected.dev <- nrow(data) * log(as.numeric(SSE / SS0))
  #Get deviance from TBF - gZellner formula
  g <- nrow(data)
  result.dev <- -(as.numeric(result) + 1/2 * log(g + 1)) / (g / (g + 1)) * 2

  expect_equal(result.dev, expected.dev)

  #Other TBF:
  pb <- lm.prior.betas("TBF")[-1]
  for (j in pb) {

    #Get TBF function
    lTBF.method <- lTBF.function(j, nrow(data), 2, devnull = 0)

    result <- meantobvs:::BF.TBF.lm(
      y = data$y,
      X = X,
      SS0 = SS0,
      lTBF.method = lTBF.method,
      n = nrow(data),
      k = 1,
      p0 = 1
    )

    expect_true(is.finite(result))
    expect_true(is.numeric(result))
    expect_true(is.double(result))
    expect_false(is.na(result))

  }

})

#Check gprior in lm
test_that("g-prior Bayes factors returns the expected expression in lm", {
  #without missing data to compare
  data <- lm.test.data()

  #Check all BF available
  pb <- lm.prior.betas()[-3] #Liangetal fails in BayesVarSel
  pb <- pb[-4] #FLS uses a different BF function

  #Build design matrix
  X.model <- lm.test.complete.design(data)[, c("(Intercept)", "x1")]
  X.fullmodel <- lm.test.complete.design(data)

  #just intercept fixed
  SS0 <- as.numeric(crossprod(lm(y ~ 1, data = data)$residuals))
  #x1 fixed
  SS1 <- as.numeric(crossprod(lm(y ~ x1, data = data)$residuals))

  #gprior except FLS:

  #Full enumeration of models to compare with BayesVarSel output
  models <- list(null = y ~ 1, x1 = y ~ x1)
  models.full <- list(x1 = y ~ x1, full = y ~ .)
  for (j in pb) {
    #check BF for null model
    gprior <- meantobvs:::BF.gprior.lm(
      y = data$y,
      X = X.model[,"(Intercept)", drop = FALSE],
      SS0 = SS0,
      prior.betas = prior.betas.gprior.lm(j),
      n = nrow(data),
      k = 0L,
      p0 = 1L
    )
    expect_equal(as.numeric(gprior), 0)

    #just intercept  fixed
    result.meantobvs <- meantobvs:::BF.gprior.lm(
      y = data$y,
      X = X.model,
      SS0 = SS0,
      prior.betas = prior.betas.gprior.lm(j),
      n = nrow(data),
      k = ncol(X.model) - 1L,
      p0 = 1L
    )

    result.other <- BayesVarSel::Btest(
      data = data,
      models = models,
      prior.betas = j
    )$BFi0[2]

    expect_true(is.finite(result.meantobvs))
    expect_true(is.numeric(result.meantobvs))
    expect_true(is.double(result.meantobvs))
    expect_false(is.na(result.meantobvs))

    expect_equal(result.meantobvs, unname(log(result.other)))

    #x1 fixed
    result.meantobvs <- meantobvs:::BF.gprior.lm(
      y = data$y,
      X = X.fullmodel,
      SS0 = SS1,
      prior.betas = prior.betas.gprior.lm(j),
      n = nrow(data),
      k = ncol(X.fullmodel) - 2L,
      p0 = 2L
    )

    result.other <- BayesVarSel::Btest(
      data = data,
      models = models.full,
      prior.betas = j
    )$BFi0[2]

    expect_true(is.finite(result.meantobvs))
    expect_true(is.numeric(result.meantobvs))
    expect_true(is.double(result.meantobvs))
    expect_false(is.na(result.meantobvs))

    expect_equal(result.meantobvs, unname(log(result.other)))
  }

  #FLS:
  #Check BF for null model
  fls <- meantobvs:::BF.FLS.lm(
    y = data$y,
    X = X.model[,"(Intercept)", drop = FALSE],
    SS0 = SS0,
    dmax = 2,
    n = nrow(data),
    k = 0L,
    p0 = 1L
  )
  expect_equal(as.numeric(fls), 0)

  #just intercept  fixed
  result.meantobvs <- meantobvs:::BF.FLS.lm(
    y = data$y,
    X = X.model,
    SS0 = SS0,
    dmax = 2,
    n = nrow(data),
    k = ncol(X.model) - 1L,
    p0 = 1L
  )

  result.other <- BayesVarSel::Btest(
    data = data,
    models = models,
    prior.betas = "FLS"
  )$BFi0[2]

  expect_true(is.finite(result.meantobvs))
  expect_true(is.numeric(result.meantobvs))
  expect_true(is.double(result.meantobvs))
  expect_false(is.na(result.meantobvs))

  expect_equal(result.meantobvs, unname(log(result.other)))

  #x1 fixed
  result.meantobvs <- meantobvs:::BF.FLS.lm(
    y = data$y,
    X = X.fullmodel,
    SS0 = SS1,
    dmax = 3,
    n = nrow(data),
    k = ncol(X.fullmodel) - 2L,
    p0 = 2L
  )

  result.other <- BayesVarSel::Btest(
    data = data,
    models = models.full,
    prior.betas = "FLS"
  )$BFi0[2]

  expect_true(is.finite(result.meantobvs))
  expect_true(is.numeric(result.meantobvs))
  expect_true(is.double(result.meantobvs))
  expect_false(is.na(result.meantobvs))

  expect_equal(result.meantobvs, unname(log(result.other)))
})

#Check prior.betas user choice in lm
test_that("checkforprior.betas.lm selects the requested methods", {
  data <- lm.test.data()
  #Build design matrix
  X <- lm.test.complete.design(data)[, c("(Intercept)", "x1")]
  #just intercept fixed
  SS0 <- as.numeric(crossprod(lm(y ~ 1, data = data)$residuals))
  n <- nrow(data)
  k <- 1L

  #BIC
  bic <- meantobvs:::checkforprior.betas.lm(
    "BIC", NULL, n, p = 1, p0 = 1, y = data$y, SS0 = SS0
  )
  expect_true(is.function(bic))

  direct.bic <- meantobvs:::BF.BIC.lm(
    y = data$y, X = X, SS0 = SS0, n = n, k = k, p0 = 1
  )
  expect_equal(bic(k = k, X = X), direct.bic)

  #TBF
  #Get TBF function
  lTBF.method <- lTBF.function("gZellner", n, p = 2, devnull = 0)

  tbf <- meantobvs:::checkforprior.betas.lm(
    "TBF", NULL, n, p = 1, p0 = 1, y = data$y, SS0 = SS0
  )
  expect_true(is.function(tbf))

  direct.tbf <- meantobvs:::BF.TBF.lm(
    y = data$y, X = X, SS0 = SS0, lTBF.method = lTBF.method, n = n, k = k, p0 = 1L
  )
  expect_equal(tbf(k = k, X = X), direct.tbf)

  #gprior
  gprior <- meantobvs:::checkforprior.betas.lm(
    "gprior", "gZellner", n, p = 1, p0 = 1L, y = data$y, SS0 = SS0
  )
  expect_length(gprior(k = k, X = X), 1)
  expect_true(is.finite(gprior(k = k, X = X)))
  expect_true(is.function(gprior))

  direct.gprior <- meantobvs:::BF.gprior.lm(
    y = data$y, X = X, SS0 = SS0, prior.betas = "gBF", n = n, k = k, p0 = 1L
  )
  expect_equal(gprior(k = k, X = X), direct.gprior)

  #expected errors
  expect.test.error(
    meantobvs:::checkforprior.betas.lm(
      "gZellner", NULL, n = 10, p = 1, p0 = 1, y = 1:10, SS0 = 1
    ),
    "Only BF approximations"
  )
  expect.test.error(
    meantobvs:::checkforprior.betas.lm(
      "TBF", "BIC", n = 10, p = 1, p0 = 1, y = 1:10, SS0 = 1
    ),
    "Prior.betas must be one of"
  )
})

## glm functions:

#Check BIC in glm
test_that("BF.BIC.glm.fit matches the BIC expression", {

  #go through all families available
  families <- glm.families() #with default links
  bas.families <- sapply(1:3,
                         function (x) family.link(glm.bas.families()[[x]]))

  for (i in names(families)) {
    family.obj <- families[[i]]
    #build test deterministic data for family i
    data <- glm.family.data(i)

    #Build design matrix
    X <- lm.test.complete.design(data)[, c("(Intercept)", "x1")]
    X2 <- lm.test.complete.design(data)[, c("(Intercept)", "x1", "x2")]

    #Fit models
    null.fit <- glm(y ~ 1, data = data, family = family.obj)
    model.fit <- glm.fit(x = X, y = data$y, family = family.obj)
    model.fit2 <- glm.fit(x = X2, y = data$y, family = family.obj)

    #Get arguments to compute BF approx
    n <- nrow(data)
    k <- 1L #number of competing variables

    result <- meantobvs:::BF.BIC.glm.fit(
      y = data$y,
      X = X,
      family = family.obj,
      devnull = null.fit$deviance,
      n = n, k = k
    )

    expect_true(is.finite(result))
    expect_true(is.numeric(result))
    expect_true(is.double(result))
    expect_false(is.na(result))

    #Check formula
    expected <- (null.fit$deviance - model.fit$deviance - k * log(n))/2
    expect_equal(as.numeric(result), expected)

    #with fixed x1
    result2 <- meantobvs:::BF.BIC.glm.fit(
      y = data$y,
      X = X2,
      family = family.obj,
      devnull = model.fit$deviance,
      n = n, k = k
    )

    expect_true(is.finite(result2))
    expect_true(is.numeric(result2))
    expect_true(is.double(result2))
    expect_false(is.na(result2))

    #Check formula
    expected <- (model.fit$deviance - model.fit2$deviance - k * log(n))/2
    expect_equal(as.numeric(result2), expected)

    #Compare with other functions, when we can
    if (family.link(family.obj) %in% bas.families) {

      family <- family.obj
      other.result <- BAS::bas.glm(data = data,
                                  formula = y ~ . ,
                                  family = family,
                                  betaprior = BAS::bic.prior(),
                                  method = "deterministic",
                                  n.models = 4)

      m1 <- sapply(1:4, function(x) identical(other.result$which[[x]], c(0L,1L)))
      m2 <- sapply(1:4, function(x) identical(other.result$which[[x]], 0L))

      expect_equal(result,
                   other.result$logmarg[which(m1)] - other.result$logmarg[which(m2)])
      #with x1 fixed
      m1 <- sapply(1:4, function(x) identical(other.result$which[[x]], c(0L,1L,2L)))
      m2 <- sapply(1:4, function(x) identical(other.result$which[[x]], c(0L,1L)))

      expect_equal(result2,
                   other.result$logmarg[which(m1)] - other.result$logmarg[which(m2)])
    }

  }
})

#Check TBF - gZellner in glm
test_that("BF.TBF.glm.fit matches the TBF expression", {

  #go through all families available
  families <- glm.families() #with default links
  bas.families <- sapply(1:3,
                         function (x) family.link(glm.bas.families()[[x]]))

  for (i in names(families)) {
    family.obj <- families[[i]]
    #build test deterministic data for family i
    data <- glm.family.data(i)

    #Build design matrix
    X <- lm.test.complete.design(data)[, c("(Intercept)", "x1")]
    X2 <- lm.test.complete.design(data)[, c("(Intercept)", "x1", "x2")]

    #Fit models
    null.fit <- glm(y ~ 1, data = data, family = family.obj)
    model.fit <- glm.fit(x = X, y = data$y, family = family.obj)
    model.fit2 <- glm.fit(x = X2, y = data$y, family = family.obj)

    #Get arguments to compute BF approx
    n <- nrow(data)
    k <- 1L #number of competing variables
    devnull <- null.fit$deviance

    #Get TBF function
    lTBF.method <- lTBF.function("gZellner", nrow(data), 2, devnull = devnull)

    result <- meantobvs:::BF.TBF.glm.fit(
      y = data$y,
      X = X,
      family = family.obj,
      n = n, k = k,
      lTBF.method = lTBF.method
    )

    expect_true(is.finite(result))
    expect_true(is.numeric(result))
    expect_true(is.double(result))
    expect_false(is.na(result))

    g <- n
    expected <- -k/2 * log(g + 1) + g/(g+1) * (devnull - model.fit$deviance)/2
    expect_equal(as.numeric(result), expected)

    #with fixed x1
    lTBF.method <- lTBF.function("gZellner", nrow(data), 2, devnull = model.fit$deviance)

    result2 <- meantobvs:::BF.TBF.glm.fit(
      y = data$y,
      X = X2,
      family = family.obj,
      n = n, k = k,
      lTBF.method = lTBF.method
    )

    expect_true(is.finite(result2))
    expect_true(is.numeric(result2))
    expect_true(is.double(result2))
    expect_false(is.na(result2))

    g <- n
    expected <- -k/2 * log(g + 1) + g/(g+1) * (model.fit$deviance - model.fit2$deviance)/2
    expect_equal(as.numeric(result2), expected)

    #Compare with other functions, when we can
    if (family.link(family.obj) %in% bas.families) {
      family <- family.obj
      other.result <- BAS::bas.glm(data = data,
                                  formula = y ~ . ,
                                  family = family,
                                  betaprior = BAS::testBF.prior(g = n),
                                  method = "deterministic",
                                  n.models = 4)

      m1 <- sapply(1:4, function(x) identical(other.result$which[[x]], c(0L,1L)))
      m2 <- sapply(1:4, function(x) identical(other.result$which[[x]], 0L))

      expect_equal(result,
                   other.result$logmarg[which(m1)] - other.result$logmarg[which(m2)])
      #with x1 fixed
      m1 <- sapply(1:4, function(x) identical(other.result$which[[x]], c(0L,1L,2L)))
      m2 <- sapply(1:4, function(x) identical(other.result$which[[x]], c(0L,1L)))

      expect_equal(result2,
                   other.result$logmarg[which(m1)] - other.result$logmarg[which(m2)])
    }

    #Other TBF:
    pb <- glm.prior.betas("TBF")[-1]
    for (j in pb) {

      #Get TBF function
      lTBF.method <- lTBF.function(j, nrow(data), 2, devnull = devnull)

      result <- meantobvs:::BF.TBF.glm.fit(
        y = data$y,
        X = X,
        family = family.obj,
        n = n, k = k,
        lTBF.method = lTBF.method
      )

      expect_true(is.finite(result))
      expect_true(is.numeric(result))
      expect_true(is.double(result))
      expect_false(is.na(result))

    }

  }
})

#Check gprior in glm
test_that("g-prior Bayes factors returns the expected expression in glm", {

  #Check all BF available
  pb <- glm.prior.betas()[-4] #Zellner-Siow does not agree?
  #go through all families available
  families <- glm.families() #Gamma does not work ok? #with default links
  bas.families <- sapply(1:3,
                         function (x) family.link(glm.bas.families()[[x]]))

  for (i in names(families)) {

    family.obj <- families[[i]]
    #build test deterministic data for family i
    data <- glm.family.data(i)

    #Build design matrix
    X <- lm.test.complete.design(data)[, c("(Intercept)", "x1")]
    X2 <- lm.test.complete.design(data)[, c("(Intercept)", "x1", "x2")]

    #Fit models
    null.fit <- glm(y ~ 1, data = data, family = family.obj, y = T, x = T)
    model.fit <- glm.fit(x = X, y = data$y, family = family.obj)
    model.fit2 <- glm.fit(x = X2, y = data$y, family = family.obj)

    #Get arguments to compute BF approx
    n <- nrow(data)
    k <- 1L #number of competing variables
    devnull <- null.fit$deviance

    for (j in pb) {

      args <- prior.betas.gprior.glm.fit(j, n = n, p = 2)

      #check BF for null model
      if (args$type == "fixed") {
        result <- meantobvs:::BF.gprior.glm.fit(
          y = data$y,
          X = X[,"(Intercept)", drop = FALSE],
          glm0 = null.fit,
          g = args$g, n = n, k = 0L
        )
      } else {
        result <- meantobvs:::BF.hypergprior.glm.fit(
          y = data$y,
          X = X[,"(Intercept)", drop = FALSE],
          glm0 = null.fit,
          prior.betas.args = args[-1],
          n = n, k = 0L
        )
      }

      expect_equal(as.numeric(result), 0)

      #just intercept fixed
      if (args$type == "fixed") {
        result <- meantobvs:::BF.gprior.glm.fit(
          y = data$y,
          X = X,
          glm0 = null.fit,
          g = args$g, n = n, k = k
        )
      } else {
        result <- meantobvs:::BF.hypergprior.glm.fit(
          y = data$y,
          X = X,
          glm0 = null.fit,
          prior.betas.args = args[-1],
          n = n, k = k
        )
      }

      expect_true(is.finite(result))
      expect_true(is.numeric(result))
      expect_true(is.double(result))
      expect_false(is.na(result))

      # args <- prior.betas.gprior.glm.fit(j, n = n, p = 1) # error?

      #with fixed x1
      if (args$type == "fixed") {
        result2.1 <- meantobvs:::BF.gprior.glm.fit(
          y = data$y,
          X = X2,
          glm0 = null.fit,
          g = args$g, n = n, k = k + 1L
        )

        result2.2 <- meantobvs:::BF.gprior.glm.fit(
          y = data$y,
          X = X,
          glm0 = null.fit,
          g = args$g, n = n, k = k
        )

        result2 <- result2.1 - result2.2

      } else {
        result2.1 <- meantobvs:::BF.hypergprior.glm.fit(
          y = data$y,
          X = X2,
          glm0 = null.fit,
          prior.betas.args = args[-1],
          n = n, k = k + 1L
        )

        result2.2 <- meantobvs:::BF.hypergprior.glm.fit(
          y = data$y,
          X = X,
          glm0 = null.fit,
          prior.betas.args = args[-1],
          n = n, k = k
        )

        result2 <- result2.1 - result2.2
      }

      expect_true(is.finite(result2))
      expect_true(is.numeric(result2))
      expect_true(is.double(result2))
      expect_false(is.na(result2))

      #changing weights and offset:
      family <- family.obj
      null.fit2 <- glm(y ~ 1, data = data, family = family,
                       weights = c(rep(1, n/2), rep(2, n/2)),
                       # offset = rep(1, n),
                       y = T, x = T)

      #just intercept fixed
      if (args$type == "fixed") {
        result3 <- meantobvs:::BF.gprior.glm.fit(
          y = null.fit2$y,
          X = X,
          glm0 = null.fit2,
          g = args$g, n = n, k = k,
          weights = null.fit2$prior.weights
        )
      } else {
        result3 <- meantobvs:::BF.hypergprior.glm.fit(
          y = null.fit2$y,
          X = X,
          glm0 = null.fit2,
          prior.betas.args = args[-1],
          n = n, k = k,
          weights = null.fit2$prior.weights
        )
      }

      expect_true(is.finite(result3))
      expect_true(is.numeric(result3))
      expect_true(is.double(result3))
      expect_false(is.na(result3))

    #Compare with other functions, when we can
    if (family.link(family.obj) %in% bas.families) {

      #direct comp
      m0 <- logmargnull.glm(null.fit)
      prior.betas <- prior.betas.gprior.glm(j, n = n, p = 2)

      other.result <- meantobvs:::BF.gprior.glm(
        y = as.numeric(data$y),
        X = X,
        family = family.obj,
        prior.betas = prior.betas,
        logmargnull = m0,
        k = k
      )

      expect_equal(result, other.result, tolerance = 10^(-5))

      #other packages
      family <- family.obj
      other.result <- BAS::bas.glm(data = data,
                                  formula = y ~ . ,
                                  family = family,
                                  betaprior = prior.betas,
                                  method = "deterministic",
                                  n.models = 4)

      m1 <- sapply(1:4, function(x) identical(other.result$which[[x]], c(0L,1L)))
      m2 <- sapply(1:4, function(x) identical(other.result$which[[x]], 0L))

      expect_equal(result,
                   other.result$logmarg[which(m1)] - other.result$logmarg[which(m2)],
                   tolerance = 10^(-5))
      #with x1 fixed
      m1 <- sapply(1:4, function(x) identical(other.result$which[[x]], c(0L,1L,2L)))
      m2 <- sapply(1:4, function(x) identical(other.result$which[[x]], c(0L,1L)))

      expect_equal(result2,
                   other.result$logmarg[which(m1)] - other.result$logmarg[which(m2)],
                   tolerance = 10^(-5))

      #changing weights and offset
      weights <- null.fit2$prior.weights
      # offset <- null.fit2$offset
      family <- family.obj
      other.result <- BAS::bas.glm(data = data,
                                   formula = y ~ . ,
                                   family = family,
                                   weights = weights,
                                   # offset = offset,
                                   betaprior = prior.betas,
                                   method = "deterministic",
                                   n.models = 4)

      m1 <- sapply(1:4, function(x) identical(other.result$which[[x]], c(0L,1L)))
      m2 <- sapply(1:4, function(x) identical(other.result$which[[x]], 0L))

      expect_equal(result3,
                   other.result$logmarg[which(m1)] - other.result$logmarg[which(m2)],
                   tolerance = 10^(-4)) #ojo!  `actual`: 1.293543 `expected`: 1.293560

      }

    }

  }

  #Grouped binomial:
  family.obj <- binomial(link = "logit")
  #build test deterministic data for family i
  data <- glm.test.data(grouped = TRUE)

  #Build design matrix
  X <- lm.test.complete.design(data)[, c("(Intercept)", "x1")]
  X2 <- lm.test.complete.design(data)[, c("(Intercept)", "x1", "x2")]

  #Fit models
  null.fit <- glm(cbind(y, control) ~ 1, data = data, family = family.obj, y = T, x = T)
  model.fit <- glm(cbind(y, control) ~ x1, data = data, family = family.obj, y = T, x = T)
  model.fit2 <- glm(cbind(y, control) ~ x1 + x2, data = data, family = family.obj, y = T, x = T)

  #Get arguments to compute BF approx
  n <- nrow(data)
  k <- 1L #number of competing variables
  devnull <- null.fit$deviance

  for (j in pb) {

    args <- prior.betas.gprior.glm.fit(j, n = n, p = 2)

    #check BF for null model
    if (args$type == "fixed") {
      result <- meantobvs:::BF.gprior.glm.fit(
        y = null.fit$y,
        X = X[,"(Intercept)", drop = FALSE],
        glm0 = null.fit,
        g = args$g, n = n, k = 0L,
        weights = null.fit$prior.weights
      )
    } else {
      result <- meantobvs:::BF.hypergprior.glm.fit(
        y = null.fit$y,
        X = X[,"(Intercept)", drop = FALSE],
        glm0 = null.fit,
        prior.betas.args = args,
        n = n, k = 0L,
        weights = null.fit$prior.weights
      )
    }

    expect_equal(as.numeric(result), 0)

    #just intercept fixed
    if (args$type == "fixed") {
      result <- meantobvs:::BF.gprior.glm.fit(
        y = null.fit$y,
        X = X,
        glm0 = null.fit,
        g = args$g, n = n, k = k,
        weights = null.fit$prior.weights
      )
    } else {
      result <- meantobvs:::BF.hypergprior.glm.fit(
        y = null.fit$y,
        X = X,
        glm0 = null.fit,
        prior.betas.args = args,
        n = n, k = k,
        weights = null.fit$prior.weights
      )
    }

    expect_true(is.finite(result))
    expect_true(is.numeric(result))
    expect_true(is.double(result))
    expect_false(is.na(result))

    #with fixed x1
    if (args$type == "fixed") {
      result2.1 <- meantobvs:::BF.gprior.glm.fit(
        y = null.fit$y,
        X = X2,
        glm0 = null.fit,
        g = args$g, n = n, k = k + 1L,
        weights = null.fit$prior.weights
      )

      result2.2 <- meantobvs:::BF.gprior.glm.fit(
        y = null.fit$y,
        X = X,
        glm0 = null.fit,
        g = args$g, n = n, k = k,
        weights = null.fit$prior.weights
      )

      result2 <- result2.1 - result2.2

    } else {
      result2.1 <- meantobvs:::BF.hypergprior.glm.fit(
        y = null.fit$y,
        X = X2,
        glm0 = null.fit,
        prior.betas.args = args,
        n = n, k = k + 1L,
        weights = null.fit$prior.weights
      )

      result2.2 <- meantobvs:::BF.hypergprior.glm.fit(
        y = null.fit$y,
        X = X,
        glm0 = null.fit,
        prior.betas.args = args,
        n = n, k = k,
        weights = null.fit$prior.weights
      )

      result2 <- result2.1 - result2.2
    }

    expect_true(is.finite(result2))
    expect_true(is.numeric(result2))
    expect_true(is.double(result2))
    expect_false(is.na(result2))

    #changing weights and offset:
    null.fit2 <- glm(cbind(y, control) ~ 1, data = data, family = family.obj,
                     weights = c(rep(1, n/2), rep(2, n/2)),
                     # offset = rep(1, n),
                     y = T, x = T)

    #just intercept fixed
    if (args$type == "fixed") {
      result3 <- meantobvs:::BF.gprior.glm.fit(
        y = null.fit2$y,
        X = X,
        glm0 = null.fit2,
        g = args$g, n = n, k = k,
        weights = null.fit2$prior.weights
      )
    } else {
      result3 <- meantobvs:::BF.hypergprior.glm.fit(
        y = null.fit2$y,
        X = X,
        glm0 = null.fit2,
        prior.betas.args = args,
        n = n, k = k,
        weights = null.fit2$prior.weights #que pesos se ponen?? cbind te genera unos
      )
    }

    expect_true(is.finite(result3))
    expect_true(is.numeric(result3))
    expect_true(is.double(result3))
    expect_false(is.na(result3))

    #Compare with other functions:
    #direct comp
    m0 <- logmargnull.glm(null.fit)
    prior.betas <- prior.betas.gprior.glm(j, n = n, p = 2)

    # other.result <- meantobvs:::BF.gprior.glm(
    #   y = null.fit$y,
    #   X = X,
    #   family = family.obj,
    #   prior.betas = prior.betas,
    #   logmargnull = m0,
    #   k = k,
    #   weights = null.fit2$prior.weights
    # )
    #
    # expect_equal(result, other.result, tolerance = 10^(-5))

    #other packages
    family <- family.obj
    other.result <- BAS::bas.glm(data = data,
                                 formula = cbind(y, control) ~ . ,
                                 family = family,
                                 betaprior = prior.betas,
                                 method = "deterministic",
                                 n.models = 4)

    m1 <- sapply(1:4, function(x) identical(other.result$which[[x]], c(0L,1L)))
    m2 <- sapply(1:4, function(x) identical(other.result$which[[x]], 0L))

    expect_equal(result,
                 other.result$logmarg[which(m1)] - other.result$logmarg[which(m2)],
                 tolerance = 10^(-2)) ### low tolerance
    #with x1 fixed
    m1 <- sapply(1:4, function(x) identical(other.result$which[[x]], c(0L,1L,2L)))
    m2 <- sapply(1:4, function(x) identical(other.result$which[[x]], c(0L,1L)))

    expect_equal(result2,
                 other.result$logmarg[which(m1)] - other.result$logmarg[which(m2)],
                 tolerance = 10^(-1))

    #changing weights and offset
    weights <- c(rep(1, n/2), rep(2, n/2))
    other.result <- BAS::bas.glm(data = data,
                                 formula = cbind(y, control) ~ .,
                                 family = family,
                                 weights = weights,
                                 # offset = null.fit2$offset,
                                 betaprior = prior.betas,
                                 method = "deterministic",
                                 n.models = 4)

    m1 <- sapply(1:4, function(x) identical(other.result$which[[x]], c(0L,1L)))
    m2 <- sapply(1:4, function(x) identical(other.result$which[[x]], 0L))

    expect_equal(result3,
                 other.result$logmarg[which(m1)] - other.result$logmarg[which(m2)],
                 tolerance = 10^(-2))

  }

})

#Check prior.betas user choice in glm
test_that("checkforprior.betas.glm selects the requested methods", {
  data <- glm.test.data()
  family <- binomial(link = "logit")
  #Build design matrix
  X <- lm.test.complete.design(data)[, c("(Intercept)", "x1")]
  #Fit models
  null.fit <- glm(y ~ 1, data = data, family = family, x = T, y = T)
  devnull <- null.fit$deviance

  n <- nrow(data)
  k <- 1L
  null.fit$offset <- rep(0,n)

  #BIC
  bic <- meantobvs:::checkforprior.betas.glm(
    "BIC", NULL, n, p = 1, p0 = 1, y = data$y, null.fit, FALSE, 0L
  )
  expect_true(is.function(bic))

  direct.bic <- meantobvs:::BF.BIC.glm.fit(
    y = data$y, X = X, family = family, devnull = null.fit$deviance, n = n, k = k
  )
  expect_equal(bic(k = k, X = X, fitstart = NULL), direct.bic)

  #TBF
  #Get TBF function
  lTBF.method <- lTBF.function("gZellner", n, p = 2, devnull = devnull)

  tbf <- meantobvs:::checkforprior.betas.glm(
    "TBF", "gZellner", n, p = 1, p0 = 1, y = data$y, null.fit, FALSE, 0L
  )
  expect_true(is.function(tbf))

  direct.tbf <- meantobvs:::BF.TBF.glm.fit(
    y = data$y, X = X, family = family, n = n, k = k, lTBF.method = lTBF.method
  )
  expect_equal(tbf(k = k, X = X, fitstart = NULL), direct.tbf)

  #gprior
  gprior <- meantobvs:::checkforprior.betas.glm(
    "gprior", "gZellner", n, p = 1, p0 = 1, y = data$y, null.fit, FALSE, 0L
  )
  expect_length(gprior(k = k, X = X, fitstart = NULL), 1)
  expect_true(is.finite(gprior(k = k, X = X, fitstart = NULL)))
  expect_true(is.function(gprior))

  direct.gprior <- meantobvs:::BF.gprior.glm.fit(
    y = data$y, X = X, glm0 = null.fit, g = n, n = n, k = k
  )
  expect_equal(gprior(k = k, X = X, fitstart = NULL), direct.gprior)

  #expected errors
  expect.test.error(
    meantobvs:::checkforprior.betas.glm(
      "gZelner", NULL, n = 10, p = 1, p0 = 1, y = 1:10, null.fit, FALSE, 0L
    ),
    "Only BF approximations"
  )
  expect.test.error(
    meantobvs:::checkforprior.betas.glm(
      "TBF", "BIC", n = 10, p = 1, p0 = 1, y = 1:10, null.fit, FALSE, 0L
    ),
    "Prior.betas must be one of"
  )
})


