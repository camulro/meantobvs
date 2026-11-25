#Make tests reproducible and non-dependent on external data

#Build small test data for linear models. If missing == TRUE, 2 NA are forced on x1
lm.test.data <- function(missing = FALSE) {
  n <- 24
  angle <- seq(0, 2 * pi, length.out = n)

  #build variables
  x1 <- seq(-1, 1, length.out = n) #24 values equally spaced from -1 to 1
  x2 <- cos(angle)
  y <- 1 + 0.8 * x1 - 0.4 * x2 + sin(angle) / 10
  #sin(angle) acts as a litle disturbance (error term determined)

  data <- data.frame(y = y, x1 = x1, x2 = x2)
  #If asked, force two missing values on x1
  if (missing) data$x1[c(4, 11)] <- NA_real_

  data
}

#Data with two numerical covariates and a three level factor
lm.testfactor.data <- function(missing = FALSE) {
  n <- 24
  angle <- seq(0, 2 * pi, length.out = n)

  #build variables
  x1 <- seq(-1, 1, length.out = n) #24 values equally spaced from -1 to 1
  x2 <- cos(angle)
  f1 <- factor(rep(c("L1", "L2", "L3"), n/3))
  f2 <- factor(rep(c("L1", "L2"), n/2), ordered = TRUE) #ordered factor

  data <- data.frame(x1 = x1, x2 = x2, f1 = f1, f2 = f2)
  data$y <- 1 + 0.6 * (f1 == "L3") + 0.8 * x1 - 0.4 * x2 + sin(angle) / 10

  if (missing) data$f1[c(4, 11)] <- NA_real_

  data
}

#Imputes the mean for variables with missing data and returns an array with it
lm.test.imputations <- function(data = lm.test.data(TRUE)) {
  imputed <- data
  #Impute observed mean for each variable with missing data
  for (j in 1:ncol(data)) {
    imputed[[j]][is.na(imputed[[j]])] <- mean(data[[j]], na.rm = TRUE)
  }

  #build an array for given imputations
  result <- array(rep(as.matrix(imputed), 1),
                  dim = c(nrow(imputed), ncol(imputed), 1),
                  dimnames = list(NULL, names(imputed), "1"))

  result
}

#Build design matrix for a linear model with intercept term
lm.test.complete.design <- function(data = lm.test.data()) {
  #The response variable has to be named "y"
  model.matrix(y ~ ., data)
}

#Get options available for prior.betas in lm
lm.prior.betas <- function(method = "gprior") {
  switch (method,
          gprior = c("gZellner",
                     "Robust",
                     "Liangetal",
                     "ZellnerSiow",
                     "FLS",
                     "intrinsic.MGC",
                     "IHG"),
          TBF = c("gZellner",
                  "Liangetal",
                  "ZellnerSiow",
                  "FLS"),
          BIC = NULL)
}

#Test that a function contains an error and it is the one expected
expect.test.error <- function(expr, pattern) {

  error <- tryCatch(force(expr), #evaluates expr and catches the error message
                    error = identity)

  #test fails if an error is not obtained
  testthat::expect_s3_class(error, "error")
  #checks if the message obtained is the expected one
  testthat::expect_match(conditionMessage(error), pattern)
}

#Build small test data for binomial models. If missing == TRUE, 2 NA are forced on x1
glm.test.data <- function(missing = FALSE, grouped = FALSE) {
  n <- 30
  angle <- seq(0, 2 * pi, length.out = n)
  x1 <- seq(-1.5, 1.5, length.out = n)
  x2 <- sin(angle)
  # y <- as.integer(seq_len(n) %% 4 %in% c(1, 2))
  eta <- 0.5 + 0.8 * x1 - 0.4 * x2

  if (!grouped) {
    y <- withr::with_seed(123,
                          rbinom(n, size = 1, prob = plogis(eta)))
    data <- data.frame(y = y, x1 = x1, x2 = x2)
  } else {
    y <- withr::with_seed(123, rbinom(n, size = 3, prob = plogis(eta)))

    control <- 3 - y
    data <- data.frame(y = y, control = control, x1 = x1, x2 = x2)
  }

  if (missing) data$x1[c(5, 16)] <- NA_real_
  data
}

#Imputes the mean for variables with missing data and returns an array with it
glm.test.imputations <- function(data) {
  imputed <- data
  imputed$x1[is.na(imputed$x1)] <- mean(data$x1, na.rm = TRUE)

  array(
    rep(as.matrix(imputed), 1),
    dim = c(nrow(imputed), ncol(imputed), 1),
    dimnames = list(NULL, names(imputed), "1")
  )
}

#available glm families
glm.families <- function() {
  list(
    "binomial" = binomial(),
    "gaussian" = gaussian(),
    "Gamma" = Gamma(link = "log"),
    "inverse.gaussian" = inverse.gaussian(),
    "poisson" = poisson(),
    # "quasi" = quasi(),
    "quasibinomial" = quasibinomial(),
    "quasipoisson" = quasipoisson()
  )
}

#available glm families that can be used through BAS package
glm.bas.families <- function() {
  list(
    "binomial" = binomial(link = "logit"),
    "Gamma" = Gamma(link = "log"),
    "poisson" = poisson(link = "log")
  )
}

glm.family.data <- function(family, missing = FALSE) {
  n <- 30
  angle <- seq(0, 2 * pi, length.out = n)

  #Build independent variables
  x1 <- seq(-1.5, 1.5, length.out = n)
  x2 <- sin(angle)

  eta <- 0.5 + 0.8 * x1 - 0.4 * x2

  #Select the kind of response needed for each family
  family.name <- switch(
    family,
    "binomial" = "binomial",
    "quasibinomial" = "binomial",
    "quasipoisson" = "poisson",
    "Gamma" = "inverse",
    "inverse.gaussian" = "inverse",
    "quasi" = "gaussian",
    family
  )

  #Build the response variable
  y <- switch(
    family.name,
    "binomial" = {
      withr::with_seed(123,
                       rbinom(n, size = 1, prob = plogis(eta)))
    },
      # as.integer(round(1.5 * x1 - 0.7 * x2) %in% c(1, 2)),
    "poisson" = {
      withr::with_seed(123,
                       rpois(n, lambda = exp(eta)))
    },
      # as.integer(round(2 + 0.8 * x1 - 0.2 * x2)),
    "inverse" = 0.5 + 1 * x1^2 + 0.5 * abs(x2),
    "gaussian" = eta + 0.25 * sin(angle),

    stop("unknown family '", family, "'", call. = FALSE)
  )

  if (family.name == "inverse") {
    data <- data.frame(y = y, x1 = x1^2, x2 = abs(x2))
  } else data <- data.frame(y = y, x1 = x1, x2 = x2)

  if (missing) data$x1[c(5, 16)] <- NA_real_
  data
}

#Get options available for prior.betas in glm
glm.prior.betas <- function(method = "gprior") {
    switch (method,
            gprior = c("gZellner",
                       "Robust",
                       "Liangetal",
                       "ZellnerSiow",
                       "FLS",
                       "intrinsic.WNC"),
            TBF = c("gZellner",
                    "Liangetal",
                    "ZellnerSiow",
                    "FLS"),
            BIC = NULL)
}

#Get combinations for families and links on glm stats:
getfamily.links <- function (family = NULL) {

  family.links <- list(
    "binomial" = c("logit", "probit", "cloglog", "cauchit", "log",
                   "identity"),
    "gaussian" = c("inverse", "log", "identity"),
    "Gamma" = c("inverse", "log", "identity"),
    "inverse.gaussian" = c("inverse", "log", "identity", "1/mu^2"),
    "poisson" = c("log", "identity", "sqrt"),
    # "quasi" = quasi(), #depends on the link and variance chosen
    "quasibinomial" = c("logit", "probit", "cloglog", "cauchit", "log",
                        "identity"),
    "quasipoisson" = c("log", "identity", "sqrt")
  )

  if (!is.null(family)) {
    checklink <- family$link %in% family.links[[family$family]]
    return(checklink)
  }

  family.links

}

