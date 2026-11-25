
#' Build design matrix (fixed and non-fixed), matrix to use -if needed- for
#' imputation step and factor objects -if needed- to use for computation.
#'
#' @keywords internal
buildmatrices <- function(formula, null.model, data, marginal.factors) {

  matrices <- list()

  ##Null model (fixed vars, including the intercept term)
  #for imputation:
  framenull <- model.frame(null.model, data = data, na.action = NULL)
  #for computation when no missing data present
  X0rdf <- X0 <- model.matrix(framenull, data)
  namesnull <- colnames(X0)

  matrices$q0 <- ncol(framenull) #covariate-factor level
  matrices$p0 <- ncol(X0) #covariate-dummy level

  ##Full model
  framefull <- model.frame(formula, data = data, na.action = NULL)

  matrices$X0 <- X0
  matrices$namesnull <- namesnull
  matrices$framenull <- framenull
  matrices$framefull <- framefull

  ##Factors' treatment: rank deficient decomposition
  if (marginal.factors) {

    X0rdf <- model.rankdef.matrix(framenull) #to remove vars from full
    X.full <- model.rankdef.matrix(framefull) #rank defficient full design

    namesx <- colnames(X.full)
    namesxnotnull <- setdiff(namesx, colnames(X0rdf)) #non-fixed var names

  } else {
    X.full <- model.matrix(formula, framefull)

    namesx <- colnames(X.full)
    namesxnotnull <- setdiff(namesx, namesnull)
  }
  #Full design matrix without fixed vars
  matrices$X.full <- X.full[, namesxnotnull, drop = FALSE]
  p <- length(namesxnotnull) #covariate-dummy level

  #Is there any variable to select from?
  if (p == 0) stop("The number of fixed variables is equal to the number of\n",
                   "regressors in the full model. No model selection can be done.\n")

  matrices$ordvars <- c(namesnull, namesxnotnull) #[X0, X.full]
  matrices$namesxnotnull <- namesxnotnull
  matrices$namesx <- namesx
  matrices$p <- p

  if (!marginal.factors) { #If no factors' treatment
    matrices$q <- p
    matrices$depvars <- namesxnotnull
    matrices$L <- 0L

    return(matrices)
  }

  ##Get covariates and factors to select from if needed
  #Info from framefull to get factors and dummy names:
  depvars <- setdiff(attr(terms(framefull), "term.labels"),
                     attr(terms(framenull), "term.labels"))

  matrices$depvars <- depvars

  factorsfull <- attr(terms(framefull), "factors")
  classesfull <- attr(terms(framefull), "dataClasses")

  #positions has number of rows equal to the number of regressors and p columns.
  #A 1 in a row denotes the position in X of a regressor (several positions for
  #the dummies of a factor).
  positions <- t(sapply(depvars,
    function(var) {

      ##OLD: (does not select interactions of factors)
      # if(is.factor(data[[var]])) {
      #   levs <- levels(data[[var]])
      #   ind <- which(namesxnotnull %in% paste0(var,levs)) #1 if the namelevel matches
      # } else ind <- which(namesxnotnull == var) #1 if the name matches

      if (classesfull[var] %in% c("factor", "ordered")) {

        #get vars corresponding to that factor (to include interactions)
        terms_var <- names(which(factorsfull[, var] > 0))
        levels_var <- lapply(terms_var,
          function(x) paste0(x, levels(data[[x]]))
        )
        #get all possible combinations for factor and levels
        gridlevs <- expand.grid(levels_var, stringsAsFactors = FALSE)

        #for interaction terms
        levelnames <- apply(gridlevs, 1, paste0, collapse = ":")

        #select dummies corresponding to each factor
        ind <- which(namesxnotnull %in% levelnames)

      } else ind <- which(namesxnotnull == var)

      #build position row for these ind
      pos <- numeric(p)
      pos[ind] <- 1
      pos
    }
  ))

  colnames(positions) <- namesxnotnull

  #Identify factors
  ifelse(dim(positions)[2] > 1,
         tmp <- colSums(positions %*% t(positions)),
         tmp <- positions)

  positionsx <- tmp == 1 #TRUE if numeric covariate
  L <- sum(!positionsx) #Number of factors to select from

  matrices$positions <- positions
  matrices$positionsx <- positionsx

  ##Factor information needed for model prior and posterior computation
  if (L > 0) {

    #1 if the dummy column variable belongs to the row factor
    positionsfac <- positions[!positionsx, , drop = FALSE]

    l <- tmp[tmp > 1] #number of levels for each factor

    #save position of the first dummy in each factor
    #used on pool estimation step to get representant over repeated models
    matrices$indf <- apply(positionsfac, 1,
      function(x) which(x == 1)[1]
    )

    #save hash for representant over saturated models
    matrices$satmodels.repr <- sapply(seq_len(L),
      function(j) {
        fac <- positionsfac[j, ]
        fac[matrices$indf[j]] <- 0
        digest::digest(fac)
      }
    )

    matrices$q <- p - sum(l) + L #Total number of factors and covariates

  } else { #if there are no factors

    positionsfac <- NULL
    l <- integer(0)
    matrices$indf <- integer(0)
    matrices$satmodels.repr <- character(0)

    matrices$q <- p
  }

  matrices$positionsfac <- positionsfac
  matrices$L <- L
  matrices$l <- l

  return(matrices)
}

#' Builds the matrix for the whole model space and adds an aditional column
#' to associated probabilites
#'
#' @keywords internal
buildmodelsmatrix <- function(q) {
  models.ord <- c(seq_len(2^q-1),0)
  zeros <- numeric(2^q)

  tmodels.mat <- sapply(models.ord,
                        function(j) num2bin.model(j, q, NULL)["bin",])

  cbind(t(tmodels.mat), zeros)
}

#' Builds models matrix when factors for each model with an active factor in
#' Btest functions.
#'
#' @keywords internal
build_ind <- function(k) {
  ind <- t(sapply(2:2^k - 1,
                  FUN = function(j) num2bin.model(j, k, NULL)["bin",]))

  rs <- rowSums(ind)
  ind[!(rs == k | (rs == (k - 1) & ind[, 1])), , drop = FALSE]
}

#' Original code from package 2.6.0 \pkg{BayesVarSel} (distributed under GPL-2),
#' by Gonzalo García-Donato and Anabel Forte.
#'
#' Adapted by Carolina Mulet to fit \pkg{meantobvs}'s code and return a matrix
#' with the binary expression of a model and the active variables with NA
#' to reduce computational burden
#'
#' @keywords internal
num2bin.model <- function(x, p, NAvars) {
  #x is the number to get its binary expression, p is the number of variables,
  #namesxnotnull is the name of the p competing vars and NAvars, the ones with NAs.
  #If NAvars = NULL, the row "bin" equals integer.base.b_C

  if (x == 0) {
    res <- numeric(p)
  } else {
    ndigits <- (floor(logb(x, base = 2)) + 1)
    res <- numeric(ndigits)
    for (i in 1:ndigits) {
      res[i] <- (x %% 2)
      x <- (x %/% 2)
    }

    res <- c(res, numeric(p - ndigits)) #variables active
  }

  resNA <- res * NAvars #variables in model with NA

  matrix(c(res, resNA), byrow = T, nrow = 2, ncol = p,
         dimnames = list(c("bin","NA"), names(NAvars)))
}

#' Original code from package 3.1-0 \pkg{lmerTest} (distributed under GPL-2, GPL-3),
#' by Alexandra Kuznetsova, Per Bruun Brockhoff and Rune Haubo Bojesen Christensen.
#'
#' Adapted by Carolina Mulet to fit \pkg{meantobvs}'s code and obtain just
#' rank defficient matrices.
#'
#' @keywords internal
model.rankdef.matrix <- function (model.frame.aux) {
  #internal function to create rank defficient matrices from a given dataframe
  #created from a model.frame call

  if (ncol(model.frame.aux) == 1) { #just the response
    Xnull.def <- cbind(`(Intercept)` = rep.int(1, nrow(model.frame.aux)))
    return(Xnull.def)
  }

  terms <- attr(terms(model.frame.aux), "term.labels")

  Xi.rdef <- sapply(terms, function(var) {

    f <- as.formula(paste0("~ 0 + ", var))
    #without intercept produces one columns per level
    model.matrix(f, data = model.frame.aux)

  }, simplify = FALSE)

  Xfull.def <- do.call(cbind, Xi.rdef)

  #Always include intercept termn:
  Xfull.def <- cbind(`(Intercept)` = rep.int(1, nrow(Xfull.def)), Xfull.def)
  return(Xfull.def)
}

#' Stable version of log(sum(exp(vector of logarithms))).
#' To avoid infite Bayes factors at the average step.
#'
#' If maximum is not infinite, it returns:
#'  log(sum(exp(x))) = max(x) + log(sum(exp(x - max(x)))),
#' a more stable computational expression for the sum in logarithmic scale.
#'
#' @keywords internal
logsumexp.stable <- function (x) {
  max.x <- max(x)
  if (!is.finite(max.x)) {

    warning("A Bayes factor in infinite.\n", immediate. = TRUE)
    return(max.x)

  } else logsum <- max.x + log(sum(exp(x - max.x)))

  return(logsum)
}

#' Choose corresponding needed functions for a given glm family to compute BFs.
#'
#' @keywords internal
addfamilyfun <- function (family) {

  fam <- family$family
  link <- family$link

  switch(fam,

         gaussian = {
           V1 <- function (mu) 0 #d V(mu)/d mu

           switch (link,

                   identity = { #canonical
                     mu2 <- function(eta) 0 #d2 mu/d eta2
                   },

                   log = {
                     mu2 <- function(eta) family$mu.eta(eta) #d2 mu/d eta2
                   },

                   inverse = {
                     mu2 <- function(eta) 2 / (eta^3) #d2 mu/d eta2
                   },

                   stop(paste0("link ", link, " non-recognized for the ", fam, " family.\n"))
           )

         },

         binomial = {
           V1 <- function (mu) 1 - 2 * mu #d V(mu)/d mu

           switch (link,

                   identity = {
                     mu2 <- function(eta) 0 #d2 mu/d eta2
                   },

                   logit = { #canonical
                     mu2 <- function(eta) {
                       p <- plogis(eta)
                       p * (1 - p) * (1 - 2*p) #d2 mu/d eta2
                     }
                   },

                   probit = {
                     mu2 <- function(eta) - eta * dnorm(eta) #d2 mu/d eta2
                   },

                   cloglog = {
                     mu2 <- function(eta) family$mu.eta(eta) * (1 - exp(eta)) #d2 mu/d eta2
                   },

                   cauchit = {
                     mu2 <- function(eta) -2/pi * eta/sqrt(1 + eta^2) #d2 mu/d eta2
                   },

                   stop(paste0("link ", link, " non-recognized for the ", fam, " family.\n"))
           )

         },

         poisson = {
           V1 <- function (mu) 1 #d V(mu)/d mu

           switch (link,

                   identity = {
                     mu2 <- function(eta) 0 #d2 mu/d eta2
                   },

                   log = { #canonical
                     mu2 <- function(eta) family$mu.eta(eta) #d2 mu/d eta2
                   },

                   sqrt = {
                     mu2 <- function(eta) 2 #d2 mu/d eta2
                   },

                   stop(paste0("link ", link, " non-recognized for the ", fam, " family.\n"))
           )

         },

         Gamma = {
           V1 <- function (mu) 2 * mu #d V(mu)/d mu

           switch (link,

                   identity = {
                     mu2 <- function(eta) 0 #d2 mu/d eta2
                   },

                   log = {
                     mu2 <- function(eta) family$mu.eta(eta) #d2 mu/d eta2
                   },

                   inverse = { #canonical
                     mu2 <- function(eta) 2 / (eta^3) #d2 mu/d eta2
                   },

                   stop(paste0("link ", link, " non-recognized for the ", fam, " family.\n"))
           )

         },

         inverse.gaussian = {
           V1 <- function (mu) 3 * mu^2 #d V(mu)/d mu

           switch (link,

                   `1/mu^2` = { #canonical
                     mu2 <- function(eta) 3/4 / (eta^(5/2)) #d2 mu/d eta2
                   },

                   identity = {
                     mu2 <- function(eta) 0 #d2 mu/d eta2
                   },

                   log = {
                     mu2 <- function(eta) family$mu.eta(eta) #d2 mu/d eta2
                   },

                   inverse = {
                     mu2 <- function(eta) 2 / (eta^3) #d2 mu/d eta2
                   },

                   stop(paste0("link ", link, " non-recognized for the ", fam, " family.\n"))
           )

         },

         quasi = {#currently no available

           switch(fam$varfun,

                  constant = {
                    V1 <- function (mu) 0 #d V(mu)/d mu
                  },

                  `mu(1-mu)` = {
                    V1 <- function (mu) 1 - 2 * mu #d V(mu)/d mu
                  },

                  mu = {
                    V1 <- function (mu) 1 #d V(mu)/d mu
                  },

                  `mu^2` = {
                    V1 <- function (mu) 2 * mu #d V(mu)/d mu
                  },

                  `mu^3` = {
                    V1 <- function (mu) 3 * mu^2 #d V(mu)/d mu
                  },

                  stop(paste0("Variance function non-recognised. ",
                              "It must be one of: constant, mu(1-mu), mu, `mu^2` and `mu^3`.\n"))
            )

           switch (link,

                   identity = {
                     mu2 <- function(eta) 0 #d2 mu/d eta2
                   },

                   logit = {
                     mu2 <- function(eta) family$mu.eta(eta) #d2 mu/d eta2
                   },

                   probit = {
                     mu2 <- function(eta) - eta * dnorm(eta) #d2 mu/d eta2
                   },

                   cloglog = {
                     mu2 <- function(eta) family$mu.eta(eta) * (1 - exp(eta)) #d2 mu/d eta2
                   },

                   cauchit = {
                     mu2 <- function(eta) -2/pi * eta/sqrt(1 + eta^2) #d2 mu/d eta2
                   },

                   inverse = {
                     mu2 <- function(eta) 2 / (eta^3) #d2 mu/d eta2
                   },

                   sqrt = {
                     mu2 <- function(eta) 2 #d2 mu/d eta2
                   },

                   `1/mu^2` = {
                     mu2 <- function(eta) 3/4 / (eta^(5/2)) #d2 mu/d eta2
                   },

                   stop(paste0("link ", link, " non-recognized for the ", fam, " family.\n"))
           )

         },

         quasibinomial = {
           V1 <- function (mu) 1 - 2 * mu #d V(mu)/d mu

           switch (link,

                   identity = {
                     mu2 <- function(eta) 0 #d2 mu/d eta2
                   },

                   logit = { #canonical
                     mu2 <- function(eta) family$mu.eta(eta) #d2 mu/d eta2
                   },

                   probit = {
                     mu2 <- function(eta) - eta * dnorm(eta) #d2 mu/d eta2
                   },

                   cloglog = {
                     mu2 <- function(eta) family$mu.eta(eta) * (1 - exp(eta)) #d2 mu/d eta2
                   },

                   cauchit = {
                     mu2 <- function(eta) -2/pi * eta/sqrt(1 + eta^2) #d2 mu/d eta2
                   },

                   stop(paste0("link ", link, " non-recognized for the ", fam, " family.\n"))
           )

         },

         quasipoisson = {
           V1 <- function (mu) 1 #d V(mu)/d mu

           switch (link,

                   identity = {
                     mu2 <- function(eta) 0 #d2 mu/d eta2
                   },

                   log = { #canonical
                     mu2 <- function(eta) family$mu.eta(eta) #d2 mu/d eta2
                   },

                   sqrt = {
                     mu2 <- function(eta) 2 #d2 mu/d eta2
                   },

                   stop(paste0("link ", link, " non-recognized for the ", fam, " family.\n"))
           )

         },

         stop(paste0(" family ", fam, " non-recognized.\n"))
  )

  #if dispersion fixed, get it; otherwise compute the Pearson estimation
  family$dispersionfun <- ifelse(
    is.na(family$dispersion),
    function (df.residual, weights, residuals) {
      if (df.residual > 0)
        sum((weights * residuals^2)[weights > 0])/df.residual
      else NA_real_},
    function (...) family$dispersion)

  #Add functions to the family object
  family$mu2 <- mu2 #second derivative of the inverse-link function with respect to eta
  family$V1 <- V1 #derivative of the variance function with respect to mu

  family
}

#' Compute observed Fisher information for a fitted glm object.
#'
#' @keywords internal
Fobs.info <- function (glmfit, family = addfamilyfun(glmfit$family)) {

  #glm fitted objects
  eta <- glmfit$linear.predictors #the linear fit on link scale
  mu <- glmfit$fitted.values
  res <- glmfit$y - mu

  #derivatives of the inverse-link function with respect to the linear predictor
  mu1 <- family$mu.eta(eta)
  mu2 <- family$mu2(eta)

  #variance as a function of the mean
  V <- family$variance(mu)
  V1 <- family$V1(mu) #derivative

  #estimated or fixed dispersion parameter
  dispersion <- family$dispersionfun(glmfit$df.residual, glmfit$weights, res / mu1)

  #Observed Fisher information:
  d <- glmfit$prior.weights / dispersion *
    (mu1^2 / V - res * mu2 / V + res * mu1^2 * V1 / V^2)

  return(d)
}

#' Obtain list of parameters for controlling the fitting process on fastglm
#' from a list control in glm.control() format. By default, method selected is
#' LLT Cholesky decomposition, faster.
#'
#' @keywords internal
glm2fastglm.control <- function (control = glm.control(), method = 2L) {

  tol <- control$epsilon
  maxit <- control$maxit

  fastglm::fastglm_control(fastmethod = method, tol, maxit)

}

#' Returns a string with the family and link function on a family object
#'
#' @keywords internal
family.link <- function(family) {
  paste(family$family, family$link, sep = "/")
}

#' Get computational time from an starting event
#'
#' @keywords internal
get.time <- function (start.time) {

  t <- difftime(Sys.time(), start.time, units = "secs")

  as.numeric(t)
}

#' Show computational time
#'
#' @keywords internal
format.time <- function(x) {
  if (!is.finite(x) || x < 0) return("--:--:--")

  x <- floor(x)

  # h <- x %/% 3600
  # m <- (x %% 3600) %/% 60
  # s <- x %% 60
  h <- floor(x / 3600)
  if (h >= 100) return("99:59:59")
  m <- floor((x - h * 3600) / 60)
  s <- x - h * 3600 - m * 60

  if (any(is.na(c(h, m, s)))) return("--:--:--")


  sprintf(sprintf("%02.0f:%02.0f:%02.0f", h, m, s))
}

#' Personalized progress bar.
#'
#' @keywords internal
upd.progress <- function(i, total, start.time, width = 20) {

  #current computation time
  elapsed <- get.time(start.time)
  #models rate computed
  mrate <- if (elapsed > 0) i / elapsed else 0
  #update each ceiling(mrate / 4) models (approx 4 times per second)
  update <- i %% max(1, ceiling(mrate / 4)) == 0 || i == total

  if (update) { #check whether or not to update progress bar
    #estimated remaining time
    eta <- if (mrate > 0) (total - i) / mrate else NA_real_
    #progress percent
    percent <- i / total * 100

    #build bar
    numcomplete <- floor(width * i / total)
    bar <- paste0(
      strrep("=", max(0, numcomplete - 1)),
      if (i < total) ">" else "=",
      strrep("-", max(0, width - numcomplete))
    )

    #write bar and info
    cat(
      sprintf(
        "\r\033[K[%s] %3.0f%% [Elapsed time: %s | Estim. time: %s | ≈ %.0f models/s]",
        bar,
        percent,
        format.time(elapsed),
        format.time(eta),
        mrate
      )
    )

    flush.console()
  }

}

#' !%in% function
#'
#' @keywords internal
"%notin%" <- function(x, table) {

  match(x, table, nomatch = 0) == 0

}
