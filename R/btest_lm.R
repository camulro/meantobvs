#' Bayes Factors and Posterior Probabilities for Linear Regression Models with Missing Data
#'
#' The model space is build from a list of linear regression models proposed to explain
#' a common response. It returns the Bayes factors and posterior probabilities computed
#' averaging over a set of imputations in the presence of missing data. If factors,
#' a weighted average Bayes factor over the submodel space is also performed to
#' account for uncertainty of factor levels.
#'
#' Given a list of competing models, the model space is made up by them, assuming that the
#' intercept term is present in every model. The simplest one M0, can be specified (\code{null.model})
#' and must be nested in the rest. Hence, the posterior distribution over the model space
#' is given through Bayes' theorem:
#'
#' Pr(Mi | \code{data}) = Pr(Mi) * AvBi / C,
#'
#' where AvBi is the average Bayes factor of Mi to M0 under missing data,
#' Pr(Mi) is the prior probability of Mi and C is the normalizing constant.
#' AvBi is the average of the \code{n.imp} data-driven Bayes factors (BF):
#'
#' AvBi = 1/\code{n.imp} * (Bi(1) + ... + Bi(\code{n.imp})),
#'
#' where Bi(j) corresponds to the BF for model Mi to M0 under the jth imputed dataset.
#' Data-driven BF can be either computed using popular g-prior choices or approximated
#' with the BIC (Schwarz, 1978), and the default choice, or the test-based BF (TBF)
#' (Held, Gravestock and Sabanés, 2015) with the \code{BF.method} argument.
#'
#' If the BF computation method chosen is \code{"gprior"}, data-driven BFs depend on
#' the prior assigned for the model-specific parameters given by \code{prior.betas}.
#' This method is faster using \pkg{BayesVarSel} functions. Choices currently available are:
#' -"Robust" and denotes the criteria-based prior of Bayarri, Berger, Forte and
#' Garcia-Donato (2012).
#' -"gZellner" is the default option and corresponds to the prior in Zellner (1986)
#' with g=n fixed.
#' -"Liangetal" prior is the hyper-g/n of Liang et al (2008) with a=3.
#' -"ZellnerSiow" is the multivariate Cauchy prior by Zellner and Siow (1980, 1984).
#' -"FLS" corresponds to the prior in Zellner (1986) with g=max(n, p*p) fixed, the
#' (benchmark) prior recommended by Fernandez, Ley and Steel (2001).
#' -"intrinsic.MGC" is the intrinsic prior derived by Moreno, Giron, Casella (2015).
#' -"IHG" corresponds to the intrinsic hyper-g prior derived in Berger, Garcia-Donato,
#' Moreno and Pericchi (2022).
#'
#' If the BF computation method chosen is \code{"TBF"}, \code{prior.betas} also determines
#' the formula of the BF as before. In this case, the choices available are: "gZellner",
#' "FLS", "Liangetal" and the adapted version of the ZellnerSiow" prior for TBF.
#'
#' The prior over the model space Pr(Mi) offers three options through \code{prior.models}:
#' -"Constant" assigns the same prior probability to every model, default one.
#' -"ScottBerger" assigns the same prior probability to every different model size.
#' -"User": if chosen, user has to provide a N dimensional vector (where N the number of
#' competing models) with the model prior probabilities of each one in \code{models}
#' through \code{priorprobs}.
#'
#' In the presence of factors, in order to account for uncertainty of factor levels,
#' Pr(Mi) is factorized following García-Donato and Paulo (2022), as follows. Each factor is
#' represented by its rank deficient dummy parametrization, making each model with a factor
#' active the representative of the model space given by the 2^l-l different submodels with
#' at least one dummy active, where l is the number of levels. This approach derives in a
#' hierarchical prior where the part over covariates and factors is a standard model
#' prior given by \code{prior.models}. The prior over the different submodels
#' assumes prior independence between factors, whereas within a factor two options are available:
#' "Constant" assigns the same prior to every submodel and "ScottBerger" to each model size,
#' which is the default choice. A direct selection over the dummy variables can
#' also be performed through \code{marginal.factors}.
#'
#' If there are missing data, a multiple imputation procedure has to be performed.
#' It can be done internally, using the imputation related arguments available, or
#' externally through an array of user-given imputed datasets with \code{imp.datasets}.
#' First path allows the user to specify the imputation model of each variable
#' with a binary named matrix \code{imp.predict.mat} for the method or methods in
#' \code{imp.mice.method}. Parallelized imputation is recommended for big \code{data}.
#' Second path is the recommended one since it allows the user a fully data-designed
#' imputation method.
#'
#' @export
#' @param data Data frame containing the data.
#' @param models List with the entertained models and their defining formulas, with one
#' nested in all the others. If the list is unnamed, default names are given.
#' @param null.model String for the name of the null model on \code{models}. By default,
#' the names of variables are used to identify the null. If provided, the string
#' must coincide with the one with the largest sum of squared errors and should
#' be the one with the smallest size.
#' @param BF.method Method used to compute or approximate data-driven Bayes factors
#' (to be literally specified). Possible choices include "BIC", "TBF" and "gprior"
#' (see details).
#' @param prior.betas Prior distribution for model coefficients (to be literally ç
#' specified). Possible choices are: "Robust", "Liangetal",
#' "gZellner", "ZellnerSiow", "FLS", "intrinsic.MGC" and "IHG" (see details).
#' @param prior.models Model prior distribution over the covariates and/or factors
#' model space (to be literally specified). Possible choices are "Constant",
#' "ScottBerger" and "User" (see details).
#' @param prior.models.dummies Prior distribution over the dummies submodel space
#' built by active factors (to be literally specified). Possible choices are
#' "Constant" and "ScottBerger" (see details).
#' @param marginal.factors Logical to indicate whether or not to marginalize factors'
#' probabilities following García-Donato and Paulo (2022). By default, it is set to TRUE.
#' @param priorprobs A N dimensional vector (being N the number of competing models)
#' defining the prior model probabilities for each one in \code{models} (if
#' \code{prior.models}= "User"; see details).
#' @param imp.mice.method If missing data, it specifies the method for \pkg{mice}'s
#' imputation. Can be either a string or a string vector for the competing variables.
#'  By default, predictive mean matching is performed.
#' @param imp.predict.mat Named binary matrix with the competing variables in rows.
#' Each entry equals 1 if the column variable is used as a predictor. Order
#' in columns defines the imputation visit sequence. By default, imputation models
#' are based on correlations.
#' @param n.imp Number of imputed datasets for model posterior computation. By default
#' it is set to 39.
#' @param maxit Number of iterations for \pkg{mice}'s imputation. By default, it is 5.
#' @param parallelmice Logical to indicate whether or not to use parallelization on
#' \code{\link[mice]{mice}}'s imputation. By default, automatically performs it if the
#' number of imputations or competing variables given by \code{formula} are big enough.
#' @param n.core Number of cores for parallel imputation.
#' @param imp.datasets Array or list for imputed datasets if given by user. By default
#' it is set to NULL and imputation is performed following other imputation arguments.
#' @param imp.seed Seed for imputation.

#' @return \code{\link[meantobvs]{btest.lm}} returns an object of type
#' \code{meantobtest} with the following elements:
#' \item{lBFi0}{Bayes factors in logaritmic scale of each model to the null}
#' \item{PostProbi}{Posterior probabilities for each model in \code{models}}
#' \item{models}{List with the entertained models.}
#' \item{nullmodel}{Name in \code{models} of the null (simplest) model}
#' \item{modelspool}{If missings, list of the combined estimates for each model
#' in \code{models} fitted by \code{\link[stats]{lm}} over the \code{n.imp} imputed
#' datasets; or \code{lm} object when there are no missings}
#' \item{positions}{Matrix with L rows and p1 * (sum_j l_j - L), where p1 is the number
#' of covariates, L the number of factors and l_j the number of levels of the jth factor,
#' with 1 if the column dummy makes up the row factor and 0 otherwise (when relevant)}
#' \item{positionsx}{Logical vector of length p indicating whether or not the
#' variable is a numerical covariate (when relevant)}
#' \item{imp.info}{List of arguments used for the imputation step and other
#' information V}
#' \item{compress.imp.array}{Compressed array of imputed datasets (when relevant)}
#' \item{BF.method}{Method used to compute data-driven Bayes factors}
#' \item{prior.betas}{Chosen \code{prior.betas} argument}
#' \item{prior.models}{Two-dimensional vector with \code{prior.models} and
#' \code{prior.models.dummies} chosen. If there are no factors or \code{marginal.factors}
#' is set to FALSE, it saves the only argument used, \code{prior.models}}
#' \item{marginal.factors}{Logical to indicate whether or not are marginalized factors'
#' model space probabilities such as García-Donato and Paulo (2022)}
#' \item{priorprobs}{Prior probabilities over the true model size}
#' \item{call}{The \code{call} to the function}
#'
#' @author Carolina Mulet and Gonzalo García-Donato
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[meantobvs]{bvs.lm}} for an exact computation
#' of the model posterior distribution (recommended when p<20).
#'
#' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
#'
#' García-Donato, G. and Forte, A. (2018) Bayesian Testing,
#' Variable Selection and Model Averaging in Linear Models using R with
#' BayesVarSel. The R Journal. 10: 329.
#'
#' Garcia-Donato, G. and Paulo, R. (2022)<DOI:10.1080/01621459.2021.1889565>
#' Variable Selection in the Presence of Factors: A Model Selection Perspective.
#' Journal of the American Statistical Association. 117. 1-27.
#'
#' Bayarri, M.J., Berger, J.O., Forte, A. and Garcia-Donato, G.
#' (2012)<DOI:10.1214/12-aos1013> Criteria for Bayesian Model choice with
#' Application to Variable Selection. The Annals of Statistics. 40: 1550-1557.
#'
#' Scott, J.G. and Berger, J.O. (2010) Bayes and empirical-Bayes multiplicity
#' adjustment in the variable-selection problem. The Annals of Statistics.
#' 38: 2587–2619.
#'
#' Zellner, A. (1986)<DOI:10.2307/2233941> On Assessing Prior Distributions and
#' Bayesian Regression Analysis with g-prior Distributions. In Bayesian
#' Inference and Decision techniques: Essays in Honor of Bruno de Finetti (A.
#' Zellner, ed.) 389-399. Edward Elgar Publishing Limited.
#'
#' Schwarz, G. (1978) Estimating the dimension of a model. The Annals of
#' Statistics. 6: 461–464.
#'
#' Held L, Sabanés Bové D, Gravestock I (2015).<DOI:10.1214/14-STS510> Approximate
#' Bayesian Model Selection with the Deviance Statistic. Statistical Science. 30.
#'
#' van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice: Multivariate Imputation
#' by Chained Equations in R. Journal of Statistical Software. 45: 1–67.
#'
#' @examples
#' \donttest{
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' data("dataS97")
#'
#' models.list <- list(
#'   M0 = gr56092 ~ 1,
#'   M1 = gr56092 ~ lifee060,
#'   M2 = gr56092 ~ gdpsh60l,
#'   M3 = gr56092 ~ p60,
#'   M7 = gr56092 ~ lifee060 + gdpsh60l + p60
#' )
#'
#' #Few imputations for simplicity, real analyses need more.
#' dataS97.meantobtest <- btest.lm(
#'   data = dataS97, models = models.list, n.imp = 2, imp.seed = 1
#' )
#'
#' #Show the results:
#' dataS97.meantobtest
#' dataS97.meantobtest$lBFi0
#' dataS97.meantobtest$PostProbi
#' }
#'

btest.lm <- function (data,
                      models,
                      null.model = NULL,
                      BF.method = ifelse(is.null(prior.betas),
                                        "BIC", "gprior"),
                      prior.betas = NULL,
                      prior.models = "Constant",
                      prior.models.dummies = "ScottBerger",
                      marginal.factors = TRUE,
                      priorprobs = NULL,
                      imp.mice.method = "pmm",
                      imp.predict.mat = NULL,
                      n.imp = 039,
                      maxit = 5,
                      parallelmice = NULL,
                      n.core = NULL,
                      imp.datasets = NULL,
                      imp.seed = runif(1,0,09011975)) {

  #N is the number of models:
  N <- length(models)

  #Check btest given arguments
  btest.args <- checkbtestarguments(models, null.model, N)

  #Define list of arguments for posterior computation
  model.context <- list(
    relax.nest = btest.args$relax.nest,
    models = btest.args$models
  )

  SSE <- numeric(N) #SSEs for each model
  Dim <- rep.int(0L, N)
  mt <- list() #list of terms for each model

  #goes through all competing models and saves some results
  covar.list <- list() #list that contains the names of the variables in each model
  compvars <- c() #name of original competing vars
  for (i in seq_len(N)) {
    f <- as.formula(btest.args$models[[i]])
    compvars <- c(compvars, attr(terms(f), "term.labels"))
    temp <- lm(formula = f, data = data, y = TRUE, x = TRUE)

    SSE[i] <- as.numeric(crossprod(temp$residuals))

    Xi <- model.rankdef.matrix(model.frame(temp))
    covar.list[[i]] <- dimnames(Xi)[[2]]
    Dim[i] <- length(covar.list[[i]])
    mt[[i]] <- temp$terms
  }
  ordered.SSE <- sort(SSE, index.return = TRUE, decreasing = TRUE)
  #Which one acts as null model:
  nullmodel.pos <- ordered.SSE$ix[1]

  model.context$nullmodel.pos <- nullmodel.pos
  model.context$covar.list <- covar.list

  #Check null model
  if (btest.args$relax.nest &&
      nullmodel.pos != btest.args$nullmodel.posuser){
      stop("The given null model does not coincide with the one with the\n",
           "largest sum of squared error (and it should).\n")
  }
  #change the string for the formula and specify models to compute BF
  null.model <- as.formula(btest.args$models[[nullmodel.pos]])

  #check null model
  if (attr(terms(null.model), "intercept") == 0) {
    stop("Null model shall contain the intercept term.",
         " Please, provide a proper null model.\n")
  }

  competing.models <- seq_len(N)[-nullmodel.pos]

  model.context$competing.models <- competing.models

  #Competing vars full formula:
  full.formula <- update(null.model, paste0(". ~ ",
                                     paste(unique(compvars), collapse = " + ")))

  #Build matrices and objects needed later on
  matrices <- buildmatrices(full.formula, null.model, data, marginal.factors)

  model.context$namesxnotnull <- matrices$namesxnotnull
  model.context$namesnull <- matrices$namesnull

  Dim <- Dim - matrices$p0 #model dimension (without fixed vars)

  #The response variable
  obsnotNA <- rownames(matrices$X0)
  y <- matrices$framenull[obsnotNA, 1] #response variable without missings
  n <- length(y)
  SS0 <- SSE[nullmodel.pos]

  #Compute model prior
  lprior.models <- priormodels.btest(prior.models, N, Dim, priorprobs)

  #Check approx method and priors chosen and define the function to be used
  lBF <- checkforprior.betas.lm(
    BF.method = BF.method, prior.betas = prior.betas,
    n = n, p = max(Dim), p0 = matrices$p0, y = y, SS0 = SS0
  )

  matrices$X.full <- matrices$X.full[obsnotNA, , drop = FALSE]

  #check for missings and define variables with NAs
  NAvars <- checkformissings(
    y = matrices$framenull[, 1], matrices$framenull[, -1], matrices$X.full
  )

  anyNAvar <- sum(NAvars) > 0
  if(!anyNAvar) n.imp <- 0 #if no NAs, no imputation needed

  #Imputation step
  if (anyNAvar) {
    if (is.null(imp.datasets)) { #if there are no given imputations, build them
      imputation <- buildimputation(
        NAvars, full.formula, data, imp.predict.mat, n.imp, maxit, n,
        matrices$q, matrices$p0, imp.mice.method, imp.seed, parallelmice,
        n.core, obsnotNA, matrices$ordvars
      )
    } else {
      imputation <- extimputation(
        full.formula, imp.datasets, n0 = dim(data)[1], matrices$framefull,
        matrices$ordvars, obsnotNA, matrices$p0, NAvars
      )
      n.imp <- imputation$n.imp
    }
  }

  if (anyNAvar && n.imp > 1) {
    #function to compute log(BFa0) for a given model as an average of BF computed
    #by BF.method over the imputed datasets
    lBF.method <- function(model) lBF.av(
      model, imputation.array = imputation$imputation.array,
      lBF = lBF, p0 = matrices$p0, n.imp = n.imp
    )
  } else {#if (anyNAvar) {
    lBF.method <- function(model) lBF(
      k = length(model),
      X = imputation$imputation.array[, c(seq_len(matrices$p0), model + matrices$p0), ]
    )
  } # else lBF.method <- function(model) lBF(
  #   k = length(model),
  #   X = cbind(matrices$X0, matrices$X.full[, model, drop = FALSE])
  # )

  lBFmodel <- function (model) lBF(k = length(model),
                                   X = cbind(matrices$X0, matrices$X.full[, model]))

  mF <- matrices$L > 0 && marginal.factors
  positionsfac <- if (mF) matrices$positionsfac else NULL
  indf <- if (mF) matrices$indf else NULL

  #Check if factors present and if marginalization of their probabilities.
  #Define model prior and BF
  lBF.comp <- BFcomp.btest(
    lprior.models, prior.models.dummies, mF, matrices, NAvars, lBF.method, lBFmodel
  )

  #Posterior computation of model space defined by models list
  posterior <- posterior.btest(model.context, lprior.models, lBF.comp)

  #Evaluate lm of each model with missings using Rubin's rule
  modelspool <- list()
  for(j in competing.models){
    namesj <- which(matrices$namesxnotnull %in% covar.list[[j]])
    if (sum(NAvars[namesj]) > 0) {
      fit <- list()
      for (i in 1:n.imp) {
        if (mF) {#remove last dummy for each factor, first q0 vars are the fixed ones
          Xi <- imputation$imputation.array[,
            c(seq_len(matrices$p0), setdiff(namesj, matrices$indf) + matrices$p0), i
          ]
        } else Xi <- imputation$imputation.array[,
          c(seq_len(matrices$p0), namesj + matrices$p0), i
        ]

        z <- lm.fit(x = Xi, y = y)
        z$terms <- mt[[j]]; class(z) <- "lm"; fit[[i]] <- z
      }
      modelspool[[j]] <- mice::pool(fit)
      modelspool[[j]]$call <- NULL #otherwise, Rstudio returns a warning trying to read modelspool[[j]]$call
    } else modelspool[[j]] <- lm(btest.args$models[[j]], data)
  }
  modelspool[[nullmodel.pos]] <- lm(null.model, data)
  names(modelspool) <- names(btest.args$models)

  result <- list()
  result$lBFi0 <- posterior$lBFi0
  result$PostProbi <- posterior$PostProbi
  result$models <- btest.args$models
  result$nullmodel <- names(btest.args$models)[nullmodel.pos]
  result$modelspool <- modelspool

  if (mF) {
    #matrix for the factors index
    result$positions <- matrices$positionsfac
    result$positionsx <- matrices$positionsx
  }

  if (anyNAvar) {
    #arguments used for imputation
    result$imp.info <- imputation$imp.info

    # save the imputed datasets for sensitivity analysis
    raw.imp.array <- serialize(imputation$imputation.array, NULL)
    result$compress.imp.array <- memCompress(raw.imp.array, type = "xz")
  }

  result$BF.method <- BF.method #method used for BF computation
  if (is.null(prior.betas) & BF.method %in% c("gprior", "TBF")) prior.betas <- "gZellner"
  if (BF.method == "BIC") prior.betas <- NULL #does not apply with BIC
  result$prior.betas <- prior.betas
  if (mF) {
    result$prior.models <- c(prior.models, prior.models.dummies)
  } else result$prior.models <- prior.models
  result$marginal.factors <- marginal.factors #whether or not factors are marginalized
  result$priorprobs <- exp(posterior$lPriorModels)
  result$call <- match.call()

  class(result) <- "meantobtest"

  return(result)
}

#' @keywords internal
posterior.btest <- function (model.context, lprior.models, lBF.comp) {
  #posterior computation of model space defined by models list

  list2env(model.context, envir = environment())
  nmodels <- names(models)

  lBFi0 <- lPriorModels <- numeric(length(nmodels))
  for (i in competing.models){
    modeli <- namesxnotnull %in% covar.list[[i]] #all dummies active if factor active
    names(modeli) <- namesxnotnull

    #check whether the null is nested in the other ones
    if (relax.nest & any(namesnull %notin% covar.list[[i]])) {
      stop("The simplest (null) model may not be nested in all the others.",
           " Please define explicitly the null model if it is the case.\n")
    }

    lBFi0[i] <- lBF.comp(modeli, i) #log(BF_a0)
    lPriorModels[i] <- lprior.models(i) #log-model prior
  }
  cat("\n")

  lPriorModels[nullmodel.pos] <- lprior.models(nullmodel.pos)
  lBFi0[nullmodel.pos] <- 0

  #Compute posterior probabilities
  lBF.PM <- lBFi0 + lPriorModels
  logC <- logsumexp.stable(lBF.PM) # C <- sum(exp(lBFi0 + lPriorModels))
  PostProbi <- exp(lBF.PM - logC)

  names(lBFi0) <-
    paste(nmodels, ".to.", nmodels[nullmodel.pos], sep = "")
  names(PostProbi) <- nmodels
  names(lPriorModels) <- nmodels

  res <- list(lBFi0 = lBFi0)
  res$PostProbi <- PostProbi
  res$lPriorModels <- lPriorModels

  return(res)
}

#' @keywords internal
BFcomp.btest <- function (lprior.models, prior.models.dummies, mF,
                          matrices, NAvars, lBF.method, lBF) {
  #Check arguments and define the functions to compute Bayes factors

  positionsfac <- matrices$positionsfac
  namesxnotnull <- matrices$namesxnotnull

  if (mF) {
    #Check model priors for dummies chosen and define the function to be used
    if (prior.models.dummies %notin% c("ScottBerger", "Constant")) {
      stop("Only priors 'ScottBerger' and 'Constant' supported.\n")
    }
    switch (prior.models.dummies, #change the string for the corresponding function
            Constant = {lprior.models.dummies <-
              function (d, df) {-sum(log(2^(df) - 1 - df))}},
            ScottBerger = {lprior.models.dummies <-
              function (d, df) {-sum(mylchoose(df, d)) - sum(log(df - 1))}}
    )

    #Build submodels for each possible factor
    submodels.matrix <- lapply(matrices$l, build_ind)
    for (i in 1:matrices$L) {
      colnames(submodels.matrix[[i]]) <- names(which(positionsfac[i,] == 1))
    }

    #Define BF to compute the marginal for the dummies
    lBF.comp <- function (modeli, i) {
      #modeli logical, all dumies active if its factor is active to build submodels
      d <- as.vector(positionsfac %*% modeli) #levels of factors
      f <- d > 0 #active factors
      df <- d[f] #levels of active factors

      areNA <- sum(modeli * NAvars) > 0
      if (any(f)) {
        #Merge submodel matrices to build submodel space
        submodels <- Reduce(function(x, y) merge(x, y, by = NULL),
                            submodels.matrix[f]) == 1
        submod.names <- colnames(submodels)

        lBF.d <- lprior.d <- numeric(nrow(submodels))
        for (j in 1:nrow(submodels)) { #go throug all submodels to compute BF
          submodj <- submodels[j,]
          dj <- as.vector(positionsfac[which(f), submod.names] %*% submodj)

          current.model <- modeli; current.model[submod.names] <- submodj

          #check if there are NAs in the model considered to save computation time
          if (areNA) {
            lBF.d[j] <- lBF.method(model = which(current.model)) #log(BF_a0)
          } else { #if there are no missings, compute the BF by the method selected
            lBF.d[j] <- lBF(model = which(current.model)) #log(BF_a0)
          }
          lprior.d[j] <- lprior.models.dummies(dj, df) #log(Pr(M_delta))
        }

        lBFi0 <- logsumexp.stable(lBF.d + lprior.d)
      } else {

        #check if there are NAs in the model considered to save computation time
        if (areNA) {
          lBFi0 <- lBF.method(model = which(modeli))
        } else { #if there are no missings, compute the BF by the method selected
          lBFi0 <- lBF(model = which(modeli))
        }

      }
      return(lBFi0)
    }

  } else { #Define the standard BF
    lBF.comp <- function (modeli, i) {
      #check if there are NAs in the model considered to save computation time
      if (sum(modeli * NAvars) > 0) {
        lBFi0 <- lBF.method(model = which(modeli))
      } else { #if there are no missings, compute the BF by the method selected
        lBFi0 <- lBF(model = which(modeli))
      }
      return(lBFi0)
    }
  }
  return(lBF.comp)
}

#' @keywords internal
priormodels.btest <- function (prior.models, N, Dim, priorprobs) {
  #Check arguments and define the functions to compute prior model probabilities

  #Check model priors chosen and define the function to be used
  if (prior.models %notin% c("ScottBerger", "Constant", "User")) {
    stop("Only priors 'ScottBerger', 'Constant' and 'User' supported.\n")
  }

  switch (prior.models, #change the string for the corresponding function
          Constant = {lprior.models <- function (modeli) -log(N)},
          ScottBerger = {lprior.models <- function (modeli)
            -log(length(unique(Dim))) - log(sum(Dim == Dim[modeli]))},
          User = {
            if (is.null(priorprobs)) stop("User prior selected but no prior probabilities provided.\n")
            if (!is.numeric(priorprobs)) stop("User prior selected but no numeric probabilities provided.\n")
            if (any(is.na(priorprobs))) stop("User prior selected but some prior probabilities not provided.\n")
            if (length(priorprobs) != N) stop("Vector of prior probabilities with incorrect length.\n")
            if (sum(priorprobs < 0) > 0) stop("Prior probabilities must be positive.\n")
            if (all(priorprobs == 0)) stop("Prior probabilities must be positive.\n")

            lprior.models <- function(modeli) log(priorprobs[modeli]/sum(priorprobs))}
  )

  return(lprior.models)
}

#' @keywords internal
checkbtestarguments <- function (models, null.model, N) {
  #check arguments
  if (!is.list(models)) stop("Argument models should be a list.\n")

  if (length(models) < 2) stop("At least two models have to be provided.\n")

  #If competing models come wihtout a name, give one by default:
  if (is.null(names(models))){
    if (!is.null(null.model)) stop("Please provide a name for the competing models.\n",
                                   "The null model must be in that list.\n")
    names(models) <- paste("model", seq_len(N), sep="")
  }

  #Check if the given null model is one of the competing models:
  if (!is.null(null.model)){

    #check null model
    if (attr(terms(as.formula(models[[null.model]])), "intercept") == 0) {
      stop("Null model shall contain the intercept term.",
           " Please, provide a proper null model.\n")
    }

    relax.nest = TRUE
    nullmodel.posuser <- which(null.model == names(models))
    if (length(nullmodel.posuser) == 0) {
      stop("The null model provided is not in the list of competing models.\n")
    }

    return(list(models = models, relax.nest = relax.nest,
                nullmodel.posuser = nullmodel.posuser))
  } else relax.nest = FALSE

  return(list(models = models, relax.nest = relax.nest))
}

#' Print an object of class \code{meantobtest}
#'
#' Print an object of class \code{meantobtest}
#' @export
#' @param meantobtest.object Object of class meantobtest.
#' @param ... Additional parameters to be passed.
#'
#' @author Gonzalo Garcia-Donato
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso See \code{\link[meantobvs]{btest.lm}},
#' \code{\link[meantobvs]{btest.glm}} and
#' \code{\link[meantobvs]{btestGD25}} for creating objects of the class
#' \code{meantobtest}.
#'
#' @examples
#' \donttest{
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#' data("dataS97")
#'
#' # Use a small list of named models and two imputations for this example.
#' models.list <- list(M0 = gr56092 ~ 1, M1 = gr56092 ~ lifee060,
#'   M2 = gr56092 ~ gdpsh60l, M3 = gr56092 ~ p60, M4 = gr56092 ~ lifee060 + p60,
#'   M5 = gr56092 ~ lifee060 + gdpsh60l, M6 = gr56092 ~ p60 + gdpsh60l,
#'   M7 = gr56092 ~ lifee060 + gdpsh60l + p60)
#' dataS97.meantobtest <- btest.lm(
#'   data = dataS97, models = models.list, n.imp = 2, imp.seed = 1
#' )
#'
#' #Show the results:
#' dataS97.meantobtest
#' }
#'
print.meantobtest <- function(meantobtest.object,...){
  if (!inherits(meantobtest.object, "meantobtest")){
    warning("An object of class meantobtest is needed.\n")
  }

  cat("-------\n")
  cat("Competing models:\n")
  print(meantobtest.object$models)
  cat("-------\n")
  cat("log Bayes factors (expressed in relation to ",
      meantobtest.object$nullmodel,")\n", sep="")
  print(round(meantobtest.object$lBFi0, 3))
  cat("-------\n")
  cat("Posterior probabilities:\n")
  print(round(meantobtest.object$PostProbi,3))
  cat("\n\n")
}
