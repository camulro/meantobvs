#' Bayes Factors and Posterior Probabilities for Generalized Linear Models with Missing Data
#'
#' The model space is build from a list of linear regression models proposed to explain
#' a common response. It returns the Bayes factors and posterior probabilities computed
#' averaging over a set of imputations in the presence of missing data. If factors,
#' a weighted average Bayes factor over the submodel space is also performed to
#' account for uncertainty of factor levels.
#'
#' GGiven a list of competing models, the model space is made up by them, assuming that the
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
#' This method corresponds to the CHIC g-prior developed by Li and Clyde (2018)
#' when the dispersion parameter is fixed, and adapted through an empirical Bayes
#' approach otherwise. Choices currently available are:
#' -"Robust" denotes the criteria-based prior of Bayarri, Berger, Forte and
#' Garcia-Donato (2012).
#' -"gZellner" is the default option and corresponds to the prior in Zellner (1986)
#' with g=n fixed.
#' -"Liangetal" prior is the adapted hyper-g/n of Liang et al (2008) with a=3.
#' -"ZellnerSiow" is the adapted multivariate Cauchy prior by Zellner and Siow (1980, 1984).
#' -"FLS" corresponds to the prior in Zellner (1986) with g=max(n, p*p) fixed, the
#' (benchmark) prior recommended by Fernandez, Ley and Steel (2001).
#' -"intrinsic.WNC" is the intrinsic prior derived by Womack, Novelo and Casella (2014).
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
#' @param family String, function or the call to a family function among
#' \code{\link[stats]{family}} to specify the error distribution and link function to be
#' used in the model.
#' @param null.model String for the name of the null model on \code{models}. By default,
#' the names of variables are used to identify the null. If provided, the string
#' must coincide with the one with the largest sum of squared errors and should
#' be the one with the smallest size.
#' @param BF.method Method used to compute or approximate data-driven Bayes factors
#' (to be literally specified). Possible choices include "BIC", "TBF" and "gprior"
#' (see details).
#' @param prior.betas Prior distribution for model coefficients (to be literally
#' specified). Possible choices are: "Robust", "Liangetal",
#' "gZellner", "FLS" and "intrinsic.WNC" (see details).
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
#' @param weights NULL or numeric vector of the same length as \code{y} to
#' specify the weights to be used in the glm fitting process.
#' @param offset NULL or a numeric vector of the same length as \code{y} to
#' specify an a priori known component included in the glm fitting process.
#' @param control List of parameters for controlling the glm fitting process.
#' It is set to \code{[stats]{glm.control()}} by default.

#' @return \code{\link[meantobvs]{btest.glm}} returns an object of type
#' \code{btest} with the following elements:
#' \item{lBFi0}{Bayes factors in logaritmic scale of each model to the null}
#' \item{PostProbi}{Posterior probabilities for each model in \code{models}}
#' \item{models}{List with the entertained models.}
#' \item{nullmodel}{Name in \code{models} of the null (simplest) model}
#' \item{modelspool}{If missings, list of the combined estimates for each model
#' in \code{models} fitted by \code{\link[stats]{glm}} over the \code{n.imp} imputed
#' datasets; or \code{glm} object when there are no missings}
#' \item{positions}{Matrix with L rows and p1 * (sum_j l_j - L), where p1 is the number
#' of covariates, L the number of factors and l_j the number of levels of the jth factor,
#' with 1 if the column dummy makes up the row factor and 0 otherwise (when relevant)}
#' \item{positionsx}{Logical vector of length p indicating whether or not the
#' variable is a numerical covariate (when relevant)}
#' \item{imp.info}{List of arguments used for the imputation step and other
#' information (when relevant)}
#' \item{compress.imp.array}{Compressed array of imputed datasets (when relevant)}
#' \item{BF.method}{Function used to compute Bayes factors}
#' \item{prior.betas}{\code{prior.betas}}
#' \item{prior.models}{Vector with \code{prior.models} and \code{prior.models.dummies}
#' chosen. If there are no factors or \code{marginal.factors} is set to \code{FALSE},
#' it saves the only argument used, \code{prior.models}}
#' \item{family}{\code{\link[stats]{family}} object used to specify the error
#' distribution and link function in the model. If \code{BF.method == "gprior"}
#' additional function used to compute Bayes factors are added}
#' \item{marginal.factors}{Logical to indicate whether or not are marginalized factors'
#' probabilities such as García-Donato and Paulo (2022)}
#' \item{priorprobs}{Prior probabilities over the true model dimension}
#' \item{call}{The \code{call} to the function}
#'
#' @author Carolina Mulet and Gonzalo García-Donato
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[meantobvs]{missingbvs.glm}} for an exact computation
#' of the model posterior distribution (recommended when p<20).
#'
#' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
#'
#' Garcia-Donato, G. and Paulo, R. (2022)<DOI:10.1080/01621459.2021.1889565>
#' Variable Selection in the Presence of Factors: A Model Selection Perspective.
#' Journal of the American Statistical Association. 117. 1-27.
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
#' Li, Y. and Clyde, M. (2018)<DOI:10.1080/01621459.2018.1469992> Mixtures
#' of g-Priors in Generalized Linear Models. Journal of the American
#' Statistical Association. 113: 1275–1287.
#'
#' Clyde, M (2025) BAS: Bayesian Variable Selection and Model Averaging using
#' Bayesian Adaptive Sampling. R package version 2.0.2
#' <https://CRAN.R-project.org/package=BAS>.
#'
#' Schwarz, G. (1978) Estimating the dimension of a model. The Annals of
#' Statistics. 6(2): 461–464.
#'
#' Held L, Sabanés Bové D, Gravestock I (2015).<DOI:10.1214/14-STS510> Approximate
#' Bayesian Model Selection with the Deviance Statistic. Statistical Science. 30.
#'
#' van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice: Multivariate Imputation
#' by Chained Equations in R. Journal of Statistical Software. 45(3): 1–67.
#'
#' @examples
#' \donttest{
#'
#' #Indian Prime Diabetes Data
#'
#' #Build the list of competing models:
#' models.list = list(
#'   M0 = Outcome ~ 1,
#'   M1 = Outcome ~ Pregnancies,
#'   M2 = Outcome ~ Glucose,
#'   M3 = Outcome ~ Insulin,
#'   M4 = Outcome ~ Pregnancies + Glucose,
#'   M5 = Outcome ~ Pregnancies + Insulin,
#'   M6 = Outcome ~ Pregnancies + Glucose + Insulin
#' )
#'
#' #Default choices are: BIC approximation and Constant prior.
#' #Few imputations for simplicity, real analyses need more.
#' diabetes.meantobtest <- btest.glm(
#'   data = diabetes, models = models.list, family = binomial(),
#'   n.imp = 2, imp.seed = 1
#' )
#'
#' #Show the results:
#' diabetes.meantobtest
#' diabetes.meantobtest$PostProbi
#' }
#'
btest.glm <- function (data,
                       models,
                       family = binomial(link = "logit"),
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
                       imp.seed = runif(1,0,09011975),
                       weights = rep.int(1, nrow(data)),
                       offset = rep.int(0, nrow(data)),
                       control = glm.control()) {

  #N is the number of models:
  N <- length(models)

  #Check Btest given arguments
  btest.args <- checkbtestarguments(models, null.model, N)

  #Check family-link argument
  checkfamily.link(family)

  #Define list of arguments for posterior computation
  model.context <- list(
    relax.nest = btest.args$relax.nest,
    models = btest.args$models
  )

  Dev <- numeric(N) #deviances for each model
  Dim <- rep.int(0L, N)
  mt <- list() #list of terms for each model

  #goes through all competing models and saves some results
  covar.list <- list() #list that contains the names of the variables in each model
  compvars <- c() #name of original competing vars
  for (i in seq_len(N)) {
    f <- as.formula(btest.args$models[[i]])
    compvars <- c(compvars, attr(terms(f), "term.labels"))
    environment(f) <- environment()
    temp <- glm(formula = f,
                data = data,
                y = TRUE, x = TRUE,
                family = family,
                weights = weights,
                offset = offset,
                control = control)

    Dev[i] <- temp$deviance

    framei <- model.frame(f, data, na.action = NULL)
    Xi <- model.rankdef.matrix(framei)
    covar.list[[i]] <- dimnames(Xi)[[2]]
    Dim[i] <- length(covar.list[[i]])
    mt[[i]] <- temp$terms
  }
  ordered.Dev <- sort(Dev, index.return = TRUE, decreasing = TRUE)
  #Which one acts as null model:
  nullmodel.pos <- ordered.Dev$ix[1]

  model.context$nullmodel.pos <- nullmodel.pos
  model.context$covar.list <- covar.list

  #Check null model
  if (btest.args$relax.nest &&
      nullmodel.pos != btest.args$nullmodel.posuser) {
      stop("The given null model does not coincide with the one with\n",
           "the largest deviance (and it should).\n")
  }
  #change the string for the formula and specify models to compute BF
  null.model <- as.formula(btest.args$models[[nullmodel.pos]])
  environment(null.model) <- environment()

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

  #Evaluate the null model:
  glmnull <- glm(formula = null.model,
                 data,
                 y = TRUE, x = TRUE,
                 family = family,
                 weights = weights,
                 offset = offset,
                 control = control)
  #correct glm arguments
  family <- glmnull$family; weights <- glmnull$prior.weights; offset <- glmnull$offset
  #add a specification for the grouped binomial to check
  family$grouped <- ifelse(paste0(null.model[[2]])[1] == "cbind", TRUE, FALSE)

  #The response variable
  y <- glmnull$y; obsnotNA <- names(y) #without missings
  n <- length(y) #observations without missings on the response
  y <- as.numeric(y) #for the C code

  #check whether or not the family chosen is available for BF.method
  useBAS <- checkforfamily(family, BF.method)

  #Compute model prior
  lprior.models <- priormodels.btest(prior.models, N, Dim, priorprobs)

  #Check approx method and priors chosen and define the function to be used
  lBF <- checkforprior.betas.glm(
    BF.method, prior.betas, n, matrices$p, matrices$p0, y, glmnull, useBAS
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
      switch (as.character(useBAS),
              `TRUE` = {lBF.method <- function(model) lBF.av(
                  model, imputation.array = imputation$imputation.array,
                  lBF = lBF, p0 = matrices$p0, n.imp = n.imp
                )
              },
              `FALSE` = {lBF.method <- function(model) lBF.av.glm.fit(
                  model, imputation.array = imputation$imputation.array,
                  lBF = lBF, p0 = matrices$p0, n.imp = n.imp, y = y, glmnull = glmnull
                )
              }
      )
  } else {
    # If n.imp == 0 or 1, we do not need fitstart argument
    # if (anyNAvar) {
      lBF.method <- function(model) lBF(
        k = length(model),
        X = imputation$imputation.array[, c(seq_len(matrices$p0), model + matrices$p0), ],
        fitstart = NULL
      )
    # } else lBF.method <- function(model) lBF(
    #   k = length(model),
    #   X = cbind(matrices$X0, matrices$X.full[, model, drop = FALSE]),
    #   fitstart = NULL
    # )
  }
  #if no NAvars active we do not need fitstart
  lBFmodel <- function (model) lBF(k = length(model),
                                   X = cbind(matrices$X0, matrices$X.full[, model]),
                                   fitstart = NULL)

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

  #Evaluate glm of each model with missings using Rubin's rule
  modelspool <- list()
  for(j in competing.models){
    namesj <- which(matrices$namesxnotnull %in% covar.list[[j]])
    if (sum(NAvars[namesj]) > 0) {
      fit <- list()
      for (i in 1:n.imp) {
        if (mF) {#remove last dummy for each factor, first q0 vars are the fixed ones
          Xi <- imputation$imputation.array[,
            c(seq_len(matrices$p0), setdiff(namesj, indf) + matrices$p0), i
          ]
        } else Xi <- imputation$imputation.array[,
          c(seq_len(matrices$p0), namesj + matrices$p0), i
        ]

        z <- glm.fit(x = Xi, y = y, family = family,
                     weights = weights, offset = offset, control = control)
        z$terms <- mt[[j]]; class(z) <- "glm"; fit[[i]] <- z
      }
      modelspool[[j]] <- mice::pool(fit)
      modelspool[[j]]$call <- NULL #otherwise, Rstudio returns a warning trying to read modelspool[[j]]$call
    } else {
      if (mF) {#remove last dummy for each factor, first q0 vars are the fixed ones
        Xi <- imputation$imputation.array[,
                                          c(seq_len(matrices$p0), setdiff(namesj, indf) + matrices$p0), i
        ]
      } else Xi <- imputation$imputation.array[,
                                               c(seq_len(matrices$p0), namesj + matrices$p0), i
      ]

      # fj <- as.formula(btest.args$models[[j]]); environment(fj) <- environment()
      modelspool[[j]] <- glm.fit(x = Xi, y = y, family = family,
                                 weights = weights, offset = offset, control = control)
    }
  }
  modelspool[[nullmodel.pos]] <- glm.fit(x = matrices$X0, y = y, family = family,
                                         weights = weights, offset = offset, control = control)
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

    #save the imputed datasets for sensitivity analysis
    raw.imp.array <- serialize(imputation$imputation.array, NULL)
    result$compress.imp.array <- memCompress(raw.imp.array, type = "xz")
  }

  #chosen methods
  result$BF.method <- BF.method #method used for BF computation
  if (is.null(prior.betas) & BF.method %in% c("gprior", "TBF")) prior.betas <- "gZellner"
  if (BF.method == "BIC") prior.betas <- NULL #does not apply with BIC
  result$prior.betas <- prior.betas

  if (mF) {
    result$prior.models <- c(prior.models, prior.models.dummies)
  } else result$prior.models <- prior.models

  #update family functions used
  if (BF.method == "gprior" & !useBAS) {
    family <- addfamilyfun(glmnull$family)
  } else family <- glmnull$family

  result$marginal.factors <- marginal.factors #whether or not factors are marginalized

  result$priorprobs <- exp(posterior$lPriorModels)
  result$call <- match.call()

  class(result) <- "meantobtest"

  return(result)
}
