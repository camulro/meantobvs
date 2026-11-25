#' Logarithm of the MC approximation of the GD25 Bayes factor in linear models
#'
#' Computes the logarithm of the Bayes factor in linear regression models when
#' covariates are random and normally distributed (García-Donato et al., 2025).
#'
#' The Bayes factor approximation is via a MonteCarlo scheme given the imputed datasets
#' and variance-covariance matrices simuled in a \code{meantobvs.imputation} object.
#' The g'BF is computed with the auxiliary function \code{lBF}.
#'
#' @param model Vector of indexes in \{1, 2, ..., p\} denoting the active variables
#' for a given model with p competing covariates.
#' @param imputation.list Object of class \code{MC.imputation} with the
#' following elements: \code{rX.imput}Array of dimension \code{n}xpx\code{nMC}
#' containing the imputed datasets; \code{rSigma}Array of dimension
#' pxpx\code{nMC} containing the corresponding covariance matrices
#' @param lBF Auxiliary function with needed fixed parameters to compute
#' the g'-Bayes factor over \code{imputation.list}.
#' @param n Number of observations.
#' @param nMC Number of samples used to approximate, by MonteCarlo, the integral
#' defining the Bayes factor.
#'
#' @return \code{lBF.miss} returns, in logarithmic scale, the MonteCarlo
#' approximation of the integral defining the Bayes factor for a given \code{model}
#' following García-Donato et al (2025).
#'
#' @author María Eugenia Castellanos and Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[meantobvs]{MC.imputation}} for computing the
#' \code{meantobvs.imputation} object used in the MonteCarlo approximation and
#' \code{\link[meantobvs]{BF.GD25}} for the Bayes factor for each MC step.
#' Use \code{\link[meantobvs]{bvs.lm}} for an exact computation
#' of the model posterior distribution in the VS problem (recommended when p<20).
#'
#' @examplesIf interactive()
#' #Daily air quality measurements in New York
#' data("airquality")
#' imp2 <- MC.imputation(X = airquality[,c("Ozone", "Wind", "Temp")], nMC = 2)
#'
#' lmnull <- lm(Solar.R ~ 1, data = airquality, y = T)
#' imp2$rX.imput <- imp2$rX.imput[-lmnull$na.action,,]
#' BF.fun <- function(X.center, Sigma11, k) meantobvs:::BF.GD25(X.center, Sigma11,
#'   y = lmnull$y, SS0 = crossprod(lmnull$residuals))
#' lBF <- meantobvs:::lBF.miss(1:3, imp2, lBF = BF.fun)
#'
#' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
#'
#' @keywords internal
lBF.miss <- function(model, imputation.list, lBF,
                     n = dim(imputation.list$rX.imput)[1], nMC = dim(imputation.list$rSigma)[3]) {
  k <- length(model)

  imputation.list.model <- list()
  imputation.list.model$rX.imput <- array(imputation.list$rX.imput[,model,],
                                          dim = c(n, k, dim(imputation.list$rSigma)[3]))
  imputation.list.model$rSigma <- array(imputation.list$rSigma[model, model,],
                                        dim = c(k, k, dim(imputation.list$rSigma)[3]))

  lBF.our <- numeric(nMC)
  for(s in 1:nMC) {
    #posterior dist. with Jeffreys independent prior
    Sigma11 <- imputation.list.model$rSigma[,,s]
    # mu <- imputation.list$rmu[,s]
    X.center <- imputation.list.model$rX.imput[,,s] #centered at imputation step

    lBF.our[s] <- lBF(X.center, Sigma11, k) #log(BFmodel0) for the sth imputation
  }

  lBF.miss <- logsumexp.stable(lBF.our) - log(nMC) #log(mean(exp(lBF.our)))
  return(lBF.miss)
}

#' Logarithm of the g'-Bayes factor for an iteration in the MC step
#'
#' Computes the logarithm of the g'-Bayes factor in the linear regression model
#' with complete data, for normally distributed covariates. It employs the imputed
#' data and covariance matrices of an \code{MC.imputation} object.
#'
#' @param X.center Matrix of dimension \code{n}x\code{k} containing the imputed data.
#' @param Sigma11 Matrix of dimension \code{k}x\code{k} containing covariance matrices.
#' @param y Response variable in the linear model.
#' @param SS0 Sum of squared error of the null model.
#' @param n Number of observations.
#' @param k Number of model-specific coefficients.
#'
#' @return \code{BF.GD25} returns, in logarithmic scale, the g'-Bayes factor
#' of García-Donato et al (2025) for imputed data \code{X.center} and covariance
#' matrix \code{Sigma11}.
#'
#' @author María Eugenia Castellanos
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[meantobvs]{MC.imputation}} for computing the
#' \code{meantobvs.imputation} object containing the \code{X.center} and \code{Sigma11}
#' matrices. Use \code{\link[meantobvs]{bvs.lm}} for an exact computation
#' of the model posterior distribution in the VS problem (recommended when p<20).
#'
#' @examplesIf interactive()
#' #Daily air quality measurements in New York
#' data("airquality")
#' imp1 <- MC.imputation(X = airquality[,c("Ozone", "Wind", "Temp")], nMC = 1)
#'
#' lmnull <- lm(Solar.R ~ 1, data = airquality, y = T)
#' lBF.imp1 <- meantobvs:::BF.GD25(imp1$rX.imput[-lmnull$na.action,,], imp1$rSigma[,,1],
#'   y = lmnull$y, SS0 = crossprod(lmnull$residuals))
#'
#' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
#'
#' @keywords internal
BF.GD25 <- function(X.center, Sigma11, y, SS0, n = length(y), k = ncol(X.center)) {
  tX.center.X.center <- crossprod(X.center)

  lBFi0 <- -.5*(n-1)*log(1-t(y) %*% X.center %*% solve((tX.center.X.center + Sigma11)) %*% t(X.center) %*% y/SS0) -
                .5*determinant(tX.center.X.center %*% solve(Sigma11) + diag(rep(1,k)), log=T)$modulus[1]

  return(lBFi0)
}

#' Logarithm of the Average Bayes factor for any regression model with missing data
#'
#' Computes the logarithm of the Average Bayes factor (AvBF) when missingness occurs,
#' for any type of covariates and regression model.
#'
#' The AvBF is computed by averaging over the \code{n.imp} imputed datasets given by a
#' \code{meantobvs.imputation} object, the \code{n.imp} data-driven BF for each given
#' imputation. These BFs are computed with the auxiliary function \code{lBF}.
#'
#' @param model Vector of indexes in \{1, 2, ..., p\} denoting the active variables
#' for a given model with p competing covariates.
#' @param imputation.array Array of dimension nx(p+\code{p0})x\code{n.imp}
#' containing the imputed datasets, where n is the number of observations.
#' @param lBF Auxiliary function with some parameters fixed to compute
#' the data-driven Bayes factor over each imputed datased in \code{imputation.list}.
#' @param p0 Number of fixed covariates (including the intercept term).
#' @param n.imp Number of imputed datasets.
#'
#' @return \code{lBF.av} returns, in logarithmic scale, the Average Bayes
#' factor over the \code{n.imp} imputed datasets in \code{imputation.array}.
#'
#' @author María Eugenia Castellanos and Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[meantobvs]{mice.imputation}} for computing the
#' \code{meantobvs.imputation} object used in the average. Use
#' \code{\link[meantobvs]{bvs.lm}} with linear models or
#' \code{\link[meantobvs]{bvs.glm}} with generalized linear models for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examplesIf interactive()
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#'
#' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#' imp2 <- mice.imputation(model.frame(f, dataS97, na.action = NULL),
#'   n.imp = 2, seed = 1)$imputation.array
#'
#' lmnull <- lm(gr56092 ~ 1, data = dataS97, y = T)
#' BF.fun <- function(X, k) meantobvs:::BF.BIC.lm(y = lmnull$y,
#'   X, SS0 = crossprod(lmnull$residuals))
#' lBF <- meantobvs:::lBF.av(1:3, imp2[-lmnull$na.action,,], lBF = BF.fun)
#'
#' @references García-Donato, G., Castellanos, M.E., Cabras, S., Quirós, A.
#' and Forte, A. (2025) Model Uncertainty and Missing Data: An Objective Bayesian
#' Perspective (with Discussion). Bayesian Analysis. 20: 1677–1778.
#'
#' van Buuren, S. and Groothuis-Oudshoorn, K. (2011) mice: Multivariate Imputation
#' by Chained Equations in R. Journal of Statistical Software. 45(3): 1–67.
#'
#' @keywords internal
lBF.av <- function(model, imputation.array, lBF, p0 = 1, n.imp = dim(imputation.array)[3]) {
  k <- length(model)

  lBF.aux <- numeric(n.imp)
  X1.array <- imputation.array[,c(1:p0, model+p0),] #first p0 columns are fixed
  for(s in 1:n.imp) {
    lBF.aux[s] <- lBF(k = k, X = X1.array[,,s]) #BF function defined previously
  }

  lBF.av <- logsumexp.stable(lBF.aux) - log(n.imp) #log(mean(exp(lBF.aux)))
  return(lBF.av)
}

#' For BF.method.glm.fit functions. Analogue to lBF.av but it previously computes
#' glm estimated coefficients for model on the first imputation in order to give
#' this value as a starting point of the IRLS algorithm -and save computation time-
#' for next iterations.
#'
#' @keywords internal
lBF.av.glm.fit <- function(model, imputation.array, lBF, p0 = 1,
                           n.imp = dim(imputation.array)[3], y, glmnull) {
  k <- length(model)

  #use first imputation to estimate coefficients (to accelerate convergence on the rest)
  X <- imputation.array[, c(1:p0, model+p0), 1]

  control <- glm2fastglm.control(glmnull$control)
  fit <- fastglm::fastglmPure(y = y, x = X,
                              family = glmnull$family,
                              weights = glmnull$prior.weights,
                              offset = glmnull$offset,
                              tol = control$tol,
                              maxit = control$maxit,
                              method = control$method)

  lBF.aux <- numeric(n.imp)
  X1.array <- imputation.array[,c(1:p0, model+p0),] #first p0 columns are fixed
  for(s in 1:n.imp) {
    #use estimated coefficients on first dataset to achieve fast convergence:
    #if datasets slightly change, coefficients are going to take close values
    lBF.aux[s] <- lBF(k = k, X = X1.array[,,s], fitstart = fit$coefficients) #BF function defined previously
  }

  lBF.av <- logsumexp.stable(lBF.aux) - log(n.imp) #log(mean(exp(lBF.aux)))
  return(lBF.av)
}

#' Logarithm of the BIC approximation of the Bayes factor in lm
#'
#' Computes the logarithm of the BIC approximation of Bayes factors for
#' complete data in linear models.
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param SS0 Sum of squared error of the null model considered.
#' @param n Number of observations.
#' @param k Number of model-specific coefficients.
#' @param p0 Number of fixed covariates (including the intercept).
#'
#' @return \code{BF.BIC.lm} returns, in logarithmic scale, the BIC
#' approximation of the Bayes factor for a given model through \code{X}.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[meantobvs]{lBF.av}} to compute the average Bayes
#' factor for missing data. Use \code{\link[meantobvs]{bvs.lm}} for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examplesIf interactive()
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#'
#' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#' imp1 <- mice.imputation(model.frame(f, dataS97, na.action = NULL),
#'   n.imp = 1, seed = 1)$imputation.array
#'
#' lmnull <- lm(gr56092 ~ 1, data = dataS97, y = T)
#' lBF <- meantobvs:::BF.BIC.lm(y = lmnull$y, X = imp1[-lmnull$na.action,,],
#'   SS0 = crossprod(lmnull$residuals))
#'
#' @references Schwarz, G. (1978) Estimating the dimension of a model. The
#' Annals of Statistics. 6(2): 461–464.
#'
#' @keywords internal
BF.BIC.lm <- function(y, X, SS0, n = length(y), k = ncol(X)-p0, p0 = 1L) {

  SSE.model <- as.numeric(crossprod(.lm.fit(y = y, x = X)$residuals))

  # BFi0 <- (SS0/SSE.model)^(n/2) / n^(k/2)
  lBFi0 <- n/2 * log(SS0/SSE.model) - k/2 * log(n) #exp((BIC0 - BICi)/2)
  return(lBFi0)
}

#' Logarithm of the test-based Bayes factor (TBF) in lm
#'
#' Computes the logarithm of the TBF approximation of Bayes factors for
#' complete data in linear models.
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param SS0 Sum of squared error of the null model considered.
#' @param lTBF.method Function to compute log-TBF, with the corresponding fixed parameters.
#' @param n Number of observations.
#' @param k Number of model-specific coefficients.
#' @param p0 Number of fixed covariates (including the intercept).
#'
#' @return \code{BF.TBF.lm} returns, in logarithmic scale, the TBF
#' approximation of the Bayes factor in linear models for a given model
#' through \code{X}.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[meantobvs]{lBF.av}} to compute the average Bayes
#' factor for missing data. Use \code{\link[meantobvs]{bvs.lm}} for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examplesIf interactive()
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#'
#' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#' imp1 <- mice.imputation(model.frame(f, dataS97, na.action = NULL),
#'   n.imp = 1, seed = 1)$imputation.array
#'
#' lmnull <- lm(gr56092 ~ 1, data = dataS97, y = T)
#' lTBF.method <- lTBF.function("gZellner", length(lmnull$y), ncol(imp1), 0)
#' lBF <- meantobvs:::BF.TBF.lm(y = lmnull$y, X = imp1[-lmnull$na.action,,],
#'   SS0 = crossprod(lmnull$residuals), lTBF.method = lTBF.method)
#'
#' @references Held, L., Sabanés Bové, D. and Gravestock, I.
#' (2015)<DOI:10.1214/14-STS510> Approximate Bayesian Model Selection with the
#' Deviance Statistic. Statistical Science, 30(2): 242–257.
#'
#' @keywords internal
BF.TBF.lm <- function(y, X, SS0, lTBF.method, n = length(y), k = ncol(X)-p0, p0 = 1L) {

  R2j <- 1 - as.numeric(crossprod(.lm.fit(y = y, x = X)$residuals))/SS0
  minuszj <- n * log(1 - R2j) #for LM

  # lBFi0 <- -k/2 * log(g + 1) + g/(g+1) * zj/2 #fixed g
  lBFi0 <- lTBF.method(k = k, dev = minuszj)
  return(lBFi0)
}

#' Logarithm of the g-prior Bayes factor in lm
#'
#' Computes the logarithm of the Bayes factors derived from a given g-prior
#' for complete data in linear models.
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param SS0 Sum of squared error of the null model considered.
#' @param prior.betas Prior distribution for model-specific coefficients in the
#' \pkg{BayesVarSel} codification. Options include "gBF", "RobustBF", "LiangBF",
#' "ZSBF", intrinsicBF" and "geointrinsicBF". See
#' \code{\link[meantobvs]{bvs.lm}} for more details.
#' @param n Integer, number of observations.
#' @param k Integer, number of model-specific coefficients.
#' @param p0 Integer, number of fixed covariates (including the intercept).
#'
#' @return \code{BF.gprior.lm} returns, in logarithmic scale, the exact
#' value of the Bayes factor derived from assigning a chosen g-prior by
#' \code{prior.betas} in linear models for a given model through \code{X}.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[meantobvs]{lBF.av}} to compute the average Bayes
#' factor for missing data. Use \code{\link[meantobvs]{bvs.lm}} for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examplesIf interactive()
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#'
#' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#' imp1 <- mice.imputation(model.frame(f, dataS97, na.action = NULL),
#'   n.imp = 1, seed = 1)$imputation.array
#'
#' lmnull <- lm(gr56092 ~ 1, data = dataS97, y = T)
#' lBF <- meantobvs:::BF.gprior.lm(y = lmnull$y, X = imp1[-lmnull$na.action,,],
#'   SS0 = crossprod(lmnull$residuals))
#'
#' @references García-Donato, G. and Forte, A. (2018) Bayesian Testing,
#' Variable Selection and Model Averaging in Linear Models using R with
#' BayesVarSel. The R Journal. 10: 329.
#'
#' Bayarri, M.J., Berger, J.O., Forte, A. and Garcia-Donato, G.
#' (2012)<DOI:10.1214/12-aos1013> Criteria for Bayesian Model choice with
#' Application to Variable Selection. The Annals of Statistics. 40: 1550-1557.
#'
#' Berger, J., Garcıa-Donato, G., Moreno, E., and Pericchi, L. (2022).
#' The intrinsic hyper-g prior for normal linear models. in preparation.
#'
#' Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger,J.O.
#' (2008)<DOI:10.1198/016214507000001337> Mixtures of g-priors for Bayesian
#' Variable Selection. Journal of the American Statistical Association.
#' 103:410-423
#'
#' Moreno, E., Giron, J. and Casella, G. (2015) Posterior model consistency
#' in variable selection as the model dimension grows. Statistical Science. 30: 228-241.
#'
#' Zellner, A. and Siow, A. (1980)<DOI:10.1007/bf02888369> Posterior Odds Ratio
#' for Selected Regression Hypotheses. In Bayesian Statistics 1 (J.M. Bernardo,
#' M. H. DeGroot, D. V. Lindley and A. F. M. Smith, eds.) 585-603. Valencia:
#' University Press.
#'
#' Zellner, A. and Siow, A. (1984). Basic Issues in Econometrics. Chicago:
#' University of Chicago Press.
#'
#' Zellner, A. (1986)<DOI:10.2307/2233941> On Assessing Prior Distributions and
#' Bayesian Regression Analysis with g-prior Distributions. In Bayesian
#' Inference and Decision techniques: Essays in Honor of Bruno de Finetti (A.
#' Zellner, ed.) 389-399. Edward Elgar Publishing Limited.
#'
#' @keywords internal
BF.gprior.lm <- function(y, X, SS0, prior.betas = "gBF",
                         n = length(y), k = as.integer(ncol(X)-p0), p0 = 1L) {

  SSE.model <- as.numeric(crossprod(.lm.fit(y = y, x = X)$residuals))

  BFi0 <- .C(prior.betas, n, k + p0, p0, as.double(SSE.model/SS0), 0.0,
             PACKAGE = "BayesVarSel")[5][[1]]
  return(log(BFi0))
}
#' Logarithm of the FLS Bayes factor in lm
#'
#' Computes the logarithm of the Bayes factors derived from the
#' Fernandez, Ley and Steel (2001) Benchmark prior for complete data in linear models.
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param SS0 Sum of squared error of the null model considered.
#' @param dmax Integer, maximum model dimension of potential regressors.
#' @param n Integer, number of observations.
#' @param k Integer, number of model-specific coefficients.
#' @param p0 Integer, number of fixed covariates (including the intercept).
#'
#' @return \code{BF.FLS.lm} returns, in logarithmic scale, the exact
#' value of the Bayes factor derived from assigning the FLS Benchmark g-prior
#' in linear models for a given model through \code{X}. See
#' \code{\link[meantobvs]{bvs.lm}} for more details.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[meantobvs]{lBF.av}} to compute the average Bayes
#' factor for missing data. Use \code{\link[meantobvs]{bvs.lm}} for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examplesIf interactive()
#' #Cross-Country Growth, from Fernández, Ley and Steel (2001)
#'
#' f <- gr56092 ~ 1 + lifee060 + gdpsh60l + p60
#' imp1 <- mice.imputation(model.frame(f, dataS97, na.action = NULL),
#'   n.imp = 1, seed = 1)$imputation.array
#'
#' lmnull <- lm(gr56092 ~ 1, data = dataS97, y = T)
#' lBF <- meantobvs:::BF.FLS.lm(y = lmnull$y, X = imp1[-lmnull$na.action,,],
#'   SS0 = crossprod(lmnull$residuals), dmax = ncol(imp1))
#'
#' @references García-Donato, G. and Forte, A. (2018) Bayesian Testing,
#' Variable Selection and Model Averaging in Linear Models using R with
#' BayesVarSel. The R Journal. 10: 329.
#'
#' Fernandez, C., Ley, E. and Steel, M.F.J.
#' (2001)<DOI:10.1016/s0304-4076(00)00076-2> Benchmark priors for Bayesian
#' model averaging. Journal of Econometrics, 100, 381-427.
#'
#' @keywords internal
BF.FLS.lm <- function(y, X, SS0, dmax,
                      n = length(y), k = as.integer(ncol(X)-p0), p0 = 1L) {

  SSE.model <- crossprod(.lm.fit(y = y, x = X)$residuals)

  BFi0 <- .C("flsBF", dmax - p0, n, k + p0, p0, as.double(SSE.model/SS0), 0.0,
             PACKAGE = "BayesVarSel")[6][[1]]
  return(log(BFi0))
}

#' Logarithm of the BIC approximation of the Bayes factor in glm
#'
#' Computes the logarithm of the BIC approximation of Bayes factors for
#' complete data in generalized linear models.
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param family String, function or the call to a family function among
#' \code{\link[stats]{family}} to specify the error distribution and link
#' function to be used in the model.
#' @param devnull Deviance of the null model considered.
#' @param n Number of observations.
#' @param k Number of model-specific coefficients.
#' @param weights NULL or numeric vector of the same length as \code{y} to
#' specify the weights to be used in the glm fitting process.
#' @param offset NULL or a numeric vector of the same length as \code{y} to
#' specify an a priori known component included in the glm fitting process.
#' @param control List of parameters for controlling the
#' \code{\link[fastglm]{fastglm}} fitting process.
#' @param fitstart Optional starting values for the coefficients in the linear
#' predictor. By default, it is \code{NULL}.
#'
#' @return \code{BF.BIC.glm.fit} returns, in logarithmic scale, the BIC
#' approximation of the Bayes factor in generalized linear models for a given
#' model through \code{X}.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[meantobvs]{lBF.av.glm.fit}} to compute the average Bayes
#' factor for missing data. Use \code{\link[meantobvs]{bvs.glm}} for
#' an exact computation of the model posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examplesIf interactive()
#' #Indian Prime Diabetes Data
#'
#' f <- Outcome ~ Pregnancies + Glucose + BloodPressure + SkinThickness + Insulin
#' imp1 <- mice.imputation(model.frame(f, diabetes, na.action = NULL),
#' n.imp = 1, seed = 1)$imputation.array[,,1]
#'
#' glmnull <- glm(Outcome ~ 1, data = diabetes, family = binomial(), y = TRUE)
#' lBF <- meantobvs:::BF.BIC.glm.fit(y = glmnull$y, X = imp1,
#'   family = binomial(), devnull = glmnull$deviance)
#'
#' @references Schwarz, G. (1978) Estimating the dimension of a model. The
#' Annals of Statistics. 6(2): 461–464.
#'
#' @keywords internal
BF.BIC.glm.fit <- function(y, X, family = binomial(link = "logit"),
                           devnull,
                           n = length(y), k = ncol(X)-1,
                           weights = rep(1, n),
                           offset = rep(0, n),
                           control = glm2fastglm.control(),
                           fitstart = NULL) { #initial values for coeffs

  fit1 <- fastglm::fastglmPure(y = y, x = X,
                               family = family,
                               start = fitstart,
                               weights = weights,
                               offset = offset,
                               tol = control$tol,
                               maxit = control$maxit,
                               method = control$method)

  #Check for fitting convergence
  if (!fit1$converged) warning("The fitting algorithm did not converge.\n",
                               immediate. = TRUE)

  lBFi0 <- (devnull - fit1$deviance - k * log(n))/2
  return(lBFi0)
}

#' Logarithm of the test-based Bayes factor (TBF) in glm
#'
#' Computes the logarithm of the TBF approximation of Bayes factors for
#' complete data in generalized linear models.
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param family String, function or the call to a family function among
#' \code{\link[stats]{family}} to specify the error distribution and link
#' function to be used in the model.
#' @param n Number of observations.
#' @param k Number of model-specific coefficients.
#' @param lTBF.method Function to compute log-TBF, with the corresponding fixed parameters.
#' @param weights NULL or numeric vector of the same length as \code{y} to
#' specify the weights to be used in the glm fitting process.
#' @param offset NULL or a numeric vector of the same length as \code{y} to
#' specify an a priori known component included in the glm fitting process.
#' @param control List of parameters for controlling the
#' \code{\link[fastglm]{fastglm}} fitting process.
#' @param fitstart Optional starting values for the coefficients in the linear
#' predictor. By default, it is \code{NULL}.
#'
#' @return \code{BF.TBF.glm.fit} returns, in logarithmic scale, the TBF
#' approximation of the Bayes factor in generalized linear models for a given
#' model through \code{X}.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[meantobvs]{lBF.av.glm.fit}} to compute the average Bayes
#' factor for missing data. Use \code{\link[meantobvs]{bvs.glm}} for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examplesIf interactive()
#' #Indian Prime Diabetes Data
#'
#' f <- Outcome ~ Pregnancies + Glucose + BloodPressure + SkinThickness + Insulin
#' imp1 <- mice.imputation(model.frame(f, diabetes, na.action = NULL),
#'   n.imp = 1, seed = 1)$imputation.array[,,1]
#'
#' glmnull <- glm(Outcome ~ 1, data = diabetes, family = binomial(), y = TRUE)
#' lTBF.method <- lTBF.glm.function("gZellner",
#'   length(glmnull$y), ncol(imp1), glmnull$deviance)
#' lBF <- meantobvs:::BF.TBF.glm.fit(y = glmnull$y, X = imp1,
#'   family = binomial(), lTBF.method = lTBF.method)
#'
#' @references Held, L., Sabanés Bové, D. and Gravestock, I.
#' (2015)<DOI:10.1214/14-STS510> Approximate Bayesian Model Selection with the
#' Deviance Statistic. Statistical Science, 30(2): 242–257.
#'
#' @keywords internal
BF.TBF.glm.fit <- function(y, X, family = binomial(link = "logit"),
                           n = length(y), k = ncol(X) - 1,
                           lTBF.method,
                           weights = rep(1, length(y)),
                           offset = rep(0, length(y)),
                           control = glm2fastglm.control(),
                           fitstart = NULL) { #initial values for coeffs

  fit1 <- fastglm::fastglmPure(y = y, x = X,
                               family = family,
                               start = fitstart,
                               weights = weights,
                               offset = offset,
                               tol = control$tol,
                               maxit = control$maxit,
                               method = control$method)

  #Check for fitting convergence
  if (!fit1$converged) warning("The fitting algorithm did not converge.\n",
                               immediate. = TRUE)

  lBFi0 <- lTBF.method(k = k, dev = fit1$deviance)
  return(lBFi0)
}

#' @keywords internal
lTBF.gfixed <- function (g, k, dev, devnull) {
  #general formula for the TBF with g fixed

  -k/2 * log(g + 1) + g/(g+1) * (devnull - dev)/2
}

#' @keywords internal
lTBF.grandom <- function (a, b, k, dev, devnull) {
  #general formula for the TBF with g ~ IncIG(a,b)

  # b^a * pgamma(b + dev/2, a + k/2) * gamma(a + k/2) *
  #   ((b + dev/2)^(a + k/2) * pgamma(b, a) * gamma(a))^(-1) * exp(dev/2)

  a*log(b) - log(pgamma(b, a)) - lgamma(a) - (a + k/2) * log(b + (devnull - dev)/2) +
  log(pgamma(b + (devnull - dev)/2, a + k/2)) + lgamma(a + k/2) + (devnull - dev)/2
}

#' @keywords internal
lTBF.hyperg <- function (k, dev, devnull) {
  #particular formula for the TBF Liang et al hyper-g version: a=1, b=0

  1 - (1 + k/2) * log((devnull - dev)/2) +
    log(pgamma((devnull - dev)/2, 1 + k/2)) + lgamma(1 + k/2) + (devnull - dev)/2
}

#' Logarithm of the g-prior (g fixed) Bayes factor in glm
#'
#' Computes the logarithm of the Bayes factors derived from a given g-prior
#' with g fixed, for complete data in generalized linear models.
#'
#' Computes the approximated expression for Bayes factors under g-priors
#' with g fixed, in logarithmic scale:
#'  lBF_10 = z1/2 + 1/2 log(Jaa0/Jaa1) - k/2 log(1+g) - Q1/(2 * (1+g)).
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param family \code{\link[stats]{family}} to specify the error distribution and link
#' function to be used in the model. Must contain derivatives needed to compute
#' observed information.
#' @param glm0 glm object derived from fitting the model with just the intercept.
#' @param g Fixed value g for g-prior.
#' @param n Number of observations.
#' @param k Number of model coefficients (without the intercept).
#' @param weights NULL or numeric vector of the same length as \code{y} to
#' specify the weights to be used in the glm fitting process.
#' @param offset NULL or a numeric vector of the same length as \code{y} to
#' specify an a priori known component included in the glm fitting process.
#' @param control List of parameters for controlling the
#' \code{\link[fastglm]{fastglm}} fitting process.
#' @param fitstart Optional starting values for the coefficients in the linear
#' predictor. By default, it is \code{NULL}.
#'
#' @return \code{BF.gprior.glm} returns, in logarithmic scale, the exact
#' value of the Bayes factor derived from assigning a chosen g-prior by
#' \code{prior.betas} in generalized linear models for a given model through
#' \code{X}.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[meantobvs]{lBF.av.glm.fit}} to compute the average Bayes
#' factor for missing data. Use \code{\link[meantobvs]{bvs.glm}} for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examplesIf interactive()
#' #Indian Prime Diabetes Data
#'
#' f <- Outcome ~ Pregnancies + Glucose + BloodPressure + SkinThickness + Insulin
#' imp1 <- mice.imputation(model.frame(f, diabetes, na.action = NULL),
#'   n.imp = 1, seed = 1)$imputation.array[,,1]
#'
#' glmnull <- glm(Outcome ~ 1, data = diabetes, family = binomial(), x = T, y = T)
#' lBF <- meantobvs:::BF.gprior.glm.fit(y = glmnull$y, X = imp1, glmnull)
#'
#' @references
#' Li, Y. and Clyde, M. (2018)<DOI:10.1080/01621459.2018.1469992> Mixtures of
#' g-priors in Generalized Linear Models. Journal of the American Statistical
#' Association. 113: 1828-1845
#'
#' @keywords internal
BF.gprior.glm.fit <- function (y, X, glm0,
                               family = addfamilyfun(glm0$family),
                               g = n, #default corresponds to gZellner
                               n = length(y), k = ncol(X)-1,#p0=1 always
                               weights = rep(1, n),
                               offset = rep(0, n),
                               control = glm2fastglm.control(),
                               fitstart = NULL) { #initial values for coeffs

  fit1 <- fastglm::fastglmPure(y = y, x = X,
                               family = glm0$family,
                               start = fitstart,
                               weights = weights,
                               offset = offset,
                               tol = control$tol,
                               maxit = control$maxit,
                               method = control$method)

  #Check for fitting convergence
  if (!fit1$converged) warning("The fitting algorithm did not converge.\n",
                               immediate. = TRUE)

  #J = X^tWX
  # w <- fit1$weights #W = diag(w_i)

  #without intercept:
  Xbeta <- X[, -1, drop = FALSE]
  beta <- fit1$coefficients[-1]

  eta <- drop(Xbeta %*% beta)

  # Jaa <- sum(w) #information associated to intercept
  D <- Fobs.info(fit1, family)
  Jaa <- sum(D)

  # Wald statistic
  # Q <- sum(w * eta^2) - sum(w * eta)^2 / Jaa
  Q <- sum(D * eta^2) - sum(D * eta)^2 / Jaa

  #Get information matrix for null model
  # Jaa0 <- sum(glm0$weights)
  D0 <- Fobs.info(glm0, family)
  Jaa0 <- sum(D0)

  #change in deviance
  z <- glm0$deviance - fit1$deviance

  #Bayes factor derived from approximate marginal likelihood (Li and Clyde, 2018)
  lBFi0 <- z/2 - log(Jaa / Jaa0)/2 - k/2 * log1p(g) - Q/(2 * (1 + g))
  return(lBFi0)
}

#' Logarithm of the hyper g-prior Bayes factor in glm
#'
#' Computes the logarithm of the Bayes factors derived from a given hyper g-prior,
#'  for complete data in generalized linear models.
#'
#' Computes the approximated expression for Bayes factors under hyper g-priors,
#' in logarithmic scale:
#'  lBF_10 = z1/2 - Q1/(v) + 1/2  log(Jaa0/Jaa1) - k/2 * v +
#     log B((a + k)/2, b/2) + log phi1(b/2, r, (a + b + k)/2, (s + Q1)/(2*v), 1 - ka) -
#     log B(a/2, b/2) - log phi1(b/2, r, (a + b)/2, s/(2*v), 1 - ka)
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param family \code{\link[stats]{family}} to specify the error distribution and link
#' function to be used in the model. Must contain derivatives needed to compute
#' observed information.
#' @param glm0 glm object derived from fitting the model with just the intercept.
#' @param prior.betas.args List of arguments corresponding to a specific hyper
#' g-prior. It must contain a, b, r, s, v and ka. By default corresponds to
#' Liangetal prior.
#' @param n Number of observations.
#' @param k Number of model coefficients (without the intercept).
#' @param weights NULL or numeric vector of the same length as \code{y} to
#' specify the weights to be used in the glm fitting process.
#' @param offset NULL or a numeric vector of the same length as \code{y} to
#' specify an a priori known component included in the glm fitting process.
#' @param fitstart Optional starting values for the coefficients in the linear
#' predictor. By default, it is \code{NULL}.
#' @param control List of parameters for controlling the
#' \code{\link[fastglm]{fastglm}} fitting process.
#'
#' @return \code{BF.gprior.glm} returns, in logarithmic scale, the exact
#' value of the Bayes factor derived from assigning a chosen g-prior by
#' \code{prior.betas} in generalized linear models for a given model through
#' \code{X}.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[meantobvs]{lBF.av.glm.fit}} to compute the average Bayes
#' factor for missing data. Use \code{\link[meantobvs]{bvs.glm}} for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examplesIf interactive()
#' #Indian Prime Diabetes Data
#'
#' f <- Outcome ~ Pregnancies + Glucose + BloodPressure + SkinThickness + Insulin
#' imp1 <- mice.imputation(model.frame(f, diabetes, na.action = NULL),
#'   n.imp = 1, seed = 1)$imputation.array[,,1]
#'
#' glm0 <- glm(Outcome ~ 1, data = diabetes, family = binomial(), x = T, y = T)
#' lBF <- meantobvs:::BF.hypergprior.glm.fit(y = glm0$y, X = imp1, glm0)
#'
#' @references
#' Li, Y. and Clyde, M. (2018)<DOI:10.1080/01621459.2018.1469992> Mixtures of
#' g-priors in Generalized Linear Models. Journal of the American Statistical
#' Association. 113: 1828-1845
#'
#' @keywords internal
BF.hypergprior.glm.fit <- function (y, X, glm0,
                                    family = addfamilyfun(glm0$family),
                                    prior.betas.args =
                                      prior.betas.gprior.glm.fit("Liangetal", n, k),
                                    #default coresponds to liangetal prior
                                    n = length(y), k = ncol(X)-1, #p0=1 always
                                    weights = rep(1, n),
                                    offset = rep(0, n),
                                    control = glm2fastglm.control(),
                                    fitstart = NULL) { #initial values for coeffs

  #Get prior arguments
  list2env(prior.betas.args, envir = environment())

  #update arguments that depend on k, if needed
  if (is.function(v)) v <- v(k)
  if (is.function(ka)) ka <- ka(k)

  fit1 <- fastglm::fastglmPure(y = y, x = X,
                               family = glm0$family,
                               start = fitstart,
                               weights = weights,
                               offset = offset,
                               tol = control$tol,
                               maxit = control$maxit,
                               method = control$method) #LLT Cholesky decomposition, faster

  #Check for fitting convergence
  if (!fit1$converged) warning("The fitting algorithm did not converge.\n",
                               immediate. = TRUE)

  #J = X^tWX
  # w <- fit1$weights #W = diag(w_i)

  #without intercept:
  Xbeta <- X[, -1, drop = FALSE]
  beta <- fit1$coefficients[-1]

  eta <- drop(Xbeta %*% beta)

  # Jaa <- sum(w) #information associated to intercept
  D <- Fobs.info(fit1, family)
  Jaa <- sum(D)

  # Wald statistic
  # Q <- sum(w * eta^2) - sum(w * eta)^2 / Jaa
  Q <- sum(D * eta^2) - sum(D * eta)^2 / Jaa

  #Get information matrix for null model
  # Jaa0 <- sum(glm0$weights)
  D0 <- Fobs.info(glm0, family)
  Jaa0 <- sum(D0)

  #change in deviance
  z <- glm0$deviance - fit1$deviance

  ## Log-factor de Bayes: Li y Clyde
  logphi.num <- BAS::phi1(b/2, r, (a + b + k)/2, s + Q/(2*v), 1 - ka,
                          log = TRUE)

  logphi.den <- BAS::phi1(b / 2, r, (a + b)/2, s/(2*v), 1 - ka,
                          log = TRUE)

  #Bayes factor derived from approximate marginal likelihood (Li and Clyde, 2018)
  lBFi0 <- z/2 - Q/(2*v) -
    log(Jaa / Jaa0)/2 - k/2 * log(v) +
    lbeta((a + k)/2, b/2) + logphi.num -
    lbeta(a/2, b/2) - logphi.den

  return(lBFi0)
}

#' Logarithm of the g-prior Bayes factor in glm
#'
#' Computes the logarithm of the Bayes factors derived from a given g-prior
#' for complete data in generalized linear models.
#'
#' @param y Response variable in the linear model.
#' @param X Full imputed covariance matrix for a particular model including
#' the fixed terms and the intercept.
#' @param family String, function or the call to a family function among
#' \code{\link[stats]{family}} to specify the error distribution and link
#' function to be used in the model. Only available the implemented
#' families in \pkg{BAS}: \code{binomial(link = "logit")},
#' \code{poisson(link = "log")} and \code{Gamma(link = "log")}.
#' @param prior.betas Prior distribution for model-specific coefficients.
#' Options include \code{\link[BAS]{g.prior}}, \code{\link[BAS]{CCH}},
#' \code{\link[BAS]{robust}} and \code{\link[BAS]{intrinsic}} among others.
#' See \code{\link[BAS]{BAS}} for more details.
#' @param logmargnull Log-marginal likelihood of the null model considered.
#' @param k Number of model-specific coefficients.
#' @param p0 Number of fixed covariates (including the intercept).
#' @param weights NULL or numeric vector of the same length as \code{y} to
#' specify the weights to be used in the glm fitting process.
#' @param offset NULL or a numeric vector of the same length as \code{y} to
#' specify an a priori known component included in the glm fitting process.
#' @param control List of parameters for controlling the glm fitting process.
#' It is set to \code{[stats]{glm.control()}} by default.
#' @param c_glm.marg Function to compute log-marginal.
#'
#' @return \code{BF.gprior.glm} returns, in logarithmic scale, the exact
#' value of the Bayes factor derived from assigning a chosen g-prior by
#' \code{prior.betas} in generalized linear models for a given model through
#' \code{X}.
#'
#' @author Carolina Mulet
#' Maintainer: <Carolina.Mulet1@@alu.uclm.es>
#'
#' @seealso Use \code{\link[meantobvs]{lBF.av}} to compute the average Bayes
#' factor for missing data. Use \code{\link[meantobvs]{bvs.glm}} for
#' an exact computation of the model  posterior distribution in the VS problem
#' (recommended when p<20).
#'
#' @examplesIf interactive()
#' #Indian Prime Diabetes Data
#'
#' f <- Outcome ~ Pregnancies + Glucose + BloodPressure + SkinThickness + Insulin
#' imp1 <- mice.imputation(model.frame(f, diabetes, na.action = NULL),
#'   n.imp = 1, seed = 1)$imputation.array[,,1]
#'
#' glmnull <- glm(Outcome ~ 1, data = diabetes, family = binomial(), y = TRUE)
#' lBF <- meantobvs:::BF.gprior.glm(y = glmnull$y, X = imp1,
#'   family = binomial(), logmargnull = logmargnull.glm(glmnull))
#'
#' @references Clyde, M (2025) BAS: Bayesian Variable Selection and Model Averaging using
#' Bayesian Adaptive Sampling. R package version 2.0.2
#' <https://CRAN.R-project.org/package=BAS>.
#'
#' Li, Y. and Clyde, M. (2018)<DOI:10.1080/01621459.2018.1469992> Mixtures of
#' g-priors in Generalized Linear Models. Journal of the American Statistical
#' Association. 113: 1828-1845
#'
#'
#' @keywords internal
BF.gprior.glm <- function(y, X, family = binomial(link = "logit"),
                          prior.betas = BAS::g.prior(length(y)),
                          logmargnull,
                          k = ncol(X)-p0, p0 = 1L,
                          weights = rep(1, length(y)),
                          offset = rep(0, length(y)),
                          control = glm.control(),
                          laplace = 0L,
                          c_glm.marg = function() utils::getFromNamespace("C_glm_deterministic", "BAS")) {

  initprob <- c(rep(1.0, p0), rep(.5, k)) #first p0 columns of X are the fixed covariates
  fit1 <- .Call(c_glm.marg(), Y = y, X = X, Roffset = offset, Rweights = weights,
                Rprobinit = initprob, Rmodeldim = 0L, modelprior = BAS::uniform(),
                betaprior = prior.betas, family = family, Rcontrol = control, Rlaplace = laplace)

  # BFi0 <- exp(fit1$logmarg - logmargnull)
  lBFi0 <- fit1$logmarg - logmargnull
  return(lBFi0)
}

#' Checks that the Bayes factor computation method given by BF.method and
#' prior.betas is implemented and returns the function to use for Bayes factor
#' computation on lm.
#'
#' @keywords internal
checkforprior.betas.lm <- function (BF.method, prior.betas, n, p, p0, y, SS0) {

  if (BF.method %notin% c("BIC", "TBF", "gprior")) {
    stop("Only BF approximations 'BIC', 'TBF' and 'gprior' supported.")
  }

  if (is.null(prior.betas)) prior.betas <- "gZellner" #default opction

  switch (BF.method,
          BIC = {
            BF.method.f <-
              function (k, X) BF.BIC.lm(y = y, X,
                                        SS0 = SS0, n = n, k, p0 = p0)
            },

          TBF = {
            #first build the function needed to compute log-TBF
            lTBF.method <- lTBF.function(prior.betas, n, p, devnull = 0)

            BF.method.f <- function (k, X) {
              BF.TBF.lm(y = y, X,
                        SS0 = SS0,
                        lTBF.method = lTBF.method,
                        n = n, k, p0 = p0)}
            },

          gprior = {
            #Get string needed to define the prior
            prior.betas <- prior.betas.gprior.lm(prior.betas)

            BF.method.f <- ifelse(
              prior.betas != "flsBF",
              function (k, X) BF.gprior.lm(y = y, X, SS0 = SS0,
                                           prior.betas = prior.betas,
                                           n = n, k, p0 = p0),
              function (k, X) BF.FLS.lm(y = y, X, SS0 = SS0, dmax = p + p0,
                                        n = n, k, p0 = p0)
            )
          }
  )
  return(BF.method.f)
}

#' Checks that the Bayes factor computation method given by BF.method and
#' prior.betas is implemented and returns the function to use for Bayes factor
#' computation on glm.
#'
#' @keywords internal
checkforprior.betas.glm <- function (BF.method, prior.betas, n, p, p0, y,
                                     glmnull, useBAS) {

  if (BF.method %notin% c("BIC", "TBF", "gprior")) {
    stop("Only BF approximations 'BIC', 'TBF' and 'gprior' supported.")
  }

  if(is.null(prior.betas)) prior.betas <- "gZellner" #default option

  devnull <- glmnull$deviance #deviance of the null model

  #BAS logmarginal computation: faster for the families available in BAS
  if (useBAS) {
    #Get functions needed to compute logmarginals and define the prior
    c_glm.marg <- function() utils::getFromNamespace("C_glm_deterministic", "BAS")
    prior.betas <- prior.betas.gprior.glm(prior.betas, n, p)

    logmargnull <- logmargnull.glm(glmnull, y, prior.betas, 0L, c_glm.marg)

    BF.method.f <- function (k, X, ...) {
      BF.gprior.glm(y = y, X,
                    family = glmnull$family,
                    prior.betas = prior.betas,
                    logmargnull = logmargnull,
                    k, p0 = p0,
                    weights = glmnull$prior.weights,
                    offset = glmnull$offset,
                    control = glmnull$control,
                    c_glm.marg = c_glm.marg)}

  } else {
    control <- glm2fastglm.control(glmnull$control)

    switch (BF.method,
            BIC = {
              BF.method.f <- function (k, X, fitstart) {
                BF.BIC.glm.fit(y = y, X,
                               family = glmnull$family,
                               devnull = devnull,
                               n = n, k,
                               weights = glmnull$prior.weights,
                               offset = glmnull$offset,
                               control = control,
                               fitstart)}
            },

            TBF = {
              #first build the function needed to compute log-TBF
              lTBF.method <- lTBF.function(prior.betas, n, p, devnull)

              BF.method.f <- function (k, X, fitstart) {
                BF.TBF.glm.fit(y = y, X,
                               family = glmnull$family,
                               n = n, k, lTBF.method = lTBF.method,
                               weights = glmnull$prior.weights,
                               offset = glmnull$offset,
                               control = control,
                               fitstart)}
            },

            gprior = {
              #Get list of arguments needed to define the prior
              prior.betas.args <- prior.betas.gprior.glm.fit(prior.betas, n, p)

              #If there are fixed variables, the log-Bayes factor of Mi respect Mnull
              #is computed as: lBinull = lBi0 - lBnull0, i. e. is the difference of
              #log-Bayes factors of Mi and Mnull model respect the one with only the intercept
              if (glmnull$rank > 1) {
                #fit model with just intercept term
                glm0 <- glm.fit(y = y, x = rep(1, n),
                                family = glmnull$family,
                                weights = glmnull$prior.weights,
                                offset = glmnull$offset,
                                control = glmnull$control)
              } else {
                glm0 <- glmnull # ¡it is already the one with just the intercept
                lBFnull0 <- 0
              }

              #Give family object needed functions to compute BFs
              #(and check if given family-link combination is ok)
              family.upd <- addfamilyfun(glm0$family)

              #select method for prior coefficients
              switch(prior.betas.args$type,
                     fixed = BF.method.glm.fit <- function (k, X, fitstart) {
                       BF.gprior.glm.fit(y = y, X, glm0 = glm0,
                                         family = family.upd,
                                         g = prior.betas.args$g,
                                         n = n, k,
                                         weights = glmnull$prior.weights,
                                         offset = glmnull$offset,
                                         control = control,
                                         fitstart)
                       },

                     `hyper-g` = BF.method.glm.fit <- function (k, X, fitstart) {
                       BF.hypergprior.glm.fit(y = y, X, glm0 = glm0,
                                              family = family.upd,
                                              prior.betas.args = prior.betas.args[-1],
                                              n = n, k,
                                              weights = glmnull$prior.weights,
                                              offset = glmnull$offset,
                                              control = control,
                                              fitstart)
                      }
                     )

              if (glmnull$rank > 1) {
                #compute log-BF for model with just intercept term
                lBFnull0 <- BF.method.glm.fit(k = ncol(glmnull$x) - 1,
                                              X = glmnull$x, fitstart = NULL)
              }

              #lBF12 = lBF10 - lBF20 by Bayes factor transitivity
              BF.method.f <- function (k, X, fitstart) {
                BF.method.glm.fit(k + p0 - 1, X, fitstart) - lBFnull0
              }
            }
    )
  }

  return(BF.method.f)
}

#' @keywords internal
lTBF.function <- function(prior.betas, n, p, devnull) {
  #Check whether or not prior is among the available ones and return the TBF
  #function needed to compute it

  switch (prior.betas,

          gZellner = {lTBF.method <- #fixed g=n
            function(k, dev) lTBF.gfixed(g = n, k, dev, devnull = devnull)},
          Liangetal = {lTBF.method <- #random g: hyper-g/n with a=3
            function(k, dev) lTBF.hyperg(k, dev, devnull = devnull)},
          ZellnerSiow = {lTBF.method <- #adapted Z-S by trG
            function(k, dev) lTBF.grandom(a = .5, b = (n+3)/2, k, dev, devnull = devnull)},
          FLS = {lTBF.method <- #fixed Benchmark prior: g=max(n, p*p)
            function(k, dev) lTBF.gfixed(g = max(n, p^2), k, dev, devnull = devnull)},

          # Robust, intrinsic.WNC and IHG non-available for TBF approximation
          stop("Prior.betas must be one of 'gZellner', 'Liangetal', 'ZellnerSiow' or 'FLS'",
               "when using TBF method.\n")
  )
  return(lTBF.method)
}

#' @keywords internal
prior.betas.gprior.lm <- function(prior.betas) {
  #Check whether or not prior is among the available ones and return the string
  #needed to compute it

  switch (prior.betas,
          # change the string for the corresponding tag in BayesVarSel code

          gZellner = {prior.betas <- "gBF"}, #fixed g=n
          Robust = {prior.betas <- "RobustBF"},
          #random g: criteria-based prior from Bayarri et al (2012)
          Liangetal = {prior.betas <- "LiangBF"}, #random g: hyper-g/n with a=3
          ZellnerSiow = {prior.betas <- "ZSBF"}, #random g: cauchy prior
          FLS = {prior.betas <- "flsBF"}, #fixed g Benchmark prior: g=max(n, p*p)
          `intrinsic.MGC` = {prior.betas <- "intrinsicBF"},
          #intrinsic prior from Moreno, Giron, Casella (2015)
          IHG = {prior.betas <- "geointrinsicBF"}, #intrinsic hyper-g prior

          stop("prior.betas must be one of 'gZellner', 'Robust', 'Liangetal', 'ZellnerSiow',\n",
               "'FLS', 'intrinsic.MGC' or  'IHG' when using gprior method.\n")
  )
  return(prior.betas)
}

#' @keywords internal
prior.betas.gprior.glm <- function(prior.betas, n, p) {
  #Check whether or not prior is among the available ones and return the function
  #needed to compute it

  switch (prior.betas,

          gZellner = {prior.betas <- BAS::g.prior(g = as.numeric(n))}, #fixed g=n
          Robust = {prior.betas <- BAS::robust(as.numeric(n))}, #random g
          Liangetal = {prior.betas <- BAS::hyper.g(alpha = 3)},
          #random g: hyper-g with a=3
          ZellnerSiow = {prior.betas <-
            BAS::CCH(alpha = 1, beta = 2, s = n + 3)}, #adapted Z-S by trG
          FLS = {prior.betas <- BAS::g.prior(g = max(n, p^2))},
          #fixed Benchmark prior: g=max(n, p*p)
          `intrinsic.WNC` = {prior.betas <- BAS::intrinsic(as.numeric(n))},
          #intrinsic prior from Womack, Novelo and Casella (2014)

          # IHG non-available for gprior
          stop("Prior.betas must be one of 'gZellner', 'Robust', 'Liangetal', 'ZellnerSiow',",
               "'FLS' or 'intrinsic.WNC' when using gprior method.\n")
  )
  return(prior.betas)
}

#' @keywords internal
prior.betas.gprior.glm.fit <- function(prior.betas, n, p) {
  #Check whether or not prior is among the available ones and return list of
  #arguments needed to compute it

  switch (prior.betas,

          gZellner = {prior.betas.args <- list(type = "fixed", g = n)}, #fixed g=n
          Robust = {prior.betas.args <- #random g
            list(type = "hyper-g",
                 a = 1, b = 2, r = 3/2, s = 0,
                 v = function (k) (n + 1)/(k + 1), ka = 1)},
          Liangetal = {prior.betas.args <-
            list(type = "hyper-g",
                 a = 1, b = 2, r = 0, s = 0, v = 1, ka = 1)},
          #random g: hyper-g/n with a=3
          ZellnerSiow = {prior.betas.args <- #adapted Z-S by trG
            list(type = "hyper-g",
                 a = 1, b = 2, r = 0, s = n + 3, v = 1, ka = 1)},
          FLS = {prior.betas.args <- list(type = "fixed", g = max(n, p^2))},
          #fixed Benchmark prior: g=max(n, p*p)
          `intrinsic.WNC` = {prior.betas.args <- #intrinsic prior from Womack, Novelo and Casella (2014)
            list(type = "hyper-g",
                 a = 1, b = 1, r = 1, s = 0,
                 v = function (k) (n + k + 1)/(k + 1),
                 ka = function (k) (n + k + 1)/n)},

          # IHG non-available for gprior
          stop("Prior.betas must be one of 'gZellner', 'Robust', 'Liangetal', 'Zellner-Siow',",
               "'FLS' or 'intrinsic.WNC' when using gprior method.\n")
  )
  return(prior.betas.args)
}

#' Logarithm of the marginal distribution of null model in glm
#'
#' @keywords internal
logmargnull.glm <- function (glmnull, y = glmnull$y,
                             prior.betas = BAS::g.prior(g = length(y)),
                             laplace = 0L,
                             c_glm.marg =
                               function() utils::getFromNamespace("C_glm_deterministic", "BAS")) {

  if (glmnull$rank == 1) { #just the intercept is fixed

    loglik0 <- as.numeric(logLik(glmnull))
    Jaa0 <- sum(glmnull$weights)

    logmargnull <- loglik0 + .5 * log(2 * pi) - .5 * log(Jaa0)

  } else {
    if (is.null(glmnull$offset)) glmnull$offset <- numeric(length(y))

    logmargnull <- BF.gprior.glm(y = y, X = glmnull$x,
                                 family = glmnull$family,
                                 prior.betas = prior.betas,
                                 logmargnull = 0L,
                                 k = ncol(glmnull$x), p0 = 0L,
                                 weights = glmnull$prior.weights,
                                 offset = glmnull$offset,
                                 control = glmnull$control,
                                 laplace = laplace,
                                 c_glm.marg = c_glm.marg)
  }

  return(logmargnull)
}
