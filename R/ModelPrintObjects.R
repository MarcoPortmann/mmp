#!2 Extract model information from different classes
#
#
ModelPrintObject <- function(x) UseMethod('ModelPrintObject')
#


ConvertNewToOldMPO <- function(x)
{
  if (!("mmpVers" %in% names(x) && x$mmpVers >= 2))
    return(x)



  x$CoefValues <- array(NA, dim = c(nrow(x$Coefficients$Main),1,6), dimnames = list(rownames(x$Coefficients$Main), 'a', c("Coefficient", "SE", "TStat", "Normal", "PValue", "special")))
  x$CoefValues[,"a", "Coefficient"] <- x$Coefficients$Main$Coef
  if (length(x$Coefficients$Main$SE) > 0) x$CoefValues[,"a", "SE"] <- x$Coefficients$Main$SE
  if (length(x$Coefficients$Main$TStat) > 0) x$CoefValues[,"a", "TStat"] <- x$Coefficients$Main$TStat
  if (length(x$Coefficients$Main$Z) > 0) x$CoefValues[,"a", "Normal"] <- x$Coefficients$Main$ZStat
  if (length(x$Coefficients$Main$PValue) > 0) x$CoefValues[,"a", "PValue"] <- x$Coefficients$Main$PValue

  # CoefTemp <- array(x$Coefficients$Main, dim = c(nrow(x$Coefficients$Main), 1, ncol(x$Coefficients$Main)),
  #                   dimnames = list(rownames(x$Coefficients$Main), "a", colnames(x$Coefficients$Main)))

  x$DependentVariable <- x$DependendVariables$Main
#  browser()
  return(x)
}



ModelPrintObject.default <- function(x)
{
  return(x)
}

ModelPrintObject.rms <- function(x)
{
  x$Statistics <- append(x$Statistics,
                         list(
                           R2 = x$stats['R2'],
                           LogLik = x$stats['Model L.R.'],
                           G = x$stats['g'],
                           Sigma = x$stats['Sigma']
                         ))
  if (!is.null(x$orig.var) & is.null(x$clusterInfo))
  {
    x$Statistics <- append(x$Statistics,
                           list(
                             StandardErrorInfo = 'SE by robcov()'
                           ))
  }
  if (!is.null(x$clusterInfo))
  {
    x$Statistics <- append(x$Statistics,
                           list(
                             StandardErrorInfo = paste('Clustered on ', x$clusterInfo$name, '; n = ',  x$clusterInfo$n, sep ='')
                           ))
  }
  if (is.null(x$CoefValues))
  {
    CoefTemp     <- cbind(Coefficient = x$coefficients,
                          SE = sqrt(diag(x$var)),                 # Std. Error
                          TStat = NA,
                          Normal= x$coefficients/sqrt(diag(x$var)),  # Normal
                          PValue = NA,
                          SpecialText = NA
    )
    CoefTemp[, 'PValue'] <- unlist(lapply(CoefTemp[, 'Normal'], GetPValueFromNormal))
    x$CoefValues <- array(CoefTemp, dim = c(nrow(CoefTemp),1,6), dimnames = list(rownames(CoefTemp), 'a', colnames(CoefTemp)))
  }
  x$DependentVariable <- all.vars(x$terms)[1]
  return(x)
}

ModelPrintObject.rdrobust <- function(x, Stats = NULL)
{
  y <- list(Coefficients = list(Main = data.frame(Coef = as.numeric(x$coef),
                                                  SE = as.numeric(x$se),
                                                  ZStat = as.numeric(x$z),
                                                  PValue = as.numeric(x$pv),
                                                  row.names = rownames(x$coef)
                                                  )))
  #y$DependentVariables <- list(Main = x$???)
  y$Statistics <- list(NObs      = sum(x$N),
                       NObsLR    = list(x = x$N,collapse = ", "),
                       EffNObsLR = list(x = x$Nh, collapse = ", "),
                       Kernel    = x$kernel,
                       BWEst     = list(x = x$bws[1, ], collapse = ", "),
                       BWBias    = list(x = x$bws[2, ], collapse = ", "),
                       Rho       = list(x = c(x$bws[1, 1]/x$bws[2, 1], x$bws[1, 2]/x$bws[2, 2]), collapse = ", ")
                      )
  y$ModelType <- "local polynomial RD point estimators with robust bias-corrected confidence intervals"

  # to do: Defaults should be moved to a "GetModelDefaults" function
  y$Defaults <- list(Coefficients = c("Main"),
                     Statistics = c("NObs", "NObsLR", "EffNObsLR", "Kernel", "BWEst", "BWBias", "Rho"))
  y$mmpVers <- 2.0
  return(y)
}


# dev m.p.
ModelPrintDescription.rdrobust <- function(x)
{
 Stats <- data.frame(Stats = c("NObs", "NObsLR", "EffNObsLRKernel", "BWEst", "BWBias", "Rho"),
                     Default = c(rep(T, 6))
                     #LabelShort = c("N. obs.", ...)
                     #LabelLong = ?
                     #Language support? LaTeX?
                     # Formatting rules like collapse and number format (if not default shall be used)
                    )
  return(list(Statistics = Stats))
}


ModelPrintObject.lrm <- function(x)
{
  x$Statistics$Freq <- x$freq
  x$Statistics <- append(x$Statistics,
                         list(
                           NObs = x$stats['Obs'],
                           Brier = x$stats['Brier'],
                           Gamma = x$stats['Gamma'],
                           MaxDeriv = x$stats['Max Deriv'],
                           P = x$stats['P'],
                           C = x$stats['C'],
                           Dxy = x$stats['Dxy'],
                           Taua = x$stats['Tau-a'],
                           Gr = x$stats['gr'],
                           Gp = x$stats['gp'],
                           Freq = x$freq,
                           Weights = ifelse(is.null(x$weights),"no", "yes")
                         ))
  x$ModelType <- 'Logistic regression model'
  NextMethod()
}

ModelPrintObject.ols <- function(x)
{
  x$Statistics <- append(x$Statistics,
                         list(
                           NObs = x$stats['n'],
                           R2 = x$stats['R2'],
                           Sigma = x$stats['Sigma'],
                           DF = x$df.residual,
                           Weights = ifelse(is.null(x$weights),"no", "yes")
                         ))

  x$Statistics$AdjR2 <- 1 - (1-x$Statistics$R2)*(x$Statistics$NObs-1)/(x$Statistics$DF)

  CoefTemp     <- cbind(Coefficient = x$coefficients,
                        SE = sqrt(diag(x$var)),                 # Std. Error
                        TStat = x$coefficients/sqrt(diag(x$var)),
                        Normal= NA,  # Normal
                        PValue = NA,
                        SpecialText = NA
  )

  CoefTemp[, 'PValue'] <- unlist(lapply(CoefTemp[, 'TStat'], GetPValueFromTStat, DF =  x$Statistics$DF ))
  x$CoefValues <- array(CoefTemp, dim = c(nrow(CoefTemp),1,6), dimnames = list(rownames(CoefTemp), 'a', colnames(CoefTemp)))
  x$DependentVariable <- all.vars(x$terms)[1]
  x$ModelType <- 'Ordinary least squares regression model'
  NextMethod()
}

ModelPrintObject.felm <- function(x)
{
  xsum <- summary(x)

  x$Statistics <- append(x$Statistics,
                         list(
                           NObs =xsum$N,
                           R2 = xsum$r.squared,
                           Sigma = xsum$sigma,
                           DF = xsum$df[[1]],
                           AdjR2 = xsum$adj.r.squared,
                           Weights = ifelse(is.null(x$weights),"no", "yes")
                         ))

  CoefTemp     <- cbind(Coefficient = xsum$coefficients[, 1],
                        SE =  xsum$coefficients[, 2],                 # Std. Error
                        TStat =  xsum$coefficients[, 3],
                        Normal= NA,  # Normal
                        PValue =  xsum$coefficients[, 4],
                        SpecialText = NA
  )
  rownames(CoefTemp) <- rownames(xsum$coefficients)

  if (!is.null(x$fe) && length(x$fe)>0)
  {
    x$FixedEffects <- names(x$fe)
    CoefTemp <- rbind(CoefTemp,  matrix(NA, nrow = length(x$FixedEffects), ncol = 6, dimnames = list(paste0('FE>>', x$FixedEffects), NULL)))
  }

  x$CoefValues <- array(CoefTemp, dim = c(nrow(CoefTemp),1,6), dimnames = list(rownames(CoefTemp), 'a', colnames(CoefTemp)))
  x$DependentVariable <- x$lhs
  if (!is.null(x$clustervar)) x$Statistics$StandardErrorInfo <- paste0('Clustered on: ', paste(names(x$clustervar), collapse = ', '))

  x$ModelType <- 'Ordinary least squares regression model'
  NextMethod()
}

ModelPrintObject.summary.felm <- function(x)
{
  #xsum <- summary(x)

  x$Statistics <- append(x$Statistics,
                         list(
                           NObs =x$N,
                           R2 = x$r.squared,
                           Sigma = x$sigma,
                           DF = x$df[[1]],
                           AdjR2 = x$adj.r.squared,
                           Weights = ifelse(is.null(x$weights),"no", "yes")
                         ))

  CoefTemp     <- cbind(Coefficient = x$coefficients[, 1],
                        SE =  x$coefficients[, 2],                 # Std. Error
                        TStat =  x$coefficients[, 3],
                        Normal= NA,  # Normal
                        PValue =  x$coefficients[, 4],
                        SpecialText = NA
  )
  rownames(CoefTemp) <- rownames(x$coefficients)

  if (!is.null(x$fe) && length(x$fe)>0)
  {
    x$FixedEffects <- names(x$fe)
    CoefTemp <- rbind(CoefTemp,  matrix(NA, nrow = length(x$FixedEffects), ncol = 6, dimnames = list(paste0('FE>>', x$FixedEffects), NULL)))
  }

  x$CoefValues <- array(CoefTemp, dim = c(nrow(CoefTemp),1,6), dimnames = list(rownames(CoefTemp), 'a', colnames(CoefTemp)))
  x$DependentVariable <- x$lhs
  if (length(grep('obust', colnames(x$coefficients)[2]))>0) x$Statistics$StandardErrorInfo <- 'Robust standard errors'

  if (!is.null(x$clustervar)) x$Statistics$StandardErrorInfo <- paste0('Clustered on: ', paste(names(x$clustervar), collapse = ', '))

  x$ModelType <- 'Ordinary least squares regression model'
  NextMethod()
}



ModelPrintObject.plm <- function(x)
{
  warning('plm objects are processed slowly and incompletely.')
  xsum <- summary(x)
  x$Statistics <- append(x$Statistics,
                         list(
                           NObs = nobs(x),
                           R2 = as.numeric(xsum$r.squared['rsq']),
                           DF = x$df.residual,
                           FStat = xsum$fstatistic,
                           AdjR2 = as.numeric(xsum$r.squared['adjrsq']),
                           Weights = ifelse(is.null(x$weights),"no", "yes")
                         ))


  CoefTemp     <- cbind(Coefficient = xsum$coefficients[,1],
                        SE = xsum$coefficients[,2],                 # Std. Error
                        TStat =  xsum$coefficients[,3],
                        Normal= NA,  # Normal
                        PValue =  xsum$coefficients[,4],
                        SpecialText = NA
  )
  rownames(CoefTemp) <- rownames(xsum$coefficients)

  x$CoefValues <- array(CoefTemp, dim = c(nrow(CoefTemp),1,6), dimnames = list(rownames(CoefTemp), 'a', colnames(CoefTemp)))

  x$DependentVariable <- all.vars(x$formula)[1]
  x$ModelType <- paste0('Panel data:', paste(x$args[unlist(lapply(x$args, function(x)(!is.null(x))))], collapse = '; '))
  NextMethod()
}




ModelPrintObject.lm <- function(x)
{

  x$Statistics <- append(x$Statistics,
                         list(
                           NObs = length(residuals(summary(x))),
                           R2 = summary(x)$r.squared,
                           AdjR2 = summary(x)$adj.r.squared,
                           Sigma = summary(x)$sigma,
                           DF = summary(x)$df[2],
                           FStat = summary(x)$fstatistic[1],
                           Weights = ifelse(is.null(x$weights),"no", "yes")
                         ))

  # The following might looks complicated. The aim is to make sure that
  # NA variables are not dropped from the model.

  CoefTemp <- as.matrix(x$coefficients)
  CoefTemp <- merge(CoefTemp, summary(x)$coefficients[,c(2,3)], by = 'row.names', all.x = T, sort = F)
  rownames(CoefTemp) <- CoefTemp$Row.names
  CoefTemp$Row.names <- NULL

  loc <- grep("^\\(Intercept\\)$", rownames(CoefTemp), perl = T)
  if (length(loc)>0) rownames(CoefTemp)[loc] <- "Intercept"

  CoefTemp <- cbind(CoefTemp, Normal = NA, PValue = NA,  SpecialText = NA)
  colnames(CoefTemp) <- c('Coefficient', 'SE', 'TStat', 'Normal', 'PValue', 'SpecialText')

  CoefTemp[, 'PValue'] <- unlist(lapply(CoefTemp[, 'TStat'], GetPValueFromTStat, DF =  x$Statistics$DF ))

  x$CoefValues <- array(unlist(CoefTemp), dim = c(nrow(CoefTemp),1,6), dimnames = list(rownames(CoefTemp), 'a', colnames(CoefTemp)))
  x$DependentVariable <- all.vars(x$terms)[1]
  #bpres <- NA
  #try(bpres <- bptest(m$terms)[4], silent = T)

  if (is.null(x$ModelType))
    x$ModelType <- 'Ordinary least squares regression model'
  NextMethod()
}

ModelPrintObject.glm <- function(x)
{
  warning('Basic glm support: Not all model details are supported.')
  x$qr





  x$ModelType <-  paste('Ordinary least squares regression model: ', x$family$family, '/', x$family$link, sep ='')

  NextMethod()
}


ModelPrintObject.selection <- function(x)
{
  #x <- a
  xx <- summary(x)
  x$Statistics <- append(x$Statistics,
                         list(
                           NObs = x$param$nObs,
                           R2 = xx$rSquared$R2,
                           Sigma = x$sigma,
                           DF = x$param$df,
                           AdjR2 =  xx$rSquared$R2adj,
                           Mills = x$param$index$Mills,
                           Weights = ifelse(is.null(x$weights),"no", "yes")
                         ))


  CoefTemp     <- cbind(Coefficient = xx$estimate[, 'Estimate'],
                        SE = xx$estimate[, 'Std. Error'],
                        TStat = xx$estimate[, 't value'],
                        Normal= NA,  # Normal
                        PValue = xx$estimate[, 'Pr(>|t|)'],
                        SpecialText = NA
  )
  CoefSelection <- CoefTemp[x$param$index$betaS, ]
  CoefTemp <- CoefTemp[x$param$index$betaO, ]


  x$CoefValues <- array(CoefTemp, dim = c(nrow(CoefTemp),1,6), dimnames = list(rownames(CoefTemp), 'a', colnames(CoefTemp)))
  x$CoefValues.Selection <- array(CoefSelection, dim = c(nrow(CoefSelection),1,6), dimnames = list(rownames(CoefSelection), 'a', colnames(CoefSelection)))

  x$DependentVariable <- NA
  x$ModelType <- 'Heckit'
  return(x)

}




ModelPrintObject.glmerMod <- function(x)
{
  warning('glmerMod is only partially supported by mmp yet. Only fixed effects are printed.')

  y <- list()
  #  y$Statistics <- append(y$Statistics,
  #                         list(
  #                              R2 = x$stats['R2'],
  #                              LogLik = x$stats['Model L.R.'],
  #                              G = x$stats['g'],
  #                              Sigma = x$stats['Sigma']
  #                             ))

  CoefTemp     <- cbind(Coefficient = fixef(x),
                        SE = sqrt(diag(vcov(x, use.hessian = T))),
                        TStat = NA,
                        Normal= fixef(x)/sqrt(diag(vcov(x, use.hessian = T))),
                        PValue = NA,
                        SpecialText = NA
  )

  loc <- grep("^\\(Intercept\\)$", rownames(CoefTemp), perl = T)
  if (length(loc)>0) rownames(CoefTemp)[loc] <- "Intercept"

  CoefTemp[, 'PValue'] <- unlist(lapply(CoefTemp[, 'Normal'], GetPValueFromNormal))
  y$CoefValues <- array(CoefTemp, dim = c(nrow(CoefTemp),1,6), dimnames = list(rownames(CoefTemp), 'a', colnames(CoefTemp)))


  y$DependentVariable <- colnames(model.frame(x))[1]
  y$ModelType <- paste('Generalized linear mixed-effects models; link:',family(x)$link, sep =' ')
  return(y)
}

ModelPrintObject.lmerMod <- function(x)
{
  warning('lmerMod is only partially supported by mmp yet. Only fixed effects are printed.')
  #require(r2glmm)
  library("piecewiseSEM")
  y <- list()
  #  y$Statistics <- append(y$Statistics,
  #                         list(
  #                              R2 = x$stats['R2'],
  #                              LogLik = x$stats['Model L.R.'],
  #                              G = x$stats['g'],
  #                              Sigma = x$stats['Sigma']
  #                             ))

  CoefTemp     <- cbind(Coefficient = fixef(x),
                        SE = sqrt(diag(vcov(x, use.hessian = F))),
                        TStat = NA,
                        Normal= fixef(x)/sqrt(diag(vcov(x, use.hessian = F))),
                        PValue = NA,
                        SpecialText = NA
  )

  loc <- grep("^\\(Intercept\\)$", rownames(CoefTemp), perl = T)
  if (length(loc)>0) rownames(CoefTemp)[loc] <- "Intercept"

  CoefTemp[, 'PValue'] <- unlist(lapply(CoefTemp[, 'Normal'], GetPValueFromNormal))
  y$CoefValues <- array(CoefTemp, dim = c(nrow(CoefTemp),1,6), dimnames = list(rownames(CoefTemp), 'a', colnames(CoefTemp)))

  y$Statistics <- append(y$Statistics,
                         list(
                           NObs = nobs(x),
                           R2 = paste0(round(rsquared(x)$Marginal, digits=2), ', ', round(rsquared(x)$Conditional, digits=2))
                         ))




  y$DependentVariable <- colnames(model.frame(x))[1]
  y$ModelType <- paste('Linear mixed-effects models; link:',family(x)$link, sep =' ')
  return(y)
}




ModelPrintObject.boot <- function(x)
{

  warning('Preliminary support for boot only.')


  CoefTemp     <- cbind(Coefficient = x$t0,
                        SE = apply(x$t, 2, sd),
                        TStat = NA,
                        Normal= x$t0/apply(x$t, 2, sd),
                        PValue = NA,
                        SpecialText = NA
  )


  CoefTemp[, 'PValue'] <- unlist(lapply(CoefTemp[, 'Normal'], GetPValueFromNormal))
  x$CoefValues <- array(CoefTemp, dim = c(nrow(CoefTemp),1,6), dimnames = list(rownames(CoefTemp), 'a', colnames(CoefTemp)))

  x$ModelType <-'boot object'
  NextMethod()
}

ModelPrintObject.tobit <- function(x)
{
  warning('Preliminary support for tobit only.')

  x$Statistics <- append(x$Statistics,
                         list(
                           Iterations = x$iter,
                           LogLik = x$loglik[2],
                           Scale = x$scale,
                           DF =x$df
                         ))

  TempSE <- sqrt(diag(x$var))[-nrow(x$var)]
  CoefTemp     <- cbind(Coefficient = x$coefficients,
                        SE = TempSE,
                        TStat = NA,
                        Normal= x$coefficients/TempSE,
                        PValue = NA,
                        SpecialText = NA
  )

  CoefTemp[, 'PValue'] <- unlist(lapply(CoefTemp[, 'Normal'], GetPValueFromNormal))
  x$CoefValues <- array(CoefTemp, dim = c(nrow(CoefTemp),1,6), dimnames = list(rownames(CoefTemp), 'a', colnames(CoefTemp)))

  x$DependentVariable <- all.vars(x$terms)[1]
  x$ModelType = 'Tobit Regression'
  NextMethod()
}

ModelPrintObject.DiscreteEffects <- function(x)
{
  CoefTemp     <- cbind(Coefficient = x$Estimate,
                        SE = x$SE,
                        TStat = NA,
                        Normal= x$Estimate/x$SE,
                        PValue = x$pvalue,
                        SpecialText = NA
  )
  y <- list()
  y$CoefValues <- array(CoefTemp, dim = c(nrow(CoefTemp),1,6), dimnames = list(rownames(CoefTemp), 'a', colnames(CoefTemp)))
  rownames(y$CoefValues) <- sub('^DE', "", rownames(x) )
  y$ModelType <- 'Discrete Effect'
  return(y)
}

ModelPrintObject.CoefTest <- function(x)
{
  warning('Untested ModelPrintObject. Probably requires output to be a list.')
  CoefTemp     <- cbind(Coefficient = x$Estimate,
                        SE = x$SE,
                        TStat = NA,
                        Normal= x$Estimate/x$SE,
                        PValue = x$pvalue,
                        SpecialText = NA
  )
  CoefTemp[, 'PValue'] <- unlist(lapply(CoefTemp[, 'Normal'], GetPValueFromNormal))
  x$CoefValues <- array(CoefTemp, dim = c(nrow(CoefTemp),1,6), dimnames = list(rownames(CoefTemp), 'a', colnames(CoefTemp)))
  rownames(x$CoefValues) <- sub('^DE', "", rownames(x) )
  x$ModelType <- 'Coef Test'
  NextMethod()
}

ModelPrintObject.SummaryStat <- function(x)
{
  warning('Untested ModelPrintObject. Probably requires output to be a list.')
  CoefTemp     <- cbind(Coefficient = x,
                        SE = NA,
                        TStat = NA,
                        Normal= NA,
                        PValue = NA,
                        SpecialText = NA
  )
  if (is.null(names(x))) rownames(CoefTemp) <- rownames(x)
  x$CoefValues <- array(CoefTemp, dim = c(nrow(CoefTemp),1,6), dimnames = list(rownames(CoefTemp), 'a', colnames(CoefTemp)))
  rownames(x$CoefValues) <- sub('^DE', "", rownames(x) )
  x$ModelType <- 'SummaryStat'
  NextMethod()
}


ModelPrintObject.data.frame <- function(x)
{
  #  warning('Untested ModelPrintObject. Probably requires output to be a list.')
  #  CoefTemp     <- cbind(Coefficient = x[,1],
  #                        SE = NA,
  #                        TStat = NA,
  #                        Normal= NA,
  #                        PValue = NA,
  #                        SpecialText = NA
  #                       )
  #  if (is.null(names(x))) rownames(CoefTemp) <- rownames(x)
  #  x$CoefValues <- array(CoefTemp, dim = c(nrow(CoefTemp),1,6), dimnames = list(rownames(CoefTemp), 'a', colnames(CoefTemp)))
  #  rownames(x$CoefValues) <- sub('^DE', "", rownames(x) )
  #  x$ModelType <- 'Generic data.frame'
  NextMethod()
}

