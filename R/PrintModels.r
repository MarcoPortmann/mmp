#! multi-model print

# Author: Marco Portmann
# Last change: 02.11.2017
#
#


#.........................................................
#
#!2 Extract model information from different classes
#
#
ModelPrintObject <- function(x) UseMethod('ModelPrintObject')
#

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


#.........................................................



GetPValue <- function(Stats, StatsType)
{
  if (is.null(Stats) || is.na(Stats) || is.null(StatsType) || is.na(StatsType)) return(NA)
  if (StatsType == -1)
  {
    return(2* (1- pnorm(abs(Stats))))
  }
  return(2* (1- pt(abs(Stats), StatsType )))
}


GetPValueFromNormal <- function(Stats)
  {
    if (is.null(Stats) || is.na(Stats)) return(NA)
      return(2* (1- pnorm(abs(Stats))))
  }

GetPValueFromTStat <- function(Stats, DF)
{
  if (is.null(Stats) || is.na(Stats)) return(NA)
    return(2* (1- pt(abs(Stats), DF)))
}





#!2 GetStars.Text
  GetStars.Text <- function(PValue, StarBreaks=c(-Inf, 0.01, 0.05, 0.1, Inf))
  {
    if (is.null(PValue) || is.na(PValue) || (!is.numeric(PValue) && !is.numeric(type.convert(PValue)))) return('')
    if (!is.numeric(PValue)) PValue <- type.convert(PValue)
    if (is.nan(PValue)) return('')
    return(paste(rep("*", max(4- as.numeric(cut(PValue, StarBreaks)), 0)), collapse=""))
  }


#!2 FormatNumber
  FormatNumber <- function(Num, Decimals = 4, IsPercent = F)
  {
    if (is.null(Num) || is.na(Num) || (!is.numeric(Num) && !is.numeric(type.convert(Num)))) return('')
    if (!is.numeric(Num)) Num <- type.convert(Num)
    if (IsPercent==T) Num <- Num * 100
      spstring <- paste("%.",Decimals,"f", sep="")
      res <- sprintf(spstring, Num)
      if (is.na(Num)==F){
        if (abs(Num)<=10^(-(Decimals-1))) {
           spstring <- paste("%.",1,"e", sep="")
           res <- sprintf(spstring, Num)
        }
      }else{res <- "NA"}

      if (IsPercent==T) res  <- paste(res, "%", sep="")
      return(res)
   }

#!2 FormatCoefficientOutput.Text

  # c('Coefficient', 'SE', 'SomeStat', 'StatIndicator', 'PValue', 'SpecialText')
  FormatCoefficientOutput.Text <- function(Coefficient = NULL, SE = NULL, PValue = NULL, SpecialText = NULL, PrintStat = 'SE', Digits = 4, IsPercent = F, CoefficientOnly = F, SpecialTextList = list('-3' = '', '-2' = 'no', '-1' = 'yes'))
  {
    if (length(Coefficient)>1)
    {
      if (!is.null(names(Coefficient)))
      {
        SE = Coefficient['SE']
        PValue = Coefficient['PValue']
        TValue = Coefficient['TValue']
        SpecialText =  Coefficient['SpecialText']

        Coefficient = Coefficient['Coefficient']
      }else
      {
        SE = Coefficient[2]
        PValue = Coefficient[5]
        # TValue = Coefficient[3]         not supported yet !
        SpecialText =  Coefficient[6]
        Coefficient = Coefficient[1]
      }
    }
    if (!is.null(SpecialText) &&  !is.na(SpecialText) && SpecialText != 0)
      return(ifelse(any(indx <- (as.character(SpecialText) == names(SpecialTextList))),  as.character(SpecialTextList[indx]), ''))


    if (is.null(Coefficient) || is.na(Coefficient)) return('')
    if ((!is.numeric(Coefficient) && !is.numeric(type.convert(Coefficient)))) return(Coefficient)
    if ((!is.numeric(Coefficient) && is.numeric(type.convert(Coefficient)))) Coefficient <- as.numeric(Coefficient)
    # Check PrintStat
    StatVal <- ''
    if (PrintStat == 'SE') {StatVal <- SE; BracketLeft <- '('; BracketRight <- ')'}
    if (PrintStat == 'PValue') {StatVal <- PValue; BracketLeft <- '['; BracketRight <- ']'}



    if (!CoefficientOnly && !is.null(StatVal) && !is.na(StatVal) && (is.numeric(StatVal) || is.numeric(type.convert(StatVal))))
    {
      return(paste(FormatNumber(Coefficient, Digits), GetStars.Text(PValue) ,'\n', BracketLeft, FormatNumber(StatVal, Digits), BracketRight, sep =''))
    }
      return(FormatNumber(Coefficient, Digits))
  }




#! Multi Model Print 3

#
#    CoefNames <- list(Var1 = ('^InAfterGroup$'),
#                      Var2 = c('Irgendwas', 'Mehr'),
#                      c('bi', 'ba'),
#                      'Inter1 * Inter2',
#                      list('Var1', 'Var2'),
#                       'Var1 x Var2' = list('Var1', 'Var2'),
#                       'Ref Day FE'  = paste('(^', 'ClosestReferendumString', '=.)|(^factor\\(', 'ClosestReferendumString', '\\).)', sep =''),
#                       '<Generic>' = '^(.*?)$'
#                      )

#' Multi-model printing
#'
#' @param AllModels List containing all model output objects. If the objects are named, the names will be used as headers.
#' @param CoefNames List with coefficient names.
#' @param Decimals
#' @param WriteToClipboard
#' @param CoefOnly
#' @param CoefTests
#' @param Factors
#' @param BottomMatter
#' @param ModelNames
#' @param ShowDependentVariable
#' @param ShowModelType
#'
#' @return
#' @export
#'
#' @examples
mmp <- function(

                 AllModels,
                 CoefNames = '<Generic>',

                Decimals=4,
                WriteToClipboard = F,

                CoefOnly = F,
                CoefTests = NULL,
              #  CoefAsPercent = F, To do: Write coefficients as percentages.
              #                            This feature makes sense e.g. for discrete effects.
                Factors = NULL,

                BottomMatter = c("r.squared", "Brier", "N"),
                ModelNames = NULL,
                ShowDependentVariable = F,
                ShowModelType = F
                #LaTeX=F,
                #landscape=F,  caption, caption.location="bottom",
                #label="tab:", append=F, factors = NA,       n.rgroup=0, rgroup=NA,
                ){


require(Hmisc)
require(lmtest)
require(abind)
LaTeX = F

  trim <- function (x) gsub("^\\s+|\\s+$", "", x)

  SplitInteractions <- function(CoefText)
  {
    return(unlist(lapply(strsplit(CoefText, '(\\*)|(:)'), trim)))
  }

  IsRegEx <- function(x)
  {
    # To do: at the moment, it is an incomplete test!
    #   But I don't want to exclude stuff that is allowed in regular regression names.
    return(length(grep('(\\^)|(\\$)', x, perl = T))>0)
  }

  PrepareCoefNames <- function(){
    #CoefNames <- c('Intercept', 'Test2')
    CoefTemp <- CoefNames

    if (is.null(CoefTemp)) CoefTemp <- list('<Generic>' = '(.*?)')
    if (!is.list(CoefTemp)) CoefTemp <- as.list(CoefTemp)

    ConcatenateRegexOrC <- function(x){
      # x <- c("Irgendwas", "Mehr")  #  x <- c("^Irgendwas$", "Me^hr")
      if ( length(x) == 1)
      {
        return(ifelse(IsRegEx(x), x, paste('^', x, '$', sep ='')))
      }
      if (IsRegEx(x)) stop('Regular expressions are not allowed when more than one string per coefficient is provided.')

      return(ifelse(IsRegEx(x), x, paste(paste('^', x, '$', sep =''), collapse='|') ))
    }

    ConcatenateFactorStrings <- function(x){
      if ( length(x) == 1)
      {                                # '(^', x, '=.         <-- sometimes there seems to be no =
     return(ifelse(IsRegEx(x), x, paste('(^', x, '.)|(^factor\\(', x, '\\).)|(^FE>>',x,'$)',  sep ='')))
      }
      if (IsRegEx(x)) stop('Regular expressions are not allowed when more than one string per coefficient is provided.')

      return(ifelse(IsRegEx(x), x, paste(paste('(^', x, '=.)|(^factor\\(', x, '\\).)|(^FE>>',x,'$)',  sep =''), collapse='|') ))
    }

    NoName <- function(x) return(is.null(names(x)) || is.na(names(x)) || names(x) == '')

    PrepareCoefName <- function(i, OriginList)
    {
      CoefItem <- OriginList[i]


      if (length(CoefItem[[1]])>1)
      {
        # Interactions - already defined as interaction
        if (is.list(CoefItem[[1]]))
        {
          if (NoName(CoefItem)) names(CoefItem) <- paste(CoefItem[[1]], collapse = ' * ')
          CoefItem[[1]] <- lapply(CoefItem[[1]],  ConcatenateRegexOrC)
        }else
        {
          if (NoName(CoefItem)) names(CoefItem) <- CoefItem[[1]][[1]]
          CoefItem[[1]] <- ConcatenateRegexOrC(CoefItem[[1]])
        }
        return(CoefItem)
      }

      if (!IsRegEx(CoefItem[[1]]))
      {
        if (CoefItem[[1]] == '<Generic>' && (NoName(CoefItem)))
        {
          names(CoefItem) <- '<Generic>'
          CoefItem[[1]] <- '(.*?)'
          return(CoefItem)
        }

        CoefItemText <- SplitInteractions(CoefItem[[1]])
        if (length(CoefItemText)>1)
        {
           if (NoName(CoefItem)) names(CoefItem) <-  paste(CoefItemText, collapse = ' * ')
          CoefItem[[1]] <- lapply(CoefItemText,  ConcatenateRegexOrC)
          return(CoefItem)
        }
      }
      if (NoName(CoefItem)) names(CoefItem) <- CoefItem[[1]][[1]]
        CoefItem[[1]] <- ifelse( (names(CoefItem) %in% names(Factors)), ConcatenateFactorStrings(CoefItem[[1]]), ConcatenateRegexOrC(CoefItem[[1]]))
      return(CoefItem)
    }



    for (i in 1:length(CoefTemp))
    {
      ItemTemp <- PrepareCoefName(i, OriginList=CoefTemp)
      CoefTemp[[i]] <- ItemTemp[[1]]
      names(CoefTemp)[i] <- names(ItemTemp)
    }


    if (length(Factors)>0)
    for (i in 1:length(Factors))
    {
      if (names(Factors)[i] %in% names(CoefTemp)) next

      ItemTemp <- PrepareCoefName(i, OriginList=Factors)
      CoefTemp <- c(CoefTemp, ItemTemp)
    }

    return(CoefTemp)
  }
#CoefNames <- PrepareCoefNames()





  GetPrettyCoefFromInputList <- function(CoefText)
  {
    #CoefText <- 'InAfterGroupXX'
    for (i in 1:length(CoefNames))
    {
      #CN2 <- CoefNames[[1]]
      if (names(CoefNames)[i] == '<Generic>') next
      if (length(CoefNames[[i]])>1) next
      if (length(grep(CoefNames[[i]],  CoefText, perl = T)) == 1)
      {
        return(names(CoefNames)[i])
      }
    }

    return(CoefText)
  }

  GetInteractionName <- function(ModelCoef)
  {
    ModelCoef <- sort(ModelCoef)
    for (i in 1:length(CoefNames))
    {
      #CN2 <- CoefNames[[1]]
      if (length(CoefNames[[i]])<2) next
      if ( length(CoefNames[[i]]) != length(ModelCoef) ) next
      Matched <- c()
      for (IEPart in CoefNames[[i]])
      {
        MatchIE <- grep(IEPart,  ModelCoef, perl = T)
        if (length(MatchIE)==0) break
        if (!all(MatchIE %in% Matched))
          Matched <- c(Matched, MatchIE[!MatchIE %in% Matched])
      }
        if (length(Matched) == length(CoefNames[[i]]))
          return(names(CoefNames)[i])
    }
    return(paste(ModelCoef, collapse = ' * '))
  }

GetPrettyCoef <- function(ModelCoef)
{
  #ModelCoef <- ModelCoefs[1]   # ModelCoef <- ModelCoefs[23]
  ModelCoef <- sort(SplitInteractions(ModelCoef))
  tmp <- list(lapply(ModelCoef, GetPrettyCoefFromInputList))
  names(tmp)[[1]] <- ifelse( length(tmp[[1]]) == 1, tmp[[1]][[1]], GetInteractionName(unlist(tmp)))
  return(tmp)
}

GetPrettyCoefList <- function(ModelCoefs)
{
  #ModelCoefs <- names(coef(RegTab501))
  #tst <- GetPrettyCoef(ModelCoefs[23])
  #c(GetPrettyCoef(ModelCoefs[1]), GetPrettyCoef(ModelCoefs[1]))
  return(unlist(lapply(ModelCoefs, GetPrettyCoef), recursive = F))
}


# summary statistics: variable names
# first col.: use this names in the function argument list
# second col.: internal name
# third col.: text which appears in the table
# fourth col.:  text which appears in a LaTeX table
# ... could add columns with long/short names


# ---------------

# The following matrix is a mess and the approach should be redisigned:
# - It is incomplete, only a fraction of stats provided by the supported models are incorporated.
# - As the number of models may grow over time and the user should have the possibility to add
#   support for more classes to ModelPrintObject, a static list might be the wrong approach.
# - At least the distinction between internal and external is unhandy.
# - There should be a more flexible approach that takes into account
#   - language support
#   - output specific needs (e.g., formatting of math expressions in LaTeX, HTML etc.)
#   - need for short/long names


# ---------------
ModelInfo <-
matrix(
c(
"Sigma"                 , "Sigma"               , "Sigma"             , "$\\sigma$"        ,
"DF"                    , "DF"                  , "DF"                , "DF"               ,
"f-statistics"          , "f-statistics"        , "F-stat."           , "F-stat."          ,
"Breusch-P."            , "Breusch-P."          , "Breusch-P."        , "Breusch-P."       ,
"N"                     , "NObs"                , "n. Obs."           , "n. Obs."          ,
"r.squared"             , "R2"                  , "R2"                , "$R^2$"            ,
"LogLik"                , "LogLik"              , "Log-Likelihood"    , "Log-Likelihood"   ,
"Brier"                 , "Brier"               , "Brier"             , "Brier"            ,
"G"                     , "G"                   , "G"                 , "G"                ,
"Gamma"                 , "Gamma"               , "Gamma"             , "Gamma"            ,
"MaxDeriv"              , "MaxDeriv"            , "MaxDeriv"          , "MaxDeriv"         ,
"P"                     , "P"                   , "P"                 , "P"                ,
"C"                     , "C"                   , "C"                 , "C"                ,
"Dxy"                   , "Dxy"                 , "Dxy"               , "Dxy"              ,
"Taua"                  , "Taua"                , "Taua"              , "Taua"             ,
"Freq"                  , "Freq"                , "Freq"              , "Freq"             ,
"adj.r.squared"         , "AdjR2"               , "Adj. R2"           , "Adj. $R^2$"       ,
"FStat"                 , "FStat"               , "FStat"             , "FStat"            ,
"Iterations"            , "Iterations"          , "Iterations"        , "Iterations"       ,
"Scale"                 , "Scale"               , "Scale"             , "Scale"            ,
"StandardErrorInfo"     , "StandardErrorInfo"   , "Standard Errors"   , "Standard Errors"  ,
"Weights"               , "Weights"             , "Weights"           , "Weights"
  ), byrow=T, ncol=4)

ModelInfo <- data.frame(ModelInfo, stringsAsFactors =F)
colnames(ModelInfo) <- c('External', 'Internal', 'Output.Text', 'Output.LaTeX')

if (length(ShowDependentVariable) == 1 && is.logical(ShowDependentVariable) &&  ShowDependentVariable)
{
 DependentVarPretty <- list()
}else
{
  if (is.list(ShowDependentVariable))
  {
    DependentVarPretty <- ShowDependentVariable
    ShowDependentVariable <- T
  }else
  {
    if (length(ShowDependentVariable)>1 || !is.logical(ShowDependentVariable))
      stop('ShowDependentVariable must either be a list or true/false.')
  }
}

#-----------------------------------------------------------------------------------------------------------------------------
#
#
#
#
#
#!2 Function starts here

# Check whether requested bottom matter variables exist, raise error otherwise.
errm <- BottomMatter[!BottomMatter %in% ModelInfo$External]
if (length(errm)>0) stop(paste('Undefined model infos supplied:', paste(errm, collapse = ', ')), sep=' ')


# ModelPrintObject extracts all data (such as coefficients and test statistics) from the model and returns it
# in a standardizied frame which can be understood by multi model print.
models <- lapply(AllModels, ModelPrintObject)

# This loop adds the model names and replaces standard errors by those provided via CoefTest.
# Model names:
#      They are assigned in the following way:  1. Names provided by the user via ModelNames.
#                                               2. If the AllModels list contains named object, those names are used.
#                                               3. Add a name based on the column number (i)
# coeftest:
#
#


for (i in 1:length(AllModels))
{

    if (!is.null(CoefTests) && length(CoefTests) >= i && !is.null(CoefTests[[i]]) && !is.na(CoefTests[i]))
    {
# NEW
      CTNames <-  gsub('^\\(Intercept\\)$', 'Intercept', rownames(CoefTests[[i]]))

      if (length(DiffNames <- unique(c(setdiff(names(models[[i]]$CoefValues[, 1, 'Coefficient']), CTNames), setdiff(CTNames, names(models[[i]]$CoefValues[, 1, 'Coefficient']))))) > 0)
        warning(paste0('Differences between variables provided in AllModels and CoefTests found for model #',i , ' found: ', paste(DiffNames, collapse=', ')))

      if(all(!c('t value', 'z value') %in% colnames(CoefTests[[i]])))
        stop('Unknown coeftest.')

      CoefTemp     <- cbind(Coefficient = CoefTests[[i]][, 'Estimate'],
                            SE = CoefTests[[i]][, 'Std. Error'],
                            TStat = ifelse( ('t value' %in% colnames(CoefTests[[i]])), CoefTests[[i]][, 't value'], NA),
                            Normal= ifelse( ('z value' %in% colnames(CoefTests[[i]])), CoefTests[[i]][, 'z value'], NA),
                            PValue = CoefTests[[i]][, grep('Pr', colnames(CoefTests[[i]])) ],
                            SpecialText = NA
                           )
      rownames(CoefTemp) <- CTNames
      models[[i]]$CoefValues <- array(CoefTemp, dim = c(nrow(CoefTemp),1,6), dimnames = list(rownames(CoefTemp), 'a', colnames(CoefTemp)))
      # check: ^DE still required?
      rownames(models[[i]]$CoefValues) <- sub('^DE', "", rownames(models[[i]]$CoefValues))
# -----------

# OLD
#      models[[i]]$CoefValues[, 1, 'Coefficient'] <- CoefTests[[i]][, 'Estimate']
#      models[[i]]$CoefValues[, 1, 'SE'] <- CoefTests[[i]][, 'Std. Error']
#      if ('t value' %in% colnames(CoefTests[[i]]))
#      {
#        models[[i]]$CoefValues[, 1, 'TStat'] <- CoefTests[[i]][, 't value']
#        models[[i]]$CoefValues[, 1, 'PValue'] <- CoefTests[[i]][, 'Pr(>|t|)']
#      }else{if ('z value' %in% colnames(CoefTests[[i]]))
#              {
#                models[[i]]$CoefValues[, 1, 'Normal'] <- CoefTests[[i]][, 'z value']
#                models[[i]]$CoefValues[, 1, 'PValue'] <- CoefTests[[i]][, 'Pr(>|z|)']
#              }else{
#                stop('Unknown coeftest.')
#              }
#           }
# --------------

      models[[i]]$Statistics$StandardErrorInfo <- 'SE by coeftest()'
#        rownames( models[[j]][[1]] ) <- rownames(CoefTests[[j]])
#        if ('bootstrap' %in% class(CoefTests[[j]]))
#          models[[j]][[1]][,4]  <- -1
      }

    if (!is.null(names(AllModels)) && length(names(AllModels))>=i) models[[i]]$ModelName <- names(AllModels)[i]
    if (!is.null(ModelNames) && length(ModelNames)>=i)  models[[i]]$ModelName <- ModelNames[i]
    if (is.null(models[[i]]$ModelName) ||  models[[i]]$ModelName == '')  models[[i]]$ModelName <- as.character(paste('\'(',i, ')', sep=''))
}


# Prepare the CoefNames list
#   For convenience, user are allowed to use shortcuts when they fill out the list.
#   Therefore a cleanup of the list is required before internal use.

if (!is.null(Factors)  && !is.list(Factors)) Factors <- as.list(Factors)
if (length(Factors)>0)
for (i in 1:length(Factors))
{
  if(is.null(names(Factors)[i]) || is.na(names(Factors)[i]) || names(Factors)[i] =='' ) names(Factors)[i] <- Factors[[i]][[1]]
}

#CoefNames <- '<Generic>'
CoefNames <- PrepareCoefNames()

# Some models store fixed effects separately from other regression coefficients.
# Variables added to the model$FixedEffects list are added to the Factors
PreDefFE <- unique(unlist(lapply(models, function(x) x['FixedEffects'])))
if (!is.null(PreDefFE)) PreDefFE <- PreDefFE[!is.na(PreDefFE)]
if (length(PreDefFE)>0)
{
  for (FEi in PreDefFE)
  {
    #FEi <- PreDefFE[2]
    if (!GetPrettyCoefList(paste0('FE>>',FEi))[[1]] %in% names(Factors))
     {
       Factors <- append(Factors, setNames(paste0('(^', FEi, '.)|(^factor\\(', FEi, '\\).)|(^FE>>',FEi ,'$)'), FEi))
       CoefNames <- append(CoefNames, setNames(paste0('(^', FEi, '.)|(^factor\\(', FEi, '\\).)|(^FE>>',FEi ,'$)'), FEi))
     }
  }
}


# Adds all variables to the list which are not explicitely mentioned by the user

# It starts with the most complete model.
#    The idea behind the ordering is that complete models are usually at the end of the table.

#MostCompleteOrder <- order( unlist(lapply(models, function(x) nrow(x[[1]]))), decreasing = T)
MostCompleteOrder <- order( unlist(lapply(models, function(x) nrow(x$CoefValues) )), decreasing = T)

for (ModI in MostCompleteOrder)
{
  #ModI <- 1     ModI <- 4

  ModCoefs <- GetPrettyCoefList(rownames(models[[ModI]]$CoefValues))
  #ModCoefs <- GetPrettyCoefList( rownames(models[[ModI]][[1]])  )
  models[[ModI]]$CoefNames <- ModCoefs
  for (j in 1:length(ModCoefs))
  {
    #j <- 1
    ModCoefJ <- names(ModCoefs)[j]
    if (any(ModCoefJ == names(CoefNames))) next
    GenPos <- grep('<Generic>',   names(CoefNames), perl = T)
    if (length(GenPos)>0)
    for (Generics in GenPos)
    {
      # Generics <- 1
      if (length(grep(CoefNames[[Generics]], ModCoefJ, perl = T)) > 0)
      {
        CoefNames <- append(CoefNames, setNames(as.list(list(paste('^',ModCoefJ, '$', sep =''))), ModCoefJ), Generics - 1)
        break
      }
    }
  }
}

if ( length(GenPos <- grep('<Generic>',   names(CoefNames), perl = T)) >0 )
  CoefNames <- CoefNames[-grep('<Generic>',   names(CoefNames), perl = T)]

#Test <<- models

 Dubl <- duplicated(names(CoefNames))
 if (any(Dubl))
 {
   stop(paste('Some coefficents are provided more then once: ', paste(names(CoefNames)[Dubl], collapse = ', '), sep = '' ))
 }



# Create the raw array that can be passed to functions that render the output in
# a specific format.

  #---- Coefficients (i.e., fixed effects in mixed model terminology)

  # AllModelsAndIndexes <- lapply(lapply....
  # This is a list that contains sub-lists of coefficient names for all models
  # list  --| - Model 1: list(with all coefficient names)
  #         |
  #         | - Model x:  list(with all coefficient names)
  #         |
  #         |
  #         | - Model n: list(with all coefficient names)
  #
  #
  # The lists of coefficent names of each model are passed to GetOutputCellIndexesFromModel as ModelCoefsJ argument.
  # GetOutputCellIndexesFromModel creates a list of all coeficients that shall be printed and checks for each element of
  # the list whether the element, i.e. the coefficient name, is present in the list ModelCoefsJ, i.e. the list of coefficients
  # from model j.
  #


  GetOutputCellIndexFromModel <- function(CoefI, ModelCoefsJ)
  {
    if (length(locations <- which(CoefI == ModelCoefsJ))==0)
      return(ifelse((CoefI %in% names(Factors)), -1, -2))

    if (length(locations) == 1 && CoefI %in% names(Factors))
      return(0)

    if (length(locations) > 1)
    {
      if (!CoefI %in% names(Factors))
        warning(CoefI, ' has multiple matches at model ', j, ', although it is not defined as factor.')
        return(0)
    }
    return(locations)
  }
  # GetOutputCellIndexFromModel returns: -2: Coefficient not found.                             => blank
  #                                      -1: Coefficient not found and coefficient is a factor. => no
  #                                       0: Coefficient found and coefficient is a factor.     => yes
  #                                      >0: Coefficient found.

  GetOutputCellIndexesFromModel <- function(ModelCoefsJ) unlist(lapply(names(CoefNames), GetOutputCellIndexFromModel, ModelCoefsJ = ModelCoefsJ))
  AllModelsAndIndexes <- lapply(lapply(models, function(x) names(x$CoefNames)), GetOutputCellIndexesFromModel)
    # GetOutputCellIndexesFromModel(ModelCoefsJ =  lapply(models, function(x) names(x$CoefNames))[[3]] )

  EmptyVals <-  array(c(rep(NA, 3 *5), -3,-2,-1), dim =c(3,1,6),  dimnames=list(1:3, 'a', c('Coefficient', 'SE', 'SomeStat', 'StatIndicator', 'PValue', 'SpecialText')))
  RawCoefArray <- mapply( function(model, indexes) abind(EmptyVals,  model$CoefValues, along = 1)[indexes+3,,, drop = F]  , model = models,  indexes = AllModelsAndIndexes, SIMPLIFY = F)
  RawCoefArray <- do.call(abind, args =list(RawCoefArray, along= 2))
  RawCoefArray <-    unlist(RawCoefArray)
  #abind(EmptyVals,  models[[1]]$CoefValues, along = 1)  [indexes+3,,, drop = F]



  #---- To do: Other stuff (e.g., random effects)

  # Bottom matter


  if (!is.null(BottomMatter) && !is.na(BottomMatter) && length(BottomMatter)>0)
  {
    # Collect all relevant names (internal, external etc.) for the variables supplied by the user.
    BottomMatter <- unique(BottomMatter)
    MInfo <- ModelInfo[ModelInfo$External %in% BottomMatter ,]

    # This function makes sure that the corresponding value for each model is return or an empty text.
    # The empty text is important. Otherwise variables which are not available for a certain model
    # would be lost when applying commans like unlist.
    ExtractModelInfoOrEmptyString <- function(x, MInf)
    {
      #MInf <- MInfo$Internal
      x$Statistics[MInfo$Internal]
      tmp <- lapply(MInf, function(y)  setNames(ifelse(is.null(unlist(x$Statistics[y])), '', x$Statistics[y]), y))
      return(data.frame(unlist(tmp)))
    }
    RawBottomMatter <- lapply(models, ExtractModelInfoOrEmptyString, MInf =MInfo$Internal)
    RawBottomMatter <- do.call(cbind, RawBottomMatter)
    DoAppendBottomMatter <- nrow(RawBottomMatter)>0
  }

# Format the output

# Text:
  # Coefficients

  FormattedCoefficients.Text <- apply(RawCoefArray, c(1, 2), FormatCoefficientOutput.Text, PrintStat = 'SE', Digits = Decimals, CoefficientOnly = CoefOnly)
  rownames(FormattedCoefficients.Text) <- names(CoefNames)

  if (ShowDependentVariable)
  {
    DepVars <- unlist(lapply(models, function(x) ifelse(is.null(x$DependentVariable), '', x$DependentVariable)))
    # To do: check DependentVarPretty for nicer names
    FormattedCoefficients.Text <- rbind('Dependent variable' = DepVars, FormattedCoefficients.Text)
  }

  if (ShowModelType)
  {
    ModTypes <- unlist(lapply(models, function(x) ifelse(is.null(x$ModelType), '', x$ModelType)))
    # To do: check DependentVarPretty for nicer names
    FormattedCoefficients.Text <- rbind('Model' = ModTypes, FormattedCoefficients.Text)
  }
  colnames(FormattedCoefficients.Text) <- lapply(models, '[[', 'ModelName')


  # BottomMatter
  if (DoAppendBottomMatter)
  {
    rownames(RawBottomMatter) <- MInfo$Output.Text

    # To do: Some variables should have a specific formatting E.g., the number of obervations is always an integer and doesn't need
    #        decimals. Ohter fields should represent non-atomic values.

    Formatted.BottomMatter.Text <- apply(RawBottomMatter, c(1,2), FormatCoefficientOutput.Text, Digits = Decimals)
    FormattedOutput.Text <- rbind(FormattedCoefficients.Text, Formatted.BottomMatter.Text)
  }else{FormattedOutput.Text <- FormattedCoefficients.Text}


# Other output types:

   # To do: LaTeX, HTML etc.



  if (WriteToClipboard)
    write.table(FormattedOutput.Text,"clipboard", sep="\t", col.names = NA , row.names = TRUE)

  return(FormattedOutput.Text)
}


