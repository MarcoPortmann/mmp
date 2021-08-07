#! multi-model print

# Author: Marco Portmann
# Last change: 07.08.2021
#
#



#' Get P-value
#'
#' @param Stats
#' @param StatsType
#'
#' @return
#' @keywords internal
#' @examples
GetPValue <- function(Stats, StatsType)
{
  if (is.null(Stats) || is.na(Stats) || is.null(StatsType) || is.na(StatsType)) return(NA)
  if (StatsType == -1)
  {
    return(2* (1- pnorm(abs(Stats))))
  }
  return(2* (1- pt(abs(Stats), StatsType )))
}


#' Get P-value from normal
#'
#' @param Stats
#'
#' @return
#' @keywords internal
#' @examples
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




#' Get stars text
#'
#' @param PValue
#' @param StarBreaks
#'
#' @return
#' @keywords internal
#' @examples
  GetStars.Text <- function(PValue, StarBreaks=c(-Inf, 0.01, 0.05, 0.1, Inf))
  {
    if (is.null(PValue) || is.na(PValue) || (!is.numeric(PValue) && !is.numeric(type.convert(PValue)))) return('')
    if (!is.numeric(PValue)) PValue <- type.convert(PValue)
    if (is.nan(PValue)) return('')
    return(paste(rep("*", max(4- as.numeric(cut(PValue, StarBreaks)), 0)), collapse=""))
  }



#' FormatNumber
#'
#' @param Num
#' @param Decimals
#' @param IsPercent
#' @keywords internal
#' @return
  #' @keywords internal
#' @examples
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


  # c('Coefficient', 'SE', 'SomeStat', 'StatIndicator', 'PValue', 'SpecialText')
#' FormatCoefficientOutput.Text
#'
#' @param Coefficient
#' @param SE
#' @param PValue
#' @param SpecialText
#' @param PrintStat
#' @param Digits
#' @param IsPercent
#' @param CoefficientOnly
#' @param SpecialTextList
#' @keywords internal
#' @return
#'
#' @examples
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
#' @return
#' @export
#' @import Hmisc lmtest abind
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
"Weights"               , "Weights"             , "Weights"           , "Weights"          ,
"NObsLR"               , "NObsLR"              , "N (l/r)"           , "N (l/r)"          ,
"EffNObsLR"             , "EffNObsLR"           , "Effective N (l/r)" , "Effective N (l/r)",
"Kernel"                , "Kernel"              , "Kernel"            , "Kernel"           ,
"BWEst"                 , "BWEst"               , "Bandwidths est."   , "Bandwidths est."  ,
"BWBias"                , "BWBias"              , "Bandwidths bias"   , "Bandwidths bias"  ,
"Rho"                   , "Rho"                 , "Rho"               , "Rho"
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
models <- lapply(models, ConvertNewToOldMPO)
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
      #MInfo <- x$Statistics[MInfo$Internal]
      CoerceLists <- function(x)
      {
        if (!is.list(x)) return(x)
        return(paste(lapply(x$x, FormatNumber), collapse = x$collapse))
      }
      Stats <- lapply(x$Statistics, CoerceLists)
      tmp <- lapply(MInf, function(y)  setNames(ifelse(is.null(unlist(Stats[y])), '', Stats[y]), y))
      return(data.frame(unlist(tmp)))
    }
    #browser()
    RawBottomMatter <- lapply(models, ExtractModelInfoOrEmptyString, MInf = MInfo$Internal)
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


