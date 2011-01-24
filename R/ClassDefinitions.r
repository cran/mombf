setClass("msPriorSpec", representation(priorType= "character", priorDistr= "character", priorPars = "vector"))

setValidity("msPriorSpec", function(object){
  msg <- NULL
  if (!any(object@priorType %in% c('coefficients','modelIndicator','nuisancePars'))) {
    msg <- "priorType must be 'coefficients', 'modelIndicator' or 'nuisancePars'"
  } else {
    if (object@priorType=='coefficients') {
      
      if (!any(object@priorDistr %in% c('pMOM','piMOM','peMOM'))) {
        msg <- "priorDistr must be 'pMOM','piMOM' or 'peMOM'"
      } else {
        if (object@priorDistr=='pMOM') {
          if (!all(c('tau','r') %in% names(object@priorPars))) msg <- "priorPars must contain elements named 'tau', 'r'"
        } else {
          if (!('tau' %in% names(object@priorPars))) msg <- "priorPars must contain an element named 'tau'"
        }
      }
      
    } else if (object@priorType=='modelIndicator') {

      if (!any(object@priorDistr %in% c('uniform','binomial'))) {
        msg <- "priorDistr must be 'uniform' or 'binomial'"
      } else {
        if (object@priorDistr=='binomial') { if (names(object@priorPars)!='p') msg <- "priorPars must contain an element named 'p'" }
      }
      
    } else {
      if (object@priorDistr!='invgamma') {
        msg <- "priorDistr must be invgamma"
      } else {
        if (!all(c('alpha','lambda') %in% names(object@priorPars))) msg <- "priorPars must contain elements named 'alpha', 'lambda'"
      }
    }
  }

  ifelse(is.null(msg),TRUE,msg)
}
)



            
