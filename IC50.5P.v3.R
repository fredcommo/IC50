IC50.5P <- function(dose, Resp, T0 = NA, Ctrl = NA, LPweight = 0.25, fixB = NA, fixT = NA, fixS = NA,
                    method = c("Both", "4PL", "5PL"),
                    showIC = .5, showSd = TRUE, output = TRUE,...)
  {
  if(any(is.na(dose))){
    Resp <- Resp[!is.na(dose)]
    dose <- dose[!is.na(dose)]
  }
  if(any(is.na(Resp))){
    dose <- dose[!is.na(Resp)]
    Resp <- Resp[!is.na(Resp)]
  }
  
  object <- cellResp(dose = dose, Resp = Resp, LPweight = LPweight)
  object@survProp <- .survProp(Resp, T0, Ctrl)
	# Optimisation step using nlm()
      init <- .initPar(.getDose(object), getSurvProp(object))
#       method <- match.arg(method)
#       if(method %in% c("Both", "4PL")){
        model4 <- nlm(f = .sce.5P, p = init,
                    x = .getDose(object), yobs = getSurvProp(object),
                    Weights = rep(1, length(resp)), LPweight = LPweight,
                    fixB = fixB, fixT = fixT, fixS = 1)
#         bestModel$goodness <- .fit(model4, .getDose(object), getSurvProp(object))
#       }

#       if(method %in% c("Both", "5PL")){
        model5 <- nlm(f = .sce.5P, p = init,
                    x = .getDose(object), yobs = getSurvProp(object),
                    Weights = rep(1, length(resp)), LPweight = LPweight,
                    fixB = fixB, fixT = fixT, fixS = fixS)
#         bestModel$goodness <- .fit(model5, .getDose(object), getSurvProp(object))
#         }
  
	# Get best model
#     if(method == "Both"){
      bestModel <- .getBestModel(object, model4, model5)
      object@goodness <- bestModel$goodness
#       }
    object@model <- bestModel$model
    Param <- bestModel$param
  
	# Estimate critical points
		  object@xCurve <- seq(min(dose, na.rm = TRUE)*1.1, max(dose, na.rm = TRUE)*1.1, length = 100)
      object@yCurve <- .LP5(Param$bottom, Param$top, Param$xmid, Param$scal, Param$s, getXcurve(object))

	    # Inflexion point using the second derivative
		    Xflex = Param$xmid + (1/Param$scal)*log10(Param$s)
				Yflex = Param$bottom + (Param$top - Param$bottom)*(Param$s/(Param$s+1))^Param$s

			# Compute simulations to estimate the IC50 conf. interval
        P <- seq(.1, .9, by = .1)
        sigma <- summary(getGoodness(object))$sigma
        estimates <- lapply(P, function(p){.estimateRange(p, sigma, Param = Param, B = 1e4)})
        estimates <- cbind.data.frame(Resp = P, do.call(rbind, estimates))
        colnames(estimates) <- c('Surv', 'Dmin', 'D', 'Dmax')
        object@estimates <- estimates

  return(object)
}
