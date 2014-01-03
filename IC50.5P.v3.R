IC50 <- function(dose, Resp, T0 = NA, Ctrl = NA, Prop = FALSE, LPweight = 0.25, fixB = NA, fixT = NA, fixS = NA,
                    method = c("Both", "4P", "5P"), showIC = .5, showSd = TRUE, output = TRUE,...)
  
  {
  method <- match.arg(method)
  if(method == "Both") cat("Running 4P vs. 5P...\n")
  if(method == "4P") cat("Running 4P...\n")
  if(method == "5P") cat("Running 5P...\n")

  if(any(is.na(dose))){
    Resp <- Resp[!is.na(dose)]
    dose <- dose[!is.na(dose)]
  }
  if(any(is.na(Resp))){
    dose <- dose[!is.na(Resp)]
    Resp <- Resp[!is.na(Resp)]
  }
  
  object <- cellResp(dose = dose, Resp = Resp, LPweight = LPweight)
  if(Prop) object@survProp <- Resp
  else object@survProp <- .survProp(Resp, T0, Ctrl)
	# Optimisation step using nlm()
      init <- .initPar(.getDose(object), getSurvProp(object))
      if(method %in% c("Both", "4P")){
        model4 <- nlm(f = .sce, p = init,
                    x = .getDose(object), yobs = getSurvProp(object),
                    Weights = rep(1, length(resp)), LPweight = LPweight,
                    fixB = fixB, fixT = fixT, fixS = 1)
        object@model <- model4
        object@goodness <- .fit(model4, .getDose(object), getSurvProp(object))
        Param <- .getPar(model4)
      }

      if(method %in% c("Both", "5P")){
        model5 <- nlm(f = .sce, p = init,
                    x = .getDose(object), yobs = getSurvProp(object),
                    Weights = rep(1, length(resp)), LPweight = LPweight,
                    fixB = fixB, fixT = fixT, fixS = fixS)
        object@model <- model5
        object@goodness <- .fit(model5, .getDose(object), getSurvProp(object))
        Param <- .getPar(model5)
        }
  
	# Get best model
    if(method == "Both"){
      bestModel <- .getBestModel(object, model4, model5)
      object@goodness <- bestModel$goodness
      object@model <- bestModel$model
      Param <- bestModel$param
      }
  
	# Estimate critical points
		  object@xCurve <- seq(min(dose, na.rm = TRUE)*1.1, max(dose, na.rm = TRUE)*1.1, length = 100)
      object@yCurve <- .LP5(Param$bottom, Param$top, Param$xmid, Param$scal, Param$s, getXcurve(object))

	    # Inflexion point using the second derivative
		    Xflex = Param$xmid + (1/Param$scal)*log10(Param$s)
				Yflex = Param$bottom + (Param$top - Param$bottom)*(Param$s/(Param$s+1))^Param$s

			# Compute simulations to estimate the IC50 conf. interval
        targets <- seq(.1, .9, by = .1)
        sigma <- summary(getGoodness(object))$sigma
        estimates <- lapply(targets, function(target){.estimateRange(target, sigma, Param = Param, B = 1e4)})
        estimates <- cbind.data.frame(Resp = targets, do.call(rbind, estimates))
        colnames(estimates) <- c('Surv', 'Dmin', 'D', 'Dmax')
        object@estimates <- estimates

    # Compute Area
        x <- getXcurve(object)
        y <- getYcurve(object)
        object@AUC <- data.frame(trapezoide = AUC(x, y), Simpson = Simpson(x, y))

  return(object)
}
