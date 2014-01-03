############################
# Required packages
if(!'ggplot2' %in% installed.packages()){
  cat('Installing ggplot2...\n')
  install.packages('ggplot2')
  }
require(ggplot2)
require(stats)

# 5P logistic functions
.survProp <- function(y, T0 = NA, Ctrl = NA){
  if(is.na(Ctrl)) Ctrl <- max(y, na.rm = TRUE)
  if(is.na(T0))
    return(y/Ctrl)
  else return((y-T0)/(Ctrl-T0))
}

.scaleResp <- function(y){
  miny <- min(y, na.rm = TRUE)
  maxy <- max(y, na.rm = TRUE)
  y <- (y - miny)/(maxy - miny)
}

.LP5 <- function(bottom, top, xmid, scal, s,  x){
  fit <-bottom+(top-bottom)/(1+10^((xmid-x)*scal))^s
  return(fit)
}

# Weighted SCE Function (Sum of Squared errors)
.sce <- function(param, x, yobs, Weights, LPweight, fixB, fixT, fixS){
  bottom <- param[1]
  top <- param[2]
  xmid <- param[3]
  scal <- param[4]
  s <- param[5]
  if(!is.na(fixB)) bottom = fixB
  if(!is.na(fixT)) top = fixT
  if(!is.na(fixS)) s = fixS
  ytheo <- .LP5(bottom, top, xmid, scal, s, x)
  residus <- yobs - ytheo
  Weights <- (1/(residus^2))^(LPweight)
  return(sum(Weights*(yobs - ytheo)^2))
}

# Get weights (not used)
.sce.diag <- function(yobs, ytheo, w) {
  sq.res <- (yobs - ytheo)^2
  weights <- 1/sq.res^w
  return(weights)
}

.initPar <- function(x, y){
  bottom.ini <- min(y, na.rm = T)
  top.ini <- max(y, na.rm = T)
  xmid.ini = (max(x, na.rm = T) + min(x, na.rm = T))/2
  z <- (y - bottom.ini)/(top.ini - bottom.ini)
  z[z<=0] <- 0.01; z[z>=1] <- 0.99
  scal.ini = coef(lm(x ~ log(z/(1-z))))[2]
  scal.ini <- as.numeric(scal.ini)
  s.ini = 1
  return(c(bottom.ini, top.ini, xmid.ini, scal.ini, s.ini))
}

.getPar <- function(LPmodel){
  bottom <- LPmodel$estimate[1]
  top <- LPmodel$estimate[2]
  xmid <- LPmodel$estimate[3]
  scal <- LPmodel$estimate[4]
  s <- LPmodel$estimate[5]
  return(cbind.data.frame(bottom = bottom, top = top, xmid = xmid, scal = scal, s = s))
}

.fit <- function(model, dose, yobs){
  Par <- .getPar(model)
  yfit <- .LP5(Par$bottom, Par$top, Par$xmid, Par$scal, Par$s, dose)
  lmLP <- lm(yfit ~ yobs)  #, weights = weights)
  return(lmLP)
}
  
.getBestModel <- function(object, model4, model5){
  yobs <- getSurvProp(object)

  fit4 <- .fit(model4, .getDose(object), yobs)
  r4 <- summary(fit4)$adj.r.squared
  fit5 <- .fit(model5, .getDose(object), yobs)
  r5 <- summary(fit5)$adj.r.squared
  
  if(r4 > r5){
    cat('The 4-parameters model looks good!\n')
    return(list(model = model4, param = .getPar(model4), goodness = fit4))
  }
  else{
    cat('The 5-parameters model looks good!\n')
    return(list(model = model5, param = .getPar(model5), goodness = fit5))
  }
}

.IClm <- function(lmModel, newy){
  res <- lmModel$residuals
  Sqres <- sum(res^2)
  yobs <- lmModel$model$yobs
  yfit <- lmModel$model$yfit
  n <- length(yobs)
  ybar <- mean(yobs, na.rm = TRUE)
  sb <- sqrt(1/(n-2)*sum(res^2)/sum((yobs-ybar)^2))
  t <- qt(.975, n-2)
  IC <- t*sqrt(1/(n-2)*Sqres*(1/n+(newy - ybar)^2/sum((newy - ybar)^2)))
  lo <- newy - IC
  hi <- newy + IC
  return(list(lo = lo, hi = hi))
}
  
.invModel <- function(Param, target){
  if(any(target>=Param$top))
    target[target>=Param$top] <- Param$top*.99
  if(any(target<Param$bottom))
    target[target<=Param$bottom] <- Param$bottom*1.01
    return(Param$xmid - 1/Param$scal*log10(((Param$top - Param$bottom)/(target - Param$bottom))^(1/Param$s)-1))
}

.estimateRange <- function(target, sigma, Param, B = 1e4){
#  Ytarget = Param$bottom + (Param$top - Param$bottom)*target
  Xtarget = .invModel(Param, target)
  if(is.na(Xtarget)) minD <- D <- maxD <- NA
  else{
    Ytmp <- target + rnorm(B, 0, sigma)
    estimate <- .invModel(Param, Ytmp)
    Q <- quantile(estimate, probs=c(.025, .5, .975), na.rm=T)
    Dmin <- signif(10^Q[1], 2)
    D <- signif(10^Q[2], 2)
    Dmax <- signif(10^Q[3], 2)
  }
  return(as.numeric(c(Dmin, D, Dmax)))
}

AUC <- function(x, y){
  auc <- lapply(2:length(x), function(i){
    da <- x[i]-x[i-1]
    db <- y[i]-y[i-1]
    y[i]*da +1/2*db*da
  })
  return(do.call(sum, auc))
}

# Simpson method
Simpson <- function(x, y){
  dx <- mean(diff(x, lag = 1), na.rm = TRUE)
  n <- length(y)
  if(n%%2 != 0){
    x <- x[-n]
    y <- y[-n]
    n <- length(x)
  }
  f1 <- y[1]
  fn <- y[n]
  fy <- y[2:(n-1)]*rep(c(4, 2), (n-2)/2)
  return(dx/3*(f1 + sum(fy) + fn))
}

# PlotResp <- function(dose, resp, estimates, newX, newY, pcol, lcol, Title, unit, showIC, showSd,...){
#   my <- sapply(unique(dose), function(d) {mean(resp[dose == d], na.rm = TRUE)})
#   mx <- unique(dose)
#   plot(my ~ mx, col = pcol, ylim = range(0, 1.1), ylab = 'Survival',...)
#   
#   if(!is.na(showIC)){
#     legend1 <- sprintf("IC%d : %.2f%s", showIC*100, estimates$D[estimates$Surv == showIC], unit)
#     legend2 <- sprintf("[%.2f, %.2f]", estimates$Dmin[estimates$Surv == showIC], estimates$Dmax[estimates$Surv == showIC])
#     legend('bottomleft', legend = c(legend1, legend2), cex = 1.5, text.col = 'steelblue4', bty = 'n')
#   }
#   
#   if(showSd){
#     Sd <- sapply(unique(dose), function(d) {sd(resp[dose == d], na.rm = TRUE)})
#     pas <- (max(mx)-min(mx))/(length(mx)-1)/5
#     lapply(1:length(Sd), function(i){
#       segments(x0 = mx[i], x1 = mx[i], y0 = my[i]-Sd[i], y1 = my[i]+Sd[i], lty = 2, lwd = 2)
#       segments(x0 = mx[i]-pas, x1 = mx[i]+pas, y0 = my[i]-Sd[i], y1 = my[i]-Sd[i], lty = 2, lwd = 2)
#       segments(x0 = mx[i]-pas, x1 = mx[i]+pas, y0 = my[i]+Sd[i], y1 = my[i]+Sd[i], lty = 2, lwd = 2)      
#     })
#   }
#   
#   lines(newY ~ newX, col = lcol,...)
#   Sub = "Weighted 5P logistic regr. (DoseResp package, version v.0)"
#   #if(LPweight==0) Sub = "Non weighted 5P logistic regr. (DoseResp package, version v.0)"
#   title (main = Title, sub = Sub, cex.sub = .75)
# }
