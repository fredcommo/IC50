# # Define class
setClass('cellResp', representation(dose = 'vector',
                                    Resp = 'vector',
                                    survProp = 'vector',
                                    initPar = 'vector',
                                    LPweight = 'numeric',
                                    yfit = 'vector',
                                    xCurve = 'vector',
                                    yCurve = 'vector',
                                    model = 'ANY',
                                    goodness = 'ANY',
                                    estimates = 'data.frame',
                                    LP4 = 'ANY',
                                    LP5 = 'ANY'))

# # Constructor
cellResp = function(dose = dose, Resp = Resp, survProp = NA, initPar = NA, LPweight = 0, yfit = NA,
                    xCurve = NA, yCurve = NA, model = NA, goodness = NA, estimates = data.frame(),
                    LP4 = ggplot(), LP5 = ggplot()){
  new('cellResp', dose = dose, Resp = Resp, survProp = survProp, initPar = initPar, LPweight = LPweight,
      yfit = yfit, xCurve = xCurve, yCurve = yCurve, model = model, goodness = goodness, estimates = estimates,
      LP4 = LP4, LP5 = LP5)
}

# # setGenerics
setGeneric(".getDose", function(object) standardGeneric(".getDose"))
setGeneric(".getDO", function(object) standardGeneric(".getDO"))
setGeneric("getSurvProp", function(object) standardGeneric("getSurvProp"))
setGeneric(".getInitPar", function(object) standardGeneric(".getInitPar"))
#setGeneric("getFitValues", function(object) standardGeneric("getFitValues"))
setGeneric("getXcurve", function(object) standardGeneric("getXcurve"))
setGeneric("getYcurve", function(object) standardGeneric("getYcurve"))
setGeneric("getModel", function(object) standardGeneric("getModel"))
setGeneric("getParam", function(object) standardGeneric("getParam"))
setGeneric("getGoodness", function(object) standardGeneric("getGoodness"))
setGeneric("getEstimates", function(object) standardGeneric("getEstimates"))
setGeneric("getLP4", function(object) standardGeneric("getLP4"))
setGeneric("getLP5", function(object) standardGeneric("getLP5"))

# # Methods
setMethod(".getDose", "cellResp", function(object) return(object@dose))
setMethod(".getDO", "cellResp", function(object) return(object@Resp))
setMethod("getSurvProp", "cellResp", function(object) return(object@survProp))
setMethod(".getInitPar", "cellResp", function(object) return(object@initPar))
#setMethod("getFitValues", "cellResp", function(object) return(object@yfit))
setMethod("getXcurve", "cellResp", function(object) return(object@xCurve))
setMethod("getYcurve", "cellResp", function(object) return(object@yCurve))
setMethod("getModel", "cellResp", function(object) return(object@model))
setMethod("getParam", "cellResp", function(object){
  model <- getModel(object)
  bottom <- model$estimate[1]
  top <- model$estimate[2]
  xmid <- model$estimate[3]
  scal <- model$estimate[4]
  s <- model$estimate[5]
  return(cbind.data.frame(bottom = bottom, top = top, xmid = xmid, scal = scal, s = s))
})
setMethod('getGoodness', 'cellResp', function(object) return(object@goodness))
setMethod('getEstimates', 'cellResp', function(object){
  estim <- object@estimates
  return(estim[order(estim$Surv, decreasing = TRUE),])
  })

setMethod("getLP4", "cellResp", function(object) return(object@LP4))
setMethod("getLP5", "cellResp", function(object) return(object@LP5))
setMethod("plot", signature = "cellResp",
          function(object, x=NA, y=NA, pcol = 'grey50', lcol = 'grey25', cex = 1.5,
                   AddIC = .5, AddSd = TRUE, unit = 'µM', Title = NA, xlab = 'Log10(Drug[c])',...){
            op <- par(no.readonly = TRUE)
            par(las = 1, cex.axis = 1.5, cex.lab = 1.75, mar = c(6.5, 5.5, 4, 2), mgp = c(3.5, 1, 0))
            dose <- .getDose(object)
            survProp <- getSurvProp(object)
            newx <- getXcurve(object)
            newy <- getYcurve(object)
          #  my <- sapply(unique(dose), function(d) {mean(survProp[dose == d], na.rm = TRUE)})
            my <- sapply(unique(dose), function(d) {props <- survProp[dose == d]; prod(props)^(1/length(props))})
            mx <- unique(dose)
            r2adj <- round(summary(getGoodness(object))$adj.r.squared, 3)
            plot(mx, my, col = pcol, cex = cex, ylim = range(0, 1.1), xlab = xlab, ylab = 'Survival',...)
            points(mx, my, pch = 1, cex = cex)
            legend('topright', legend = paste('Goodness of fit:', r2adj), bty = 'n', cex = 1.5)
            
            if(!is.na(AddIC)){
              sigma <- summary(getGoodness(object))$sigma
              estim <- .estimateRange(AddIC, sigma, getParam(object), B = 1e4)
              legend1 <- sprintf("IC%d : %.2f%s", AddIC*100, estim[2], unit)
              legend2 <- sprintf("[%.2f, %.2f]", estim[1], estim[3])
              legend('bottomleft', legend = c(legend1, legend2), cex = 1.5, text.col = 'steelblue4', bty = 'n')
            }
            
            if(AddSd){
              bounds <- .IClm(getGoodness(object), getYcurve(object))
              xx <- c(newx, rev(newx))
              yy <- c(bounds$lo, rev(bounds$hi))
              polygon(xx, yy, border = NA, col = rgb(.8,.8,.8, .3))
              
              Sd <- sapply(unique(dose), function(d) {sd(survProp[dose == d], na.rm = TRUE)})
#               Sd <- sapply(unique(dose), function(d) {
#                 props <- survProp[dose == d]
#                 n <- length(props)
#                 p <- prod(props)^(1/n)
#                 sqrt(p*(1-p)/n)
#                 })
              pas <- (max(mx)-min(mx))/(length(mx)-1)/10
              lapply(1:length(Sd), function(i){
                segments(x0 = mx[i], x1 = mx[i], y0 = my[i]-Sd[i], y1 = my[i]+Sd[i], lty = 2, lwd = 2)
                segments(x0 = mx[i]-pas, x1 = mx[i]+pas, y0 = my[i]-Sd[i], y1 = my[i]-Sd[i], lty = 2, lwd = 2)
                segments(x0 = mx[i]-pas, x1 = mx[i]+pas, y0 = my[i]+Sd[i], y1 = my[i]+Sd[i], lty = 2, lwd = 2)      
              })
            }
            
            lines(newy ~ newx, col = lcol,...)
            if(object@LPweight != 0)
              Sub = "Weighted 5P logistic regr. (DoseResp package, version v.0)"
            else Sub = "Non weighted 5P logistic regr. (DoseResp package, version v.0)"
            title (main = Title, sub = Sub, cex.sub = .75)
            par(op)
          }
)

setMethod('show', signature = 'cellResp',
          function(object){return(getEstimates(object))}
)
