# Demo
require('devtools')
require('rGithubClient')

getFilesList <- function(git, tag = ''){
  flist <- git@tree$path
  return(flist[grep(tag, flist)])
}

git <- getRepo('fredcommo/IC50')
Rlist <- getFilesList(git, '[^Demo].R')
sourceRepoFile(git, Rlist)

op <- par(no.readonly = TRUE)
ss <- function(i){set.seed(123345+i)}

# From ODs
x <- seq(log10(0.001), log10(10), len = 8)
y <- lapply(1:3, function(i){
  ss(i)
  rev(.LP5(0.1, 2.8, median(x), 1.5, 1, x) + rnorm(length(x), sd = .25))
  })
Resp <- do.call(c, y) + rnorm(3*8, .5, .1)
dose <- rep(signif(10^x, 2), 3)
dose <- log10(dose)
test <- IC50.5P(dose, Resp, T0 = .5, Ctrl = max(Resp))
plot(test)
getEstimates(test)

# From props of control
x <- seq(log10(0.001), log10(10), len = 8)
y <- lapply(1:3, function(i){
  ss(i)
  rev(.LP5(0, 1, median(x), 1.7, .25, x) + rnorm(length(x), .05, sd = .1))
})
Resp <- do.call(c, y)
dose <- rep(signif(10^x, 2), 3)
dose <- log10(dose)
test <- IC50.5P(dose, Resp)
plot(test)
getEstimates(test)
