[Example]: https://github.com/fredcommo/IC50/blob/master/plot1.png
# IC50
[Example]

### Compute IC50 with 4-P/5-P logistic regressions, using optic densities (ODs) or proportion of control, and drug concentrations.

### Demo

```
# If not installed yet
install.packages("devtools")
install_github("rGithubClient", "brian_bot")
```

```
require("devtools")
require("rGithubClient")

getFilesList <- function(git, tag = ''){
  flist <- git@tree$path
  return(flist[grep(tag, flist)])
}

git <- getRepo('fredcommo/IC50')
Rlist <- getFilesList(git, '[^Demo].R')
sourceRepoFile(git, Rlist)
op <- par(no.readonly = TRUE)
ss <- function(i){set.seed(123345+i)}
```
```
# From ODs
x <- seq(log10(0.001), log10(10), len = 8)
y <- lapply(1:3, function(i){
  ss(i)
  rev(.LP5(0.1, 2.8, median(x), 1.5, 1, x) + rnorm(length(x), sd = .25))
  })
Resp <- do.call(c, y) + rnorm(3*8, .5, .1)
dose <- rep(signif(10^x, 2), 3)
dose <- log10(dose)
test <- IC50(dose, Resp, T0 = .5, Ctrl = max(Resp))
plot(test)
getEstimates(test)
```

```
# On the previous example, use 'method' to specify a model (Both, 4P, 5P).
  # if Both (default), 4P and 5P are compared, and the best one (goodness of fit) is returned.
test <- IC50(dose, Resp, T0 = .5, Ctrl = max(Resp), method = "4P")
plot(test)
getEstimates(test)
```

```
# From props of control
x <- seq(log10(0.001), log10(10), len = 8)
y <- lapply(1:3, function(i){
  ss(i)
  rev(.LP5(0, 1, median(x), 1.7, .25, x) + rnorm(length(x), .05, sd = .1))
})
Resp <- do.call(c, y)
dose <- rep(signif(10^x, 2), 3)
dose <- log10(dose)
test <- IC50(dose, Resp)
plot(test)
getEstimates(test)
```
