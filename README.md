# RobustSpline

R package *RobustSpline*. Thin-plate spline methods for functional data: 
1. regression for functional predictors and scalar response, and
2. robust (multivariate) functional location estimation.
   
The random functions are supposed to be discretely observed, but their domains are allowed to be multi-dimensional. Both robust and non-robust fits are implemented.

To install the package *RobustSpline* in Windows, you need to have *RTools* installed on your computer, see 

https://cran.r-project.org/bin/windows/Rtools/

Then it is enough to run

```R
#!R

library(devtools)
install_github("NagyStanislav/RobustSpline")

library(RobustSpline)
help(package="RobustSpline")
```

All major functions in the package are implemented both in R and in C++ for computational efficacy.

- Ioannis Kalogridis and Stanislav Nagy. (2025). [Robust functional regression with discretely sampled predictors.](https://arxiv.org/abs/2311.13291) _Computational Statistics and Data Analysis,_ to appear.
- Ioannis Kalogridis and Stanislav Nagy. (2025). Robust multidimensional location estimation from discretely sampled functional data.

## Regression for functional predictors and scalar response

Scalar-on-function regression with discretely observed functional data. Main functions are:

* `IRLS` for the general *Iteratively Reweighted Least Squares* procedure for a (possibly penalized and robust) estimation in a linear model;
* `ridge` for fast *Ridge Regression* with a given penalty matrix;
* `ts_reg` for fitting robust thin-plate splines regression;
* `ts_ridge` for fast fitting of non-robust thin-plate splines regression.

### Thin-plate spline regression: One-dimensional functional data

Generate the scalar-on-function regression data. The functional data are observed on an irregular domain of length `p.used`.

```R
errordist = function(x) rt(x,df=2)   # distribution of the errors: a t_2 distribution
Xdist = rnorm                        # distribution used for constructing the regressors

n = 500                         # sample size
snr = 0.2                       # signal-to-noise ratio
p = 100                         # length of the complete discretization grid
grid = seq(1/p, 1-1/p, len = p) # complete discretization grid
p.used = 50                     # number of observed points on the grid

# irregular grid where the random functions are observed
sub.p <- sort(sample(1:p, p.used)) 
sub.grid <- grid[sub.p]     

# generating the data
x <- matrix(0, n, p)
f1 <- -sin(5*grid/1.2)/0.5-1 # the true regression function beta0
for(i in 1:n){ # generate regressors X
  x[i, ] <- sqrt(2)*(1*pi-pi/2)^{-1}*Xdist(1)*
    sapply(grid, FUN= function(x) sin((1-1/2)*pi*x)  )
  for(j in 2:50){
    x[i, ] <- x[i, ] + (j*pi-pi/2)^{-1}*
      Xdist(1)*sqrt(2)*sapply(grid, FUN= function(x) sin((j-1/2)*pi*x))
  }
}
y0 <- x%*%f1/p                    # simulate response, without noise
y <- y0 + snr*sd(y0)*errordist(n) # add noise
tobs = matrix(sub.grid,ncol=1)    # observation points in the domain
```
For loss functions that are not power functions, a scale estimate is necessary to ensure approximate scale equivariance of the estimates.

```R
# preliminary scale estimator
lambda0 = exp(-15)
tspr = ts_preprocess(X = x[, sub.p],tobs,m=2)
Z = tspr$Z; H = tspr$H; 
res_sc = IRLS(Z,Y=y,lambda=lambda0,H,type="absolute") 
sc = RobStatTM::scaleM(res_sc$resids)
```

Fit thin-plate spline regression estimators, using both non-robust and robust methods.

```R
# non-robust estimator, square loss function
system.time(fit.ls <- ts_ridge(X = x[, sub.p], Y = y, tobs = tobs, 
                               m=2, jcv = "all"))
# > user  system elapsed 
# > 0.80    0.02    0.84

# robust estimator, Huber loss function
system.time(fit.huber <- ts_reg(X = x[, sub.p], Y = y, tobs = tobs, m=2,
                                type="Huber", jcv = "all", sc=sc))
# > user  system elapsed 
# > 0.84    0.03    0.89
```

Functions `ts_reg` and `ts_ridge` use cross-validation to find a value of the penalization parameter lambda. The criteria are always based on the residuals (`resids`) and hat values (`hats`) in the fitted models. Six different cross-validation methods (plus a custom method that can be provided to the functions) are implemented:
1. "AIC" Akaike's information criterion given by `mean(resids^2)+2*mean(hats)`, where `n` is the length of both `resids` and `hats`.
2. "GCV" Leave-one-out cross-validation criterion given by `mean((resids^2)/((1-hats)^2))`.
3. "GCV(tr)" Modified leave-one-out cross-validation criterion given by `mean((resids^2)/((1-mean(hats))^2))`.
4. "BIC" Bayes information criterion given by `mean(resids^2)+log(n)*mean(hats)`.
5. "rGCV" A robust version of GCV where mean is replaced by a robust M-estimator of scale of `resids/(1-hats)`, see [scaleTau2](https://search.r-project.org/CRAN/refmans/robustbase/html/scaleTau2.html) for details.
6. "rGCV(tr)" Modified version of a rGCV given by a robust M-estimator of scale of `resids/(1-mean(hats))`.

Compare the resulting estimates of the regression function beta. We take the cross-validation method `jcv=3` ("GCV(tr)").

```R
fit = fit.ls    # non-robust fit
jcv = 3         # cross-validation method selected

par(mfrow=c(1,2),mar=c(4,4,.5,.5))
plot(f1~grid, type="l", xlab = "Domain", ylab=expression(beta[0]), lwd=2)
lines(fit$beta_hat[,jcv]~tobs, col=2, lwd=3)
rug(tobs)
legend("topleft", c("true", "estimated"), lwd=2, col=1:2)
#
plot(y~fit$fitted[,jcv], xlab = "Fitted values", ylab="Responses")
points(y0~fit$fitted[,jcv], col="orange", pch=16)
legend("topleft", c("observed", "noiseless"), pch=c(1,16), col=c(1,"orange"))
```
The non-robust fit, square loss function. 
* Left: The true slope function (black) and its estimate (red),
* Right: The fitted values (x-axis) against the observed response values (y-axis, black) and the noiseless responses (y-axis, orange).

<img width="943" height="571" alt="fit ls" src="https://github.com/user-attachments/assets/c20d2ac5-664d-499d-a868-4a5941a1fd2b" />

The same plots with the robust fit, Huber loss function.

<img width="943" height="571" alt="fit huber" src="https://github.com/user-attachments/assets/8c46b3ad-d35c-4a8e-9c59-9479b0a87aa0" />

### Thin-plate spline regression: Multi-dimensional functional data

A real data example: data `ozone2` from R package [fields](https://www.rdocumentation.org/packages/fields/versions/16.3/topics/ozone2). The task is to predict the ozone concentration at a given station (southern Chicago), based on a sample of `n=89` bivariate discretely observed random functions (ozone concentration in the midwestern US over the period June 3,1987 through August 31, 1987, 89 days). 

Preprocessing and setting up the data.

```R
library(RobustSpline)
library(fields)
library(viridis)
library(akima)

data(ozone2)

n = nrow(ozone2$y)       # 89 days
p = nrow(ozone2$lon.lat) # 153 stations

### Missing data imputation: simply average of the station
y = ozone2$y
dim(y)
apply(y,1,function(x) sum(is.na(x)))
y = t(apply(y,1,function(x){ x[is.na(x)] = mean(x,na.rm=TRUE);
return(x)}))
apply(y,1,function(x) sum(is.na(x)))

### Data preprocessing
j = 45      # response will be the ozone concentration at j-th station
tun = 1.345 # tuning for the Huber estimator
Y = y[,j] 
tobs = ozone2$lon.lat[-j,]  # covariates observed at all but the j-th station
X = y[,-j]                  # covariates: all but the j-th station

dim(X)    # n-times-(p-1)
# > 89 152
dim(tobs) # (p-1)-times-2
# > 152   2
length(Y) # n
# > 89
```
Preprocessing for thin-plate regression.

```R
m = 2
tspr = ts_preprocess(X,tobs,m)
attach(tspr)
```

Huber estimator: Preliminary scale estimation with an M-scale of an undersmoothned L1 estimator.

```R
lambda0 = exp(-08)
res_sc = IRLS(Z,Y,lambda=lambda0,H,type="absolute",imax = 2000)
res_sc$last_check # convergence criterion
(sc = RobStatTM::scaleM(res_sc$resids))
```

Huber estimator of the regression function.

```R
type="Huber"
system.time({
  res.h <- ts_reg(X=X, Y=Y, tobs=tobs, m=2, type=type, sc=sc)
})

# >  user  system elapsed 
# > 257.41    1.73  280.93

cvj = 6  # scaletau(tr) as the cross-validation criterion to use
bhat.h = res.h$beta_hat[,cvj]
```
Least squares regression estimator of the regression function (square loss function, non-robust).

```R
### Least-squares estimator
type="square"
system.time({
  res.ls <- ts_reg(X=X, Y=Y, tobs=tobs, m=2, type=type, sc=sc)
})

# > user  system elapsed 
# > 1.30    0.00    1.42

cvj = 3 # GCV(tr) as the cross-validation criterion to use
bhat.ls = res.ls$beta_hat[,cvj]
```
Plot the resulting estimates. 

a) Least squares estimator.

```R
nlevels <- 1000
cols = terrain.colors(nlevels)
par(oma = c(0,0,0,0), mar = c(2.5, 2.5, 0.5, 0.5), cex = 1.5)

plot(ozone2$lon.lat,cex=.7,pch=16, type = "n")
s.ls <- interp(ozone2$lon.lat[-j, 1], ozone2$lon.lat[-j, 2], bhat.ls, 
               extrap = FALSE, xo = seq(min(ozone2$lon.lat[, 1]), 
                                        max(ozone2$lon.lat[, 1]),
                                        len = 400), 
               yo= seq(min(ozone2$lon.lat[, 2]), max(ozone2$lon.lat[, 2]), 
                       len = 400 ),linear = TRUE)

.filled.contour(x=s.ls$x, y= s.ls$y, z=s.ls$z, 
                levels = seq(pmin(min(bhat.h), min(bhat.ls)),
                             pmax(max(bhat.h), max(bhat.ls)), length = nlevels),
                col= cols)
US(add = TRUE, lwd = 2)
cityNames <- c("Chicago","Indianapolis","Milwaukee",
               "St Louis")
cityLonLatMat <- rbind(c(-87.6298,41.8781),
                       c(-86.1581,39.7694),
                       c(-87.9065,43.0389),
                       c(-90.1994,38.6270))
for (iCity in 1:length(cityNames))
{
  text(cityLonLatMat[iCity,1] + 0.15,cityLonLatMat[iCity,2],
       cityNames[iCity],adj = 0, cex = 0.7) }
points(ozone2$lon.lat, cex=.55, pch = 16 )
points(ozone2$lon.lat[j, 1], ozone2$lon.lat[j, 2], cex=0.9, pch=16, 
       col = "yellow")
```

<img width="1086" height="931" alt="LS1_new" src="https://github.com/user-attachments/assets/f9b413d7-7609-4e40-83b1-8a450c511314" />

b) Huber estimator.

```R
plot(ozone2$lon.lat,cex=.7, pch=16, type = "n")
s.r <- interp(ozone2$lon.lat[-j, 1], ozone2$lon.lat[-j, 2], bhat.h, 
              extrap = FALSE, 
              xo = seq(min(ozone2$lon.lat[, 1]), max(ozone2$lon.lat[, 1]),
                       len = 400), linear = TRUE,
              yo= seq(min(ozone2$lon.lat[, 2]), max(ozone2$lon.lat[, 2]), 
                      len = 400 ))

.filled.contour(x=s.r$x, y= s.r$y, z=s.r$z, 
                levels = seq(pmin(min(bhat.h), min(bhat.ls)),
                             pmax(max(bhat.h), max(bhat.ls)), 
                             length = nlevels), col= cols)
US(add = TRUE, lwd = 2)
cityNames <- c("Chicago","Indianapolis","Milwaukee",
               "St Louis")
cityLonLatMat <- rbind(c(-87.6298,41.8781),
                       c(-86.1581,39.7694),
                       c(-87.9065,43.0389),
                       c(-90.1994,38.6270))
for (iCity in 1:length(cityNames))
{
  text(cityLonLatMat[iCity,1] + 0.15,cityLonLatMat[iCity,2],
       cityNames[iCity],adj = 0, cex = 0.7) }
points(ozone2$lon.lat, cex=.55, pch=16)
# Outliers are depicted as red triangles
points(ozone2$lon.lat[j, 1], ozone2$lon.lat[j, 2], cex=0.9, pch=16, 
       col = "yellow")
```

<img width="1086" height="931" alt="H1_new0" src="https://github.com/user-attachments/assets/97dc63f3-7a00-4713-9032-30b878b40fd4" />



