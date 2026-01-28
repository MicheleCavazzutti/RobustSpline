#' Pre-processing raw data for thin-plate multivariate location estimation
#'
#' Given the discretely observed functional data as an input, pre-processes the 
#' dataset to a form suitable for fitting thin-plate spline functional location
#' estimation.
#'
#' @param Y Matrix of observed values of functional data \code{Y} of size 
#'  \code{p}-times-\code{n}, one column per functional observation, 
#'  rows corresponding to the positions in the rows of \code{tobs}. If matrix
#'  \code{Y} does not contain any missing values \code{NA}, it is also 
#'  possible to be supplied as a long vector of length \code{p*n}, stacked by
#'  columns of the matrix \code{Y}. If \code{Y} contains missing values 
#'  \code{NA}, it must be a matrix.
#'  
#' @param tobs Domain locations for the observed points of \code{X}. Matrix
#'  of size \code{p}-times-\code{d}, one row per domain point.
#' 
#' @param r Order of the thin-plate spline, positive integer.
#' 
#' @return A list of values:
#' \itemize{
#'  \item{"Y"}{ A vector of length \code{m} of all the reponses that are not 
#'  missing. If \code{Y} was a matrix without missing values, this is simply 
#'  \code{c(Y)}; if \code{Y} was a vector, this is directly \code{Y}. In case
#'  \code{Y} contained missing values, this is the elements of \code{c(Y)} that
#'  are not \code{NA}.}
#'  \item{"Z"}{ A matrix of size \code{m}-times-\code{p} for fitting the 
#'  penalized robust regression model for thin-plate coefficients.}
#'  \item{"H"}{ A penalty matrix of size \code{p}-times-\code{p} used for
#'  fitting the location estimation model with thin-plate splines.}
#'  \item{"Q"}{ A matrix containing the null space of the rows of Phi of size
#'  \code{p}-times-\code{p-M}, where \code{M} is the number of monomials
#'  used for the construction of the thin-plate spline. This matrix is used to
#'  pass from parameters gamma to parametrization chi.}
#'  \item{"Omega"}{ A matrix of size \code{p}-times-\code{p} containing the 
#'  \link{eta}-transformed matrix of inter-point distances in \code{tobs}.}
#'  \item{"Phi"}{ A matrix of size \code{p}-times-\code{M} corresponding to the
#'  monomial part of the thin-plate spline.}
#'  \item{"w"}{ A vector of weights of length \code{m}. Each weight corresponds
#'  to \code{1/(n*m[i])}, where \code{m[i]} is the number of observations for 
#'  the \code{i}-th function.}
#'  \item{"degs"}{ A matrix of size \code{M}-times-\code{d} with degrees of the
#'  monomials in each dimension per row. Used for the construction of 
#'  \code{Phi}.}
#'  \item{"tobs"}{ The same as the input parameter \code{tobs}, for later use.}
#'  \item{"M"}{ Number of all monomials used, is equal to 
#'  \code{choose(m+d-1,d)}.}
#'  \item{"r"}{ Order of the spline, positive integer.}
#'  \item{"p"}{ Number of observed time points, positive integer.}
#'  \item{"d"}{ Dimension of domain, positive integer.}
#'  \item{"n"}{ Sample size, positive integer.}
#' }
#'
#' @examples
#' d = 1  # dimension of the domain
#' m = 50 # number of observation times per function
#' tobs = matrix(runif(m*d), ncol=d) # matrix of observation times
#' n = 20 # sample size
#' truemeanf = function(x) 10+15*x[1]^2 # true mean function
#' truemean = apply(tobs,1,truemeanf) # discretized values of the true mean
#' Y = replicate(n, truemean + rnorm(m)) # matrix of discrete functional data, size p*n
#' tsp = ts_preprocess_location(Y, tobs, 2) # preprocessing matrices 
#' 
#' lambda0 = 1e-5 # regularization parameter
#' res_IRLS = IRLS(Z = tsp$Z, Y = tsp$Y, lambda = lambda0, H = tsp$H, type = "square")
#' res_ridge = ridge(Z = tsp$Z, Y = tsp$Y, lambda = lambda0, H = tsp$H)
#' # resulting estimates of the parameters theta, using IRLS and ridge
#' 
#' # testing that ridge and IRLS (square) give the same results
#' all.equal(res_IRLS$theta_hat, res_ridge$theta_hat)
#' 
#' res = res_ridge
#' resf = transform_theta_location(res$theta_hat, tsp)
#' 
#' plot(rep(tobs,n), c(Y), cex=.2, pch=16)
#' points(tobs, truemean, pch=16, col="orange")
#' points(tobs, resf$beta_hat, col=2, pch=16)
#' 
#' @importFrom MASS Null
#' @export

ts_preprocess_location = function(Y, tobs, r){
  
  # Y Observed values of functional data of size \code{p}-times-\code{n}
  # tobs Domain locations for the observed points of size \code{p}-times-\code{d}
  # r Order of the thin-plate spline, positive integer.
  
  p = nrow(tobs)    # number of distinct observation points
  d = ncol(tobs)    # dimension of the domain
  
  ### preparing the vector of effective values of the response
  if(sum(is.na(Y))==0){ 
    # if the design of observation is full, also vector Y is accepted
    if(is.matrix(Y)){
      if(nrow(Y)!=p) 
        stop("If Y is matrix, its number of rows must equal 
             the number of columns of tobs.")
      n = ncol(Y)
      Y = c(Y) # stack into a large vector of length p*n
    } else {
      n = length(Y)/p
      if(n%%1!=0) stop("The length of Y is not divisible by the number
                       of columns of tobs.")
    }
    w = rep(p,n) # vector of numbers of observations per function (length n)
    m = p*n      # total number of available observations
    inds = rep(1:p,n) # list of indices of active (observed) points, concatenated
  } else { # if the design of observation is not full, 
    if(!is.matrix(Y)) 
      stop("If Y contains NA's, only matrix form of Y is accepted.")
    
    if(nrow(Y)!=p) 
      stop("If Y is matrix, its number of rows must equal 
             the number of columns of tobs.")
    
    inds = apply(Y,2,function(x) which(!is.na(x))) 
    # list of indices of active (observed) points, one list per function
    inds = unlist(inds) # concatenated list of indices, length m
    w = apply(Y, 2, function(x) sum(!is.na(x))) 
    # vector of numbers of observations per function (length n)
    Y = Y[,w>0,drop=FALSE]
    # deleting functions that are never observed
    n = ncol(Y)
    m = sum(w) # total number of available observations
    Y = c(Y)   # stacking Y by columns
    Y = Y[!is.na(Y)] # vector of length m of clean observed values
  }
  
  M = choose(r+d-1,d)
  if(p-M<=0) stop(paste("p must be larger than",M))
  if(2*r<=d) stop(paste("r must be larger than",ceiling(d/2)))
  
  # Monomials phi of order <r
  allcom = expand.grid(replicate(d, 0:(r-1), simplify=FALSE))
  degs = as.matrix(allcom[rowSums(allcom)<r,,drop=FALSE])
  M = nrow(degs)
  if(M!=choose(r+d-1,d)) stop("Error in degrees of polynomials")
  
  # Fast matrix of Euclidean distances
  Em = matrix(sqrt(rowSums(
    apply(tobs,2,function(x) outer(x,x,"-")^2))),nrow=nrow(tobs))
  
  # Matrix Omega for the penalty term
  Omega = eta(Em,d,r)
  
  # Matrix Phi, monomials evaluated at tobs
  Phi = apply(degs,1,function(x) apply(t(tobs)^x,2,prod))
  # Phi = matrix(nrow=p,ncol=M)
  # for(i in 1:p) for(j in 1:M) Phi[i,j] = prod(tobs[i,]^degs[j,])
  
  # Null space of the rows of Phi p-(p-M)
  Q = MASS::Null(Phi)
  
  # Test that t(Phi)%*%Q = 0
  if(max(abs(t(Phi)%*%Q))>1e-4) 
    stop("Problem in computing the null space")
  if(ncol(Q)!=p-M) 
    stop("Problem in computing the null space")
  
  A = Omega%*%Q
  H = matrix(0, nrow=p, ncol=p)
  H[1:(p-M),1:(p-M)] = t(Q)%*%A
  ## NOTE: H is numerically not-symmetric (errors in the order 10^-10), but should be theoretically: therefore, I enforce symmetry
  H <- 0.5 * (H + t(H))
  Z = cbind(Omega%*%Q, Phi) # matrix p-times-p
  Z = Z[inds,] # matrix m-times-p
  # Z = do.call("rbind", replicate(n, Z, simplify = FALSE)) # 
  # replicates the matrix Z n-times, and then builds a final matrix
  
  # weights attached to observations (length m)
  w = rep(1/(w*n),w)
  
  return(list(Y=Y, Z=Z, H=H, Q=Q, 
              Omega=Omega, Phi=Phi, w=w, 
              degs=degs, tobs=tobs, M=M, r=r, p=p, d=d, n=n))          
}

#' Prediction for thin-plate multivariate location estimation
#'
#' Given the estimated thin-plate location spline from function 
#' \link{ts_location}, this function predicts the spline function values
#' at new observation points.
#'
#' @param tobs Domain locations for the observed points from which the location 
#' spline was estimated. Matrix of size \code{p1}-times-\code{d}, one row per 
#' domain point.
#' 
#' @param tobsnew Domain locations for the new points where the location 
#' spline is to be estimated. Matrix of size \code{p2}-times-\code{d}, one 
#' row per domain point.
#' 
#' @param gamma Vector \code{gamma} determining the first batch of parameters
#' of the thin-plate spline. Typically outcome of \code{gamma_hat} from 
#' functions \link{ts_location} or \link{transform_theta_location}.
#' 
#' @param delta Vector \code{delta} determining the second batch of parameters
#' of the thin-plate spline. Typically outcome of \code{delta_hat} from 
#' functions \link{ts_location} or \link{transform_theta_location}.
#' 
#' @param r Order of the thin-plate spline, positive integer.
#' 
#' @return A numerical vector of thin-plate spline values of length \code{p2}, 
#' corresponding to the rows of the matrix \code{tobsnew}.
#'
#' @examples
#' d = 1   # dimension of the domain
#' m = 10 # number of points per curve
#' tobs = matrix(runif(m*d), ncol=d)  # location of obsevation points
#' n = 500                            # sample size
#' truemeanf1 = function(x)   # true location function
#'   cos(4*pi*x[1])
#' truemean = apply(tobs,1,truemeanf1) # discretized values of the true location
#' Y = replicate(n, truemean + rnorm(m)) # a matrix of functional data, size m*n
#' 
#' # introduce NAs
#' obsprob = 0.5 # probability of observing a point
#' B = matrix(rbinom(n*m,1,obsprob),ncol=n)
#' for(i in 1:m) for(j in 1:n) if(B[i,j]==0) Y[i,j] = NA
#' 
#' # thin-plate spline fitting
#' res = ts_location(Y, tobs=tobs, r=2, type="square", method="ridge")
#' 
#' jcv = 3 # cross-validation criterion chosen
#' 
#' plot(rep(tobs,n), c(Y), cex=.2, pch=16, ann=FALSE)
#' title("True/estimated location function")
#' points(tobs, truemean, pch=16, col="orange")
#' fullt = matrix(seq(0,1,length=101),ncol=1)
#' lines(fullt, apply(fullt,1,truemeanf1), col="orange", lwd=2)
#' points(tobs, res$beta_hat[,jcv], col=2, pch=16)
#' 
#' # prediction
#' tobsnew = matrix(seq(0,1,length=501),ncol=1)
#' preds = ts_predict_location(tobs, tobsnew, res$gamma_hat[,jcv], 
#'   res$delta_hat[,jcv], r=2)
#' lines(tobsnew, preds, col=2, lwd=2)
#' legend("topleft",c("data","true","estimated"),
#'        pch=16, col=c(1,"orange",2))
#' 
#' @export

ts_predict_location = function(tobs, tobsnew, gamma, delta, r){
  
  # tobs Domain locations for the observed points of size \code{p}-times-\code{d}
  # r Order of the thin-plate spline, positive integer.
  
  p = nrow(tobs)    # number of distinct observation points
  d = ncol(tobs)    # dimension of the domain
  
  if(d!=ncol(tobsnew)) stop("Dimension mismatch in tobs and tobsnew.")
  if(length(gamma)!=p) stop("The length of gamma does not equal nrow(tobs).")
  
  M = choose(r+d-1,d)
  if(length(delta)!=M) stop("The length of delta does not match.")
  if(p-M<=0) stop(paste("p must be larger than",M))
  if(2*r<=d) stop(paste("r must be larger than",ceiling(d/2)))
  
  # Monomials phi of order <r
  allcom = expand.grid(replicate(d, 0:(r-1), simplify=FALSE))
  degs = as.matrix(allcom[rowSums(allcom)<r,,drop=FALSE])
  M = nrow(degs)
  if(M!=choose(r+d-1,d)) stop("Error in degrees of polynomials")
  
  # Fast matrix of Euclidean distances
  Em = distnAB(tobsnew, tobs, nrow(tobsnew), nrow(tobs), ncol(tobs))
  
  # Matrix Omega for the penalty term
  Omega = eta(Em,d,r)
  
  # Matrix Phi, monomials evaluated at tobs
  Phi = apply(degs,1,function(x) apply(t(tobsnew)^x,2,prod))
  
  return(Omega%*%gamma + Phi%*%delta)        
}

#' Transform the vector of estimated coefficients for location estimation
#'
#' Splits the vector of raw estimated coefficients (output of functions 
#' \link{IRLS}, \link{ridge} or \link{HuberQp}) after performing \link{ts_preprocess_locationi} 
#' into parts interpretable in the setup of thin-plate spline location 
#' estimation.
#'
#' @param theta Output vector of raw results of length \code{p} from function
#' \link{IRLS}, \link{ridge} or \link{HuberQp}.
#'
#' @param tspr Output of \link{ts_preprocess_location}.
#'
#' @return A list of estimated parameters:
#' \itemize{
#'  \item{"xi_hat"}{ Estimate of \code{xi}, the part of the parameters that
#'  correspond to matrix \code{Omega} when transformed by \code{Q}. A vector
#'  of length \code{(p-M)}.}
#'  \item{"delta_hat"}{ Estimate of \code{delta}, the part of the parameters that
#'  correspond to matrix \code{Phi}. A vector of length \code{M}, the number of
#'  monomials used for the construction of the thin-plate spline.}
#'  \item{"gamma_hat"}{ Estimate of \code{gamma}, the part of the parameters that
#'  correspond to matrix \code{Omega}. It holds true that 
#'  \code{gamma_hat = Q*xi_hat}. A vector of length \code{p}.}
#'  \item{"beta_hat"}{ Estimate of the location parameter \code{mu} 
#'  evaluated at the \code{p} points from \code{tobs}, where \code{Y} was
#'  observed.}
#' }
#'
#' @examples
#' d = 1  # dimension of the domain
#' m = 50 # number of observation times per function
#' tobs = matrix(runif(m*d), ncol=d) # matrix of observation times
#' n = 20 # sample size
#' truemeanf = function(x) 10+15*x[1]^2 # true mean function
#' truemean = apply(tobs,1,truemeanf) # discretized values of the true mean
#' Y = replicate(n, truemean + rnorm(m)) # matrix of discrete functiona data, size p*n
#' tsp = ts_preprocess_location(Y, tobs, 2) # preprocessing matrices 
#' 
#' lambda0 = 1e-5 # regularization parameter
#' res_IRLS = IRLS(Z = tsp$Z, Y = tsp$Y, lambda = lambda0, H = tsp$H, type = "square")
#' res_ridge = ridge(Z = tsp$Z, Y = tsp$Y, lambda = lambda0, H = tsp$H)
#' # resulting estimates of the parameters theta, using IRLS and ridge
#' 
#' # testing that ridge and IRLS (square) give the same results
#' all.equal(res_IRLS$theta_hat, res_ridge$theta_hat)
#' 
#' res = res_ridge
#' resf = transform_theta_location(res$theta_hat, tsp)
#' 
#' plot(rep(tobs,n), c(Y), cex=.2, pch=16)
#' points(tobs, truemean, pch=16, col="orange")
#' points(tobs, resf$beta_hat, col=2, pch=16)
#' 
#' @export

transform_theta_location = function(theta, tspr){
  p = tspr$p
  M = tspr$M
  Q = tspr$Q
  Omega = tspr$Omega
  Phi = tspr$Phi
  xihat=theta[1:(p-M)]
  dhat=theta[(p-M+1):p]
  ghat = Q%*%xihat
  return(list(xi_hat=xihat,
              delta_hat=dhat,
              gamma_hat=ghat,
              beta_hat = c(Omega%*%ghat + Phi%*%dhat)))
}

#' Robust thin-plate splines location estimation for functional data
#'
#' Provides a (potentially robust) thin-plates spline location estimator for 
#' discretely observed functional data. The functional data do not need to be 
#' observed on a common grid (that is, \code{NA} values in the matrix \code{Y}
#' below are allowed). The tuning parameter \code{lambda} is selected using a 
#' specified cross-validation criterion.
#'
#' @param Y Matrix of observed values of functional data \code{Y} of size 
#'  \code{p}-times-\code{n}, one column per observation, rows corresponding to 
#'  the positions in the rows of \code{tobs}. Functions observed on different 
#'  grids can be provided by including missing values \code{NA} in the matrix.
#'
#' @param tobs Domain locations for the observed points of \code{Y}. Matrix
#'  of size \code{p}-times-\code{d}, one row per domain point.
#' 
#' @param r Order of the thin-plate spline, positive integer.
#'
#' @param type The type of the loss function used in the minimization problem.
#' Accepted are \code{type="absolute"} for the absolute loss \code{rho(t)=|t|/2};
#' \code{type="quantile"} for the (asymmetric) quantile loss 
#' \code{rho(t)=t(alpha-I[t<0])} (\code{absolute} loss with \code{alpha=1/2});
#' \code{type="square"} for the square loss \code{rho(t)=t^2}; 
#' \code{type="Huber"} for the Huber loss \code{rho(t)=t^2/2} if 
#' \code{|t|<tuning} and \code{rho(t)=tuning*(|t|-tuning/2)} otherwise; and 
#' \code{type="logistic"} for the logistic loss 
#' \code{rho(t)=2*t + 4*log(1+exp(-t))-4*log(2)}.
#' 
#' @param alpha The order of the quantile if \code{type="quantile"}. By default
#' taken to be \code{alpha=1/2}, which gives the absolute loss 
#' (\code{type="absolute"}).
#' 
#' @param tuning A non-negative tuning constant for the Huber/quantile loss 
#' function (that is,  \code{type="Huber"}, \code{type="absolute"} or \code{type="quantile"}). 
#' If left to NULL is defaulted to standard values. See \link{IRLS} and \link{HuberQp} for further details.
#' 
#' @param jcv A numerical indicator of the cross-validation method used to 
#' select the tuning parameter \code{lambda}. The criteria are always 
#' based on the residuals (\code{resids}) and hat values (\code{hats}) in
#' the fitted models. Possible values are:
#' \itemize{
#'  \item{"all"}{ All the criteria below are considered.}
#'  \item{"AIC"}{ Akaike's information criterion given by 
#'  \code{mean(resids^2)+log(n)*mean(hats)}, where \code{n} is the length of
#'  both \code{resids} and \code{hats}.}
#'  \item{"GCV"}{ Leave-one-out cross-validation criterion given by
#'  \code{mean((resids^2)/((1-hats)^2))}.}
#'  \item{"GCV(tr)"}{ Modified leave-one-out cross-validation criterion 
#'  given by \code{mean((resids^2)/((1-mean(hats))^2))}.}
#'  \item{"BIC"}{ Bayes information criterion given by 
#'  \code{mean(resids^2)+2*mean(hats)}.}
#'  \item{"rGCV"}{ A robust version of \code{GCV} where mean is replaced
#'  by a robust M-estimator of scale of \code{resids/(1-hats)}, see 
#'  \link[robustbase]{scaleTau2} for details.}
#'  \item{"rGCV(tr)"}{ Modified version of a \code{rGCV} given by 
#'  a robust M-estimator of scale of \code{resids/(1-mean(hats))}.}
#'  \item{"custom"}{ The custom criterion given by function \code{custfun}. 
#'  Works only if \code{custfun} is part of the input.}
#'  }
#'  
#' @param vrs Version of the algorhitm to be used in function \link{ridge}; 
#' either \code{vrs="C"} for the \code{C++} version, or \code{vrs="R"} for the 
#' \code{R} version. Both should give (nearly) identical results, see 
#' \link{IRLS}.
#' 
#' @param method A method for estimating the fit. Possible options are 
#' \code{"IRLS"} for the IRLS algorithm, \code{"ridge"} for ridge regression, or
#' \code{"HuberQp"} for Huber regression based on quadratic programming.
#' Ridge is applicable only if \code{type="square"}; this method is much faster,
#' but provides only a non-robust fit. HuberQp is applicable only if \code{type="huber"};
#' this method is faster than the IRLS and also provides a robust fit. However, is not 
#' as fast as ridge.
#' 
#' @param plotCV Indicator of whether a plot of the evaluated cross-validation 
#' criteria as a function of \code{lambda} should be given.
#' 
#' @param lambda_grid An optional grid for select \code{lambda} from. By default
#' this is set to be an exponential of a grid of \code{lambda_length} 
#' equidistant values in the interval from -28 to -1.
#' 
#' @param lambda_length Number of elements in the grid of values \code{lambda}. 
#' By default chosen to be 51. 
#'
#' @param custfun A custom function combining the residuals \code{resids} and
#' the hat values \code{hats}. The result of the function must be numeric, see 
#' \link{GCV_crit}.
#' 
#' @return The output differs depending whether \code{jcv="all"} or 
#' not. If a specific cross-validation method is selected (that is, 
#' \code{jcv} is not \code{"all"}), a list is returned:
#'  \itemize{
#'  \item{"lambda"}{ The selected tuning parameter \code{lambda} that minimizes
#'  the chosen cross-validation criterion.}
#'  \item{"fitted"}{ A vector of fitted values using the tuning 
#'  parameter \code{lambda}. The length of this vector equals \code{m}, the 
#'  number of non-missing values in the matrix \code{Y}. If \code{Y} does not
#'  contain missing values, the length of this vector is \code{n*p}.}
#'  \item{"residuals"}{ A vector or a matrix of the same structure as 
#'  \code{fitted} with the residuals \code{Y - fitted} in each column.} 
#'  \item{"theta_hat"}{ A numerical matrix of size \code{p}-times-\code{1} of 
#'  estimated regression coefficients from \link{IRLS}.}
#'  \item{"beta_hat"}{ The final location estimate \code{beta} 
#'  evaluated at the \code{p} points from \code{tobs}, where \code{Y} was
#'  observed. A numerical vector of length \code{p}.}
#'  \item{"hat_values"}{ Diagonal terms of the (possibly penalized and weighted)
#'  hat matrix of the form \code{Z*solve(t(Z)*W*Z+n*lambda*H)*t(Z)*W}, 
#'  where \code{W} is the diagonal weight matrix in the final iteration 
#'  of \link{IRLS}. A numerical vector of length \code{m}.}
#'  \item{"weights"}{ The vector of weights given to the observations in the 
#'  final iteration of \link{IRLS}. A numerical vector of length \code{m}.}
#'  \item{"converged"}{ Indicator whether the \link{IRLS} procedure succefully 
#'  converged. Takes value 1 if IRLS converged, 0 otherwise.}
#'  \item{"obj_fun"}{ Objective function evaluated in the estimated theta hat.}
#' }
#' In case when \code{jcv="all"}, all these values are given for each 
#' cross-validation method considered. For \code{lambda} and \code{converged} 
#' provides a list of length 6 or 7 (depending on whether \code{custfun} is 
#' specified); for \code{fitted}, \code{beta_hat},
#' \code{hat_values}, and \code{weights} it gives a matrix with 6 or 7 
#' columns, each corresponding to one cross-validation method. 
#'
#' @references
#' Ioannis Kalogridis and Stanislav Nagy. (2025). Robust multidimensional 
#' location estimation from discretely sampled functional data. 
#' \emph{Under review}.
#'
#' @examples
#' d = 1                              # dimension of domain
#' m = 50                             # number of observation points
#' tobs = matrix(runif(m*d), ncol=d)  # location of observation points
#' n = 20                             # sample size
#' truemeanf = function(x)            # true location function 
#'   cos(4*pi*x[1])
#' truemean = apply(tobs,1,truemeanf) # discretized values of the true location
#' Y = replicate(n, truemean + rnorm(m)) # a matrix of functional data, size m*n
#' 
#' # introduce NAs
#' obsprob = 0.2                      # probability of a point being observed
#' B = matrix(rbinom(n*m,1,obsprob),ncol=n)
#' for(i in 1:m) for(j in 1:n) if(B[i,j]==0) Y[i,j] = NA
#' 
#' res = ts_location(Y, tobs=tobs, r=2, type="square", plotCV=TRUE)
#' 
#' @importFrom graphics abline par points plot title legend
#' @importFrom stats median
#' @export

ts_location = function(Y, tobs, r, type, alpha=1/2, tuning = NULL,
                       jcv = "all", vrs="C", method="IRLS",
                       plotCV=FALSE, lambda_grid=NULL,
                       lambda_length = 51, custfun=NULL,
                       resids.in = rep(1,length(Y)),
                       toler=1e-7, imax=1000,
                       tolerGCV=toler, imaxGCV=imax,
                       echo = FALSE){
  
  method = match.arg(method,c("IRLS", "ridge", "HuberQp", "QuantileQp"))
  type = match.arg(type,c("square","absolute","quantile","Huber","logistic"))
  if(method=="ridge" & type!="square") 
    stop("method 'ridge' available only for type 'square'.")
  if(method=="HuberQp" & type!="Huber") 
    stop("method 'HuberQp' available only for type 'Huber'.")
  if(method=="QuantileQp" & (type!="quantile" & type!="absolute" )) 
    stop("method 'QuantileQp' available only for type 'quantile' or 'absolute'.")
  
  jcv = match.arg(jcv,c("all", "AIC", "GCV", "GCV(tr)", "BIC", "rGCV", 
                        "rGCV(tr)", "custom"))
  if(jcv=="all") jcv = 0
  if(jcv=="AIC") jcv = 1
  if(jcv=="GCV") jcv = 2
  if(jcv=="GCV(tr)") jcv = 3
  if(jcv=="BIC") jcv = 4
  if(jcv=="rGCV") jcv = 5
  if(jcv=="rGCV(tr)") jcv = 6
  if(jcv=="custom") jcv = 7
  if(jcv==7 & is.null(custfun)) stop("With custom cross-validation, 
                                    cusfun must be provided.")
  
  # pre-processing for thin-plate splines
  tspr = ts_preprocess_location(Y, tobs, r)
  # attach(tspr)
  Y = tspr$Y; 
  Z = tspr$Z; H = tspr$H; Q = tspr$Q; Omega = tspr$Omega; Phi = tspr$Phi;
  w = tspr$w;
  degs = tspr$degs; M = tspr$M; r = tspr$r; 
  p = tspr$p; d = tspr$d; n = tspr$n;
  
  if(is.null(lambda_grid)){
    # define grid for search for lambda
    rho1 = -12  # search range minimium exp(rho1)
    rho2 = -1   # search range maximum exp(rho2)
    if(is.null(lambda_length)) lambda_length = 20
    lambda_grid = exp(seq(rho1,rho2,length=lambda_length-1))
  } else {
    if(!is.numeric(lambda_grid)) 
      stop("Grid for lambda values must contain numeric values.")
    if(any(lambda_grid<0)) 
      stop("Grid for lambda values must contain non-negative 
           values.")
    lambda_length = length(lambda_grid)
  }
  
  GCVfull <- Vectorize(
    function(x) GCV_location(x,
                             Z = Z, Y = Y, H = H, type=type, tuning = tuning, alpha=alpha, w=w, vrs=vrs,
                             method = method,
                             custfun = custfun,
                             resids.in = resids.in,
                             toler=tolerGCV, imax=imaxGCV))(lambda_grid)
  ncv = nrow(GCVfull)-2
  GCVconverged = GCVfull[ncv+1,]
  GCVic = GCVfull[ncv+2,]
  
  if(echo) print(paste(c("Numbers of iterations in IRLS", 
                         GCVic), collapse=", "))
  
  GCVfull = GCVfull[1:ncv,]
  cvnames = c("AIC","GCV","GCV(tr)","BIC","rGCV","rGCV(tr)",
              "custom")
  if(jcv==0) rownames(GCVfull) = cvnames[1:ncv] # if all the criteria are used 
  
  if(plotCV){
    if(jcv == 0){
      par(mfrow=c(3,2))
      for(i in 1:ncv){
        plot(log(GCVfull[i,])~log(lambda_grid),type="l",
             lwd=2,
             xlab=expression(log(lambda)),
             ylab="CV criterion",
             main = rownames(GCVfull)[i])
        points(log(GCVfull[i,])~log(lambda_grid),
               cex = 1-GCVconverged, col="red", pch=16)
        abline(v=log(lambda_grid[which.min(GCVfull[i,])]),lty=2)
      }
      par(mfrow=c(1,1))  
    } else {
      plot(log(GCVfull[jcv,])~log(lambda_grid),type="l",
           lwd=2, xlab=expression(log(lambda)),
           ylab="CV criterion",
           main = cvnames[jcv])
      points(log(GCVfull[jcv,])~log(lambda_grid),
             cex = 1-GCVconverged, col="red", pch=16)
      abline(v=log(lambda_grid[which.min(GCVfull[jcv,])]),lty=2)
    }
  }
  
  lopt = lambda_grid[apply(GCVfull,1,which.min)]
  
  if(jcv>0){
    lambda = lopt[jcv] # lambda parameter selected
    #
    if(method=="IRLS"){ 
      res = IRLS(Z,Y,lambda,H,type=type,alpha=alpha,w=w,tuning=tuning,vrs=vrs,sc=1, 
                 resids.in = resids.in, 
                 toler=toler, imax=imax)
      obj_fun_eval <- evaluate_objective(res$theta_hat, Y, Z, lambda*H,type=type,alpha=alpha,tuning=tuning)
    }
    if(method=="ridge"){
      res = ridge(Z,Y,lambda,H,w=w,vrs=vrs)
      obj_fun_eval <- evaluate_objective(res$theta_hat, Y, Z, lambda*H,type=type,alpha=alpha,tuning=1.345)
    }
    if(method=="HuberQp"){
      res = HuberQp(Z,Y,lambda,H,w=w,vrs=vrs,tuning=tuning)
      obj_fun_eval <- evaluate_objective(res$theta_hat, Y, Z, lambda*H,type=type,alpha=alpha,tuning=tuning)
    }
    if(method=="QuantileQp"){
      res = QuantileQp(Z,Y,lambda,H,alpha = alpha,w=w,vrs=vrs)
      obj_fun_eval <- evaluate_objective(res$theta_hat, Y, Z, lambda*H,type=type,alpha=alpha,tuning=1.345)
    }
    res_ts = transform_theta_location(res$theta_hat,tspr)
    if(method=="IRLS") return(
      list(lambda = lambda,
           fitted = res$fitted, 
           residuals = Y-res$fitted,
           theta_hat = res$theta_hat,
           beta_hat = res_ts$beta_hat,
           gamma_hat = res_ts$gamma_hat,
           delta_hat = res_ts$delta_hat,
           hat_values = res$hat_values,
           weights = res$weights, 
           converged = res$converged,
           obj_fun = obj_fun_eval))
    if(method=="ridge") return(
      list(lambda = lambda,
           fitted = res$fitted, 
           residuals = Y-res$fitted,
           theta_hat = res$theta_hat,
           beta_hat = res_ts$beta_hat,
           gamma_hat = res_ts$gamma_hat,
           delta_hat = res_ts$delta_hat,
           hat_values = res$hat_values,
           obj_fun = obj_fun_eval))
    if(method=="HuberQp") return(
      list(lambda = lambda,
           fitted = res$fitted, 
           residuals = Y-res$fitted,
           theta_hat = res$theta_hat,
           beta_hat = res_ts$beta_hat,
           gamma_hat = res_ts$gamma_hat,
           delta_hat = res_ts$delta_hat,
           hat_values = res$hat_values,
           obj_fun = obj_fun_eval))
    if(method=="QuantileQp") return(
      list(lambda = lambda,
           fitted = res$fitted, 
           residuals = Y-res$fitted,
           theta_hat = res$theta_hat,
           beta_hat = res_ts$beta_hat,
           gamma_hat = res_ts$gamma_hat,
           delta_hat = res_ts$delta_hat,
           hat_values = res$hat_values,
           obj_fun = obj_fun_eval))
  } else {
    n = length(Y)
    fitted = matrix(nrow=n,ncol=ncv)
    betahat = matrix(nrow=nrow(tobs), ncol=ncv)
    gammahat = matrix(nrow=nrow(tobs),ncol=ncv)
    deltahat = matrix(nrow=choose(d+r-1,d),ncol=ncv)
    thetahat = matrix(nrow=nrow(tobs), ncol=ncv)
    hatvals = matrix(nrow=n, ncol=ncv)
    weights = matrix(nrow=n, ncol=ncv)
    converged = rep(NA,ncv)
    obj_fun = rep(NA,ncv)
    for(jcv in 1:ncv){
      lambda = lopt[jcv] # lambda parameter selected
      #
      if(method=="IRLS"){
        res = IRLS(Z,Y,lambda,H,type=type,tuning=tuning,alpha=alpha,w=w,vrs=vrs,sc=1, 
                   resids.in = resids.in, 
                   toler=toler, imax=imax)
        obj_fun_eval <- evaluate_objective(res$theta_hat, Y, Z, lambda*H,type=type,alpha=alpha,tuning=tuning)
      }
      if(method=="ridge"){
        res = ridge(Z,Y,lambda,H,w=w,vrs=vrs)
        obj_fun_eval <- evaluate_objective(res$theta_hat, Y, Z, lambda*H,type=type,alpha=alpha,tuning=1.345)
      }
      if(method=="HuberQp"){
        res = HuberQp(Z,Y,lambda,H,w=w,vrs=vrs,tuning=tuning)
        obj_fun_eval <- evaluate_objective(res$theta_hat, Y, Z, lambda*H,type=type,alpha=alpha,tuning=tuning)
      }
      if(method=="QuantileQp"){
        res = QuantileQp(Z,Y,lambda,H,alpha=alpha,w=w,vrs=vrs)
        obj_fun_eval <- evaluate_objective(res$theta_hat, Y, Z, lambda*H,type=type,alpha=alpha,tuning=1.345)
      }
      res_ts = transform_theta_location(res$theta_hat,tspr)
      fitted[,jcv] = res$fitted
      thetahat[,jcv] = res$theta_hat
      betahat[,jcv] = res_ts$beta_hat
      gammahat[,jcv] = res_ts$gamma_hat
      deltahat[,jcv] = res_ts$delta_hat
      hatvals[,jcv] = res$hat_values
      if(method=="IRLS"){
        weights[,jcv] = res$weights
        converged[jcv] = res$converged
      }
      obj_fun[jcv] = obj_fun_eval
    }
    return(list(lambda = lopt,
                fitted = fitted,
                residuals = Y-fitted,
                theta_hat = thetahat,
                beta_hat = betahat,
                gamma_hat = gammahat,
                delta_hat = deltahat,
                hat_values = hatvals,
                weights = weights, 
                converged = converged,
                obj_fun = obj_fun))
  }
}
