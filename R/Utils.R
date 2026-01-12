#' Generate a dataset for scalar-on-function linear regression
#'
#' Generates the data according to the model 
#' \code{Y = alpha0 + int beta0(t)*X(t) dt + epsilon}, where \code{X} are the 
#' regressor functions, \code{alpha0} is a real parameter, \code{beta0} is a 
#' functional parameter, and \code{epsilon} are independent error terms. 
#' The domain of \code{X} and \code{beta0} can be one-dimensional, or also
#' higher-dimensional. The functional data \code{X} are assumed to be observed
#' only in an incomplete grid of \code{p} points in its domain.
#'
#' @param alpha0 A real value of the intercept parameter.
#'
#' @param beta0 A function parameter. A real function whose arguments are
#' \code{d}-dimensional vectors.
#'
#' @param n Sample size.
#'
#' @param d Dimension of the domain of \code{X} and \code{beta0}.
#' 
#' @param p The size of the discretization grid where the functions \code{X} are
#' observed. The observation points are sampled randomly from the full grid of
#' all \code{p1^d} points in the domain.
#' 
#' @param bfX A list of basis functions for \code{X}. Each element of \code{bfX}
#' contains the list of values of a basis function when evaluated in a complete
#' grid of observations, that is \code{p1^d} real values.
#' 
#' @param bcX A matrix of basis coefficients that are combined with the basis
#' functions from \code{bfX} to get the regressor functions \code{X}. The matrix
#' must be of dimensions \code{n}-times-\code{K}, where \code{n} is the sample
#' size, and \code{K} is the length of the list \code{bfX}.
#' 
#' @param sd.noiseEps Standard deviation for the noise component \code{eps}. 
#' Noise is generated to be centered Gaussian, independent of \code{X}.
#' 
#' @param obs_only Indictar of whether to generate the response \code{Y}
#' the functions \code{X} and \code{beta0} should be used in the complete 
#' grid of size \code{p1^d} (\code{obs_only=FALSE}), or only in the grid 
#' of size \code{p} where \code{X} is observed (\code{obs_only=FALSE}). By
#' default set to \code{FALSE}.
#' 
#' @param p1 The complete size of the discretization grid where the 
#' (practically unobservable) complete version of the function \code{X} is 
#' known. \code{p1} corresponds to the complete grid in dimension 1; if 
#' \code{d>1}, the same grid is replicated also in the other dimensions and the
#' resulting size of the domain is \code{p1^d}. By default, 
#' \code{p1=101}.
#' 
#' @param sd.noiseX Standard deviation of an additional noise component that is
#' added to the regressors \code{X}. The final regressor functions are obtained
#' by taking a linear combination of basis functions from \code{bfX} weighted
#' by the coefficients from \code{bcX}, and then adding indepenent centered 
#' normal noise to each of the \code{p1^d} of its values with standard
#' deviation \code{sd.noiseX}. By default, \code{sd.noiseX = 0}.
#' 
#' @return A list composed of several elements. Some of these elements are
#' assumed to be observed; some are given for post-estimation diagnostic 
#' purposes. The elements that are assumed to be observed are:
#' \itemize{
#'  \item{"tgrid"}{ A complete grid of observation points in one dimension, 
#'  replicated in each dimension if \code{d>1}. An equidistant grid of 
#'  \code{p1} points in the interval [0,1].}
#'  \item{"X"}{ Matrix of observed values of \code{X} of size 
#'  \code{n}-times-\code{p}, one row per observation, columns corresponding to the 
#'  positions in the rows of \code{tobs}.}
#'  \item{"Y"}{ A vector of observed responses, size \code{n}.}
#'  \item{"tobs"}{ Domain locations for the observed points of \code{X}. Matrix
#'  of size \code{p}-times-\code{d}, one row per domain point.}
#' }
#' Additional data for diagnostic purposes that is not assumed to be observed:
#' \itemize{
#'  \item{"Xfull"}{ A \code{d+1}-dimensional array of completely observed 
#'  values of functions \code{X} of size 
#'  \code{n}-times-\code{p1}-...-times-\code{p1}. The first dimension is for 
#'  observations, additional ones for dimensions of the
#'  domain, i.e. \code{Xfull[2,4,5]} stands for the function value of function 2
#'  at the domain points \code{(tgrid[4],tgrid[5])} for \code{d=2}.}
#'  \item{"betafull"}{ A \code{d}-dimensional array of the values of the 
#'  regression function \code{beta0} of size \code{p1}-times-...-times-\code{p1}.
#'  Interpretation of values as for the functions of \code{Xfull}.}
#'  \item{"tfull"}{ A flattened matrix of the full grid of all observation 
#'  points. A matrix of size \code{p1^d}-times-\code{d}, each row 
#'  corersponding to the domain vector of a single point where \code{Xfull} and
#'  \code{betafull} are evaluated.}
#'  \item{"betafull0"}{ A vector of all evaluated values of the function 
#'  \code{beta0} in the complete grid of \code{p1^d} observation points. A 
#'  numerical vector of length \code{p1^d}, each value corresponding to a 
#'  row of \code{tfull}. Idential to a flattened version of \code{betafull}.}
#'  \item{"truemean"}{ A vector of length \code{n} of the true conditional mean 
#'  values of Y given the i-th regressor from \code{X}, without the intercept
#'  \code{alpha0}.}
#'  \item{"eps"}{ Vector of erros used in the regression. Naturally, 
#'  \code{Y = truemean + eps + alpha0}.}
#'  \item{"d"}{ Dimension of the domain of regressors.}
#' }
#'
#' @examples
#' d = 1       # dimension of domain
#' n = 50      # sample size
#' p = 10      # number of observed points in domain
#' p1 = 101 # size of the complete grid
#' alpha0 = 3  # intercept
#' beta0 = function(t) t^2 # regression coefficient function
#' sd.noiseEps = 0.1 # noise standard deviation
#' tgrid = seq(0,1,length=p1) # full grid of observation points
#' 
#' # Generating the basis functions for X
#' basis_fun = function(t,k) return(sin(2*pi*k*t)+t^k)
#' 
#' K = 5 # number of basis functions used in the expansion
#' bfX = list()
#' for(k in 1:K) bfX[[k]] = basis_fun(tgrid,k)
#' bcX = matrix(rnorm(n*K,mean=3,sd=5),ncol=K)  # basis coefficients
#' 
#' gen = generate(alpha0, beta0, n, d, p, bfX, bcX, sd.noiseEps)
#' 
#' # plot the first regressor function
#' plot(gen$Xfull[1,]~tgrid,type="l",main="Predictor X",lwd=2)
#' 
#' # plot the true mean values against the responses
#' plot(gen$truemean+alpha0, gen$Y)
#' abline(a=0,b=1)
#' 
#' # Check that the response is modelled as alpha0 + truemean + eps
#' all.equal(gen$Y,alpha0 + gen$truemean + gen$eps)
#' 
#' # Check that the truemean is an approximation to the integral of X and beta0
#' all.equal(apply(gen$Xfull,1,function(x) mean(x*gen$betafull)),gen$truemean)
#' 
#' # Check that betafull0 is flattened version of betafull
#' all.equal(c(gen$betafull), gen$betafull0)
#' 
#' @export

generate = function(alpha0, beta0, n, d, p,  
                    bfX, bcX, sd.noiseEps, obs_only = FALSE, p1 = 101, sd.noiseX=0){
  
  pfull = p1^d  # complete grid of observations
  tgrid = seq(0,1,length=p1)
  if(pfull<p1) stop("Less points in full discretization that p")
  #
  K = length(bfX) # number of basis functions in the expansion of X
  
  # generate full covariates X
  if(d==1){
    ifull = 1:p1
    tfull = matrix(tgrid,ncol=1)
    Xfull = array(rnorm(pfull*n,sd=sd.noiseX),
                  dim=c(n,p1))
    for(i in 1:n) for(k in 1:K) Xfull[i,] = 
      Xfull[i,] + bfX[[k]]*bcX[i,k]
    # predictors are iid normally distributed noise plus signal
    # plot(Xfull[1,]~tgrid,type="l",main="Predictor X",lwd=2)
  }
  if(d==2){
    ifull = expand.grid(1:p1,1:p1)
    tfull = expand.grid(tgrid,tgrid)
    Xfull = array(rnorm(pfull*n,sd=sd.noiseX),
                  dim=c(n,p1,p1)) # predictors
    for(i in 1:n) for(k in 1:K) Xfull[i,,] = 
      Xfull[i,,] + bfX[[k]]*bcX[i,k]
    # image(Xfull[1,,],main="Predictor X")
  }
  if(d>2){
    if(K>0) warning("For d>2 the basis in X is ignored, X is just noise.")
    ifull = expand.grid(replicate(d,1:p1,simplify=FALSE))
    tfull = expand.grid(replicate(d,tgrid,simplify=FALSE))
    Xfull = array(rnorm(pfull*n,sd=sd.noiseX),
                  dim=c(n,rep(p1,d))) # predictors
  }
  
  if(d==1){
    betafull = Vectorize(beta0)(tgrid)
    betafull0 = betafull
  }
  if(d==2){
    betafull0 = apply(tfull,1,beta0)
    betafull = matrix(betafull0,nrow=p1)
  }
  if(d>2){
    betafull0 = apply(tfull,1,beta0)
    betafull = array(betafull0,dim=rep(p1,d))
  }
  
  ###
  # Generating discretely observed data
  ###
  tind = sort(sample(1:pfull,p))
  tobs = tfull[tind,,drop=FALSE] # observed points
  
  # Observed values of X, only at the p points in tobs
  if(d==1){
    X = Xfull[,tind]
    # beta = betafull[tind]
  }
  if(d==2){
    Xff = matrix(Xfull,nrow=n)  # flattened X matrix
    X = Xff[,tind]              # observed points of X
    # beta = betafull0[tind]      # beta0 at observed points
  }
  if(d>2){
    Xff = matrix(Xfull,nrow=n)
    X = Xff[,tind]
    # beta = betafull0[tind]
  }
  
  # plot(apply(Xfull,1,function(x) mean(x*betafull))~
  #        apply(X,1,function(x) mean(x*beta)),
  #      xlab="Approximated integrals from discrete functions",
  #      ylab="Exact integrals")
  # abline(a=0,b=1)
  
  ###
  # model
  ###
  eps = rnorm(n,sd=sd.noiseEps)          # random errors
  
  if(obs_only){
    truemean = apply(X,1,function(x) mean(x*betafull0))
  } else {
    truemean = apply(Xfull,1,function(x) mean(x*betafull))
  }
  Y = alpha0 + truemean + eps
  
  return(list(tgrid=tgrid, X=X,Y=Y,tobs=tobs,
              Xfull=Xfull,betafull=betafull,betafull0=betafull0,
              tfull=tfull,
              truemean=truemean,eps=eps,d=d))
}

#' Reconstruct and visualise the slope function beta0 and its estimate
#'
#' Based on the thin-splines estimates of \code{beta0} obtained from 
#' functions \link{IRLS} or \link{ridge}, this function reconstructs the 
#' complete estimator of the regression function \code{beta0} on its full 
#' domain. It also visualises the results. 
#'
#' @param ts_prep Output of the function \link{ts_preprocess}.
#' 
#' @param theta Estimated coefficients, output of \link{IRLS} or \link{ridge}
#' after passed through \link{transform_theta}. Either \code{theta} or
#' \code{lambda} must be provided. If \code{theta} is given, \code{lambda}
#' is disregarded.
#'
#' @param lambda Tuning parameter, a non-negative real number. Either \code{theta} or
#' \code{lambda} must be provided. If \code{theta} is not given, \link{IRLS}
#' is fitted with parameter \code{lambda}.
#'
#' @param Y Vector of responses of length \code{n}.
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
#' @param p1 Number of equidistant discretization points in the interval [0,1]
#' in each dimension \code{d} of the domain. The total number of discretization
#' points at which the estimator of \code{beta0} will be generated is 
#' \code{p1^d}.
#' 
#' @param betafull The true function \code{beta0} to be plotted along its
#' estimator. A \code{d}-dimensional array of the values of the 
#' regression function \code{beta0} of size \code{p1}-times-...-times-\code{p1}.
#' Interpretation of values as for the functions of \code{Xfull} in function
#' \link{generate}. Can be provided directly as the output of \link{generate}.
#' 
#' @param main Title of the displayed plot.
#' 
#' @return A plot with the true and the estimated function \code{beta0}, and
#' the corresponding output of the function \link{IRLS}.
#'
#' @examples
#' d = 1       # dimension of domain
#' n = 50      # sample size
#' p = 10      # number of observed points in domain
#' p1 = 101 # size of the complete grid
#' alpha0 = 3  # intercept
#' beta0 = function(t) t^2 # regression coefficient function
#' sd.noiseEps = 0.1 # noise standard deviation
#' tgrid = seq(0,1,length=p1) # full grid of observation points
#' basis_fun = function(t,k) return(sin(2*pi*k*t)+t^k)
#' K = 5 # number of basis functions used in the expansion
#' bfX = list()
#' for(k in 1:K) bfX[[k]] = basis_fun(tgrid,k)
#' bcX = matrix(rnorm(n*K,mean=3,sd=5),ncol=K)  # basis coefficients
#' 
#' # Generate the raw data
#' gen = generate(alpha0, beta0, n, d, p, bfX, bcX, sd.noiseEps)
#' 
#' # Preprocess the data
#' X = gen$X; Y = gen$Y; tobs = gen$tobs
#' m = 2   # order of thin-plate splines
#' tspr = ts_preprocess(X,tobs,m)
#' 
#' # Reconstruction by specifying regression method
#' reconstruct(tspr, theta = NULL,lambda = 1e-3,Y,"absolute",
#'   p1=p1, betafull=gen$betafull, main="L1 regression")
#'   
#' # Reconstruction by using theta from previously fitted model
#' res = IRLS(tspr$Z,Y,1e-3,tspr$H,"absolute")
#' reconstruct(ts_prep = tspr, theta = transform_theta(res$theta_hat,tspr), 
#' lambda=NULL, Y = Y, type = "absolute", 
#' p1 = p1, betafull = gen$betafull)
#' 
#' @export

reconstruct = function(ts_prep = NULL, 
                       theta=NULL, lambda=NULL,
                       Y=NULL, type=NULL, alpha=1/2,
                       p1 = 101,
                       betafull=NULL,main=NULL){
  tgrid = seq(0,1,length=p1)
  # plots the complete estimated function beta
  if(is.null(theta) & is.null(lambda)) 
    stop("Either theta or lambda must be given.")
  if(is.null(ts_prep)) stop("ts_prep must be given.")
  if(is.null(theta)){
    if(any(is.null(ts_prep$Z),is.null(Y),is.null(ts_prep$H),is.null(type)))
      stop("Z, Y, H, and type must be provided if theta is to be estimated.")
    # estimate directly
    res = IRLS(ts_prep$Z,Y,lambda,ts_prep$H,type, alpha)
    theta = transform_theta(res$theta_hat,ts_prep)
  } else res = NULL
  
  d = ts_prep$d
  if(d>2) stop("Visualisation only for d<=2.")
  #
  tfull = expand.grid(replicate(d,tgrid,simplify=FALSE))
  #
  tobs=ts_prep$tobs
  m = ts_prep$m
  degs = ts_prep$degs
  
  # Fast matrix of Euclidean distances
  if(d==1) Emfull = abs(outer(tfull[,1],tobs[,1],"-"))
  if(d==2) Emfull = sqrt(outer(tfull[,1],tobs[,1],"-")^2 + 
                           outer(tfull[,2],tobs[,2],"-")^2)
  
  # Matrix Omega for the penalty term
  Omegafull = eta(Emfull,d,m)
  # Matrix Phi, monomials evaluated at t
  Phifull = apply(degs,1,function(x) apply(t(tfull)^x,2,prod))
  #
  b_hatfull0 = c(Omegafull%*%theta$gamma_hat + Phifull%*%theta$delta_hat)
  b_hatfull = matrix(b_hatfull0,nrow=p1)
  
  if(d==1){
    if(!is.null(betafull)){
      plot(c(b_hatfull,betafull)~c(tgrid,tgrid),type="n",
           xlab="t",ylab=expression(beta[0]),main=main)
      lines(betafull~tgrid,lwd=2)
      lines(b_hatfull~tgrid,col=2,lwd=2)
      legend("topleft",c("True","Estimated"),col=1:2,lwd=2)
    } else {
      plot(c(b_hatfull)~c(tgrid),type="n",
           xlab="t",ylab=expression(beta[0]),main=main)
      lines(b_hatfull~tgrid,col=2,lwd=2)
      legend("topleft",c("Estimated"),col=2,lwd=2)
    }
    rug(tobs)
  }
  if(d==2){
    if(!is.null(betafull)) op = par(mfrow=c(1,2))
    image(b_hatfull,main=main)
    if(!is.null(betafull)) image(betafull,main="True beta")
    if(!is.null(betafull)) par(op)
  }
  return(res)
}

#' Interpolating functional data using thin-plate splines
#'
#' Given a discretely (and possibly) irregularly observed sample of functional
#' data, this function performs thin-plate spline interpolation of each 
#' functional datum. As an output, a matrix or function values of densely
#' observed interpolated functional data in a common grid of the domain 
#' is obtained.
#' 
# p.out = number of function values in the interpolated function in each d
# I = c(a,b) for a complete domain of the form of the square [a,b]^d
#'
#'
#' @param Xtobs A list of length \code{n}, each element of \code{Xtobs} is a 
#' list with two elements corresponding to the \code{i}-th function: 
#' (a) \code{X} A numerical vector of length \code{p} of observed values of 
#' the \code{i}-th function, where \code{p} is the number of observations for
#' this function. (b) \code{tobs} Domain locations for the observed points of 
#' the \code{i}-th function corresponding to \code{X}. Matrix of size 
#' \code{p}-times-\code{d}, one row per domain point. Here, \code{d} is the 
#' dimension of the domain of \code{X}.
#' 
#' @param r Order of the thin-plate spline, positive integer.
#' 
#' @param p.out Number of equi-distant points in which the densely observed
#' thin-plate spline interpolation of each function should be evaluated. In
#' dimension \code{d}, a total of \code{p.out^d} points are evaluated in the 
#' hyper-cube \code{[I[1],I[2]]^d}, and at these points, 
#' each function is evaluated.
#' 
#' @param I Two numbers \code{I[1]<I[2]} determining the box 
#' \code{[I[1],I[2]]^d} where all the interpolated functions are evaluated. 
#' 
#' @param solve.method Indicator of which solver of systems of linear equations
#' to use for the inversion of the matrices in the interpolation procedure. 
#' Possible options are \code{"C"} for the solver from \code{Armadillo} library
#' in \code{C++}, or \code{R} for the function \link[base]{solve} from \code{R}. 
#' 
#' @return A list of values:
#' \itemize{
#'  \item{"X"}{ A matrix of size \code{n}-times-\code{p.out^d} of the function
#'  values of the densely evaluated interpolated data. One function per row.}
#'  \item{"tobs"}{ A matrix of size \code{p.out^d}-times-\code{d} with all the
#'  coodinates of the points where the function values are evaluated.
#'  }
#'  }
#'
#' @examples
#' 
#' # Construct the irregularly observed functional data
#' Xtobs = list()
#' 
#' # Number of observations and dimension of the domain
#' n = 10
#' d = 1
#'
#'# Generating functions as random noise
#'for(i in 1:n){
#'  p = 10 + rpois(1, 5) # Number of observed points for the i-th function
#'  tobs = matrix(runif(p*d),ncol=d) # Points in the domain
#'  X = rnorm(p,sd=0.05) # Observed function values
#'  Xtobs[[i]] = list(X=X, tobs=tobs)
#'}
#'  
#'# Thin-plate spline interpolation of order r=4
#'intr = ts_interpolate(Xtobs, r = 4)
#'  
#'# Visalization of the i-th function
#'i = 1
#'if(d==1){
#'  plot(Xtobs[[i]]$X~Xtobs[[i]]$tobs[,1])  
#'  lines(intr$X[i,]~intr$tobs[,1])
#'}
#'if(d==2){
#'  rgl::plot3d(Xtobs[[i]]$tobs[,1], Xtobs[[i]]$tobs[,2], Xtobs[[i]]$X)
#'  rgl::points3d(intr$tobs[,1], intr$tobs[,2], intr$X[i,],col="red",add=TRUE)
#'}
#'
#'@export

ts_interpolate = function(Xtobs, r, p.out = 101, I=c(0,1), 
                          solve.method=c("C","R")){
  
  solve.method = match.arg(solve.method)
  n = length(Xtobs)
  d = ncol(Xtobs[[1]]$tobs)
  
  tobs.out1 = seq(I[1],I[2],length=p.out)
  tobs.out = as.matrix(expand.grid(replicate(d, tobs.out1, simplify=FALSE)))
  p.out = nrow(tobs.out)  
  res = matrix(nrow=n, ncol=p.out)
  
  for(i in 1:n){
    X = Xtobs[[i]]$X
    tobs = Xtobs[[i]]$tobs
    p = length(X)
    if(d!=ncol(tobs))
      stop("Not all functions have the same dimension of the domain")
    
    if(p!=nrow(tobs)) 
      stop("Number of columns of X must equal number of rows of t")
    
    M = choose(r+d-1,d)
    if(p-M<=0) stop(paste("p must be larger than",M))
    if(2*r<=d) stop(paste("r must be larger than",ceiling(d/2)))
    
    # Monomials phi of order <r
    allcom = expand.grid(replicate(d, 0:(r-1), simplify=FALSE))
    degs = as.matrix(allcom[rowSums(allcom)<r,,drop=FALSE])
    M = nrow(degs)
    if(M!=choose(r+d-1,d)) stop("Error in degrees of polynomials")
    
    # Fast matrix of Euclidean distances between tobs
    Em = matrix(sqrt(rowSums(
      apply(tobs,2,function(x) outer(x,x,"-")^2))),nrow=nrow(tobs))
    
    # Matrix Omega for the penalty term
    Omega = eta(Em,d,r)
    
    # Matrix Phi, monomials evaluated at tobs
    Phi = apply(degs,1,function(x) apply(t(tobs)^x,2,prod))
    
    A = rbind(cbind(Omega, Phi),cbind(t(Phi), matrix(0, nrow=M, ncol=M)))
    b = c(X,rep(0,M))
    if(solve.method=="R") coefs = solve(A,b, tol = 1e-21)
    if(solve.method=="C") coefs = c(solveC(A,matrix(b,ncol=1)))
    #
    gamma = coefs[1:p]
    delta = coefs[-(1:p)]
    
    # Matrices Omega and Phi for evaluating interpolated function
    Em.out = distnAB(tobs.out,tobs,p.out,p,d)
    Omega.out = eta(Em.out,d,r)
    Phi.out = apply(degs,1,function(x) apply(t(tobs.out)^x,2,prod))
    #
    res[i,] = Omega.out%*%gamma + Phi.out%*%delta
  }
  return(list(X=res, tobs = tobs.out))
}

#' Quantile regression based on FPCA of Kato (2012)
#'
#' Fits a (robust) quantile regression estimator in a scalar-on-function 
#' regression problem with discretely observed predictors. The estimator is 
#' based on quantile regression applied in the context of functional principal
#' components. The tuning parameter \code{J} (the number of principal 
#' compontents chosen) is selected using cross-validation.
#'
#' @param X Matrix of observed values of \code{X} of size 
#'  \code{n}-times-{p}, one row per observation, columns corresponding to the 
#'  positions in the rows of \code{t}.
#'
#' @param Y Vector of responses of length \code{n}.
#'  
#' @param t Domain locations for the observed points of \code{X}. A numerical
#'  vector of length \code{p}.
#' 
#' @param tau Order of the quantile to be used in the estimator. By default, 
#' \code{tau=0.5}, which corresponds to the regression median.
#' 
#' @param Kmax Maximum number of functional principal components to be 
#' considered. By default, set to \code{Kmax=10}.
#' 
#' @param n_fine Number of equi-spaced points in the domain for the complete
#' function. By default, set to \code{n_fine=200}.
#'
#' @return A list composed of
#' \itemize{
#'  \item{"fitted_Y"}{ A vector of \code{n} fitted values.}
#'  \item{"alpha_hat"}{ Estimate of \code{alpha0}, a numerical value.}
#'  \item{"beta_hat"}{ Estimate of the regression function \code{beta0} 
#'  evaluated at the \code{p} points from \code{tobs}, where \code{X} was
#'  observed.}
#'  \item{"gcv_vals"}{ Observed values of the generalized cross-validation 
#'  criteria for the selection of the dimension. A numerical vector of length
#'  \code{Kmax}.}
#'  \item{"best_m"}{ Selected dimension using generalized cross-validation.}
#' }
#'
#' @references
#' Kato, K. (2012). Estimation in functional linear quantile regression. 
#' \emph{The Annals of Statistics}, 40(6), 3108-3136.
#'
#' @examples
#' n = 50      # sample size
#' p = 10      # dimension of predictors
#' X = matrix(rnorm(n*p),ncol=p) # design matrix
#' Y = X[,1]   # response vector
#' tobs = sort(runif(p)) 
#' 
#' (res = qr_fpca_piecewise_fine(X, Y, tobs))
#' 
#' @importFrom quantreg rq
#' @export

qr_fpca_piecewise_fine <- function(X, Y, t, tau = 0.5, K_max = 10, n_fine = 200) {
  
  # Piecewise constant interpolating function 
  
  piecewise_constant <- function(t_obs, x_obs) {
    
    # t_obs is a vector containing the discretization points
    # x_obs is a vector containing the values to be interpolated
    
    function(t) {
      sapply(t, function(tt) {
        if(tt >= t_obs[length(t_obs)]) {
          return(x_obs[length(x_obs)])
        } else {
          j <- max(which(t_obs <= tt))
          return(x_obs[j])
        }
      })
    }
  }
  
  n <- nrow(X) # number of observations
  p <- ncol(X) # number of discretization points
  
  # X is the n x p predictor matrix
  # t is a vector containing the discretization points (the grid is assumed to be common)
  # tau is the quantile to be estimated (by default 0.5 corresponding to the median)
  # K_max is the number of candidate eigenfunctions to keep
  # n_fine is the number of points that the random functions will be evaluated on after interpolation
  
  # Create piecewise constant functions for each one of the X_i
  X_func <- lapply(1:n, function(i) piecewise_constant(t, X[i,]))
  t_fine <- seq(min(t), max(t), length.out = n_fine) # equispaced discretization
  
  dt <- diff(c(t_fine, t_fine[n_fine] + (t_fine[n_fine]-t_fine[n_fine-1]))) # adjust the rightmost point
  
  # Evaluate curves on fine grid
  X_fine <- t(sapply(X_func, function(f) f(t_fine)))  # n x n_fine
  
  # Mean function on fine grid
  mean_fun_fine <- colMeans(X_fine)
  
  # Covariance matrix (kernel) on fine grid
  C <- matrix(0, n_fine, n_fine)
  for(i in 1:n) {
    Xc <- X_fine[i,] - mean_fun_fine
    C <- C + (Xc %*% t(Xc)) * outer(dt, dt)
  }
  C <- C / n
  
  # Eigen-decomposition
  eig <- eigen(C, symmetric = TRUE)
  lambda <- eig$values # eigenvalues in descending order
  phi_fine_all <- eig$vectors  # all eigenfunctions
  phi_fine_all <- phi_fine_all / sqrt(dt) # renormalize for orthonormality in L2
  
  
  # Riemann sums for all FPCA scores up to K_max
  xi_all <- matrix(0, n, K_max)
  for(i in 1:n) {
    for(k in 1:K_max) {
      xi_all[i,k] <- sum((X_fine[i,] - mean_fun_fine) * phi_fine_all[,k] * dt)
    }
  }
  
  # GCV selection of number of eigenfunctions
  
  gcv_vals <- numeric(K_max)
  for(m in 1:K_max) {
    xiK <- xi_all[,1:m, drop=FALSE]
    fit <- rq(Y ~ xiK, tau = tau) # Quantile regression with the rq function of the quantreg package
    fitted_Y <- fitted(fit) # fitted values
    rho <- function(r, tau) r * (tau - as.numeric(r < 0)) # check loss function
    gcv_vals[m] <- sum(rho(Y - fitted_Y, tau)) / (n - (m + 1)) # GCV criterion, see Section 4 Kato (2012)
  }
  
  best_m <- which.min(gcv_vals) # minimizer of the GCV criterion
  
  # Final fit using best_m
  xiK <- xi_all[,1:best_m, drop=FALSE]
  fit <- rq(Y ~ xiK, tau = tau)
  beta_hat <- coef(fit)[-1] # exclude intercept
  alpha_hat <- coef(fit)[1] # intercept
  
  # Functional coefficient on fine grid (linear combination of the eigenfunctions)
  beta_fun_fine <- phi_fine_all[,1:best_m] %*% beta_hat
  
  # Convert to original t
  beta_fun <- sapply(t, function(tt) beta_fun_fine[max(which(t_fine <= tt))])
  mean_fun <- sapply(t, function(tt) mean_fun_fine[max(which(t_fine <= tt))])
  phi <- matrix(0, nrow=length(t), ncol=best_m)
  for(k in 1:best_m) {
    phi[,k] <- sapply(t, function(tt) phi_fine_all[max(which(t_fine <= tt)),k])
  }
  
  # Fitted scalar responses
  fitted_Y <- alpha_hat + xiK %*% beta_hat
  
  return(list(
    mean_fun = mean_fun,
    eigenfunctions = phi,
    eigenvalues = lambda[1:best_m],
    scores = xiK,
    qr_fit = fit,
    fitted_Y = fitted_Y,
    alpha_hat = alpha_hat,
    beta_fun = beta_fun,
    t = t,
    gcv_vals = gcv_vals,
    best_m = best_m
  ))
}