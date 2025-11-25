#' Cross-Validation for the IRLS procedure
#'
#' Provides the cross-validation indices from \link{GCV_crit} in
#' conjunction with the \link{IRLS} function directly as an argument of the
#' parameter \code{lambda}. 
#'
#' @param lambda A candidate parameter value; non-negative real number.
#'
#' @param Z Data matrix of dimension \code{n}-times-\code{p}, where \code{n} is
#' the number of observations, \code{p} is the dimension.
#'
#' @param Y Vector of responses of length \code{n}.
#'
#' @param H Penalty matrix of size \code{p}-times-\code{p} that
#' is used inside the quadratic term for penalizing estimated parameters.
#' 
#' @param type The type of the loss function used in the minimization problem.
#' Accepted are \code{type="absolute"} for the absolute loss \code{rho(t)=|t|}; 
#' \code{type="square"} for the square loss \code{rho(t)=t^2}; 
#' \code{type="Huber"} for the Huber loss \code{rho(t)=t^2/2} if 
#' \code{|t|<tuning} and \code{rho(t)=tuning*(|t|-tuning/2)} otherwise; and 
#' \code{type="logistic"} for the logistic loss 
#' \code{rho(t)=2*t + 4*log(1+exp(-t))-4*log(2)}.
#' 
#' @param sc Scale parameter to be used in the IRLS. By default \code{sc=1}, 
#' that is no scaling is performed.
#' 
#' @param vrs Version of the algorhitm to be used in function \link{IRLS}; 
#' either \code{vrs="C"} for the \code{C++} version, or \code{vrs="R"} for the 
#' \code{R} version. Both should give (nearly) identical results, see 
#' \link{IRLS}.
#' 
#' @param custfun A custom function combining the residuals \code{resids} and
#' the hat values \code{hats}. The result of the function must be numeric, 
#' see \link{GCV_crit}.
#' 
#' @param resids.in Initialization of the vector of residuals used to launch 
#' the IRLS algorithms. Optional.
#' 
#' @param toler A small positive constant specifying the tolerance level for 
#' terminating the algorithm. The prcedure stops if the maximum absolute 
#' distance between the residuals in the previous iteration and the new 
#' residuals drops below \code{toler}.
#' 
#' @param imax Maximum number of allowed iterations of IRLS. 
#'
#' @details Function \code{custfun} has two arguments 
#' corresponding to \code{resids} and \code{hats}. The output of the function 
#' must be numeric. 
#' 
#' @return A named numerical vector of values. The length of the vector depends
#' on the input. The vector contains the values:
#' \itemize{
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
#' }
#'
#' @examples
#' n = 50      # sample size
#' p = 10      # dimension of predictors
#' Z = matrix(rnorm(n*p),ncol=p) # design matrix
#' Y = Z[,1]   # response vector
#' lambda = 1  # tuning parameter for penalization
#' H = diag(p) # penalty matrix
#' type = "absolute" # absolute loss
#' 
#' # Run with the IRLS procedure
#' res = IRLS(Z, Y, lambda, H, type)
#' with(res,GCV_crit(resids,hat_values))
#' with(res,GCV_crit(resids,hat_values,custfun = function(r,h) sum(r^2)))
#'     
#' GCV(lambda,Z,Y,H,type,custfun = function(r,h) sum(r^2))

GCV <- function(lambda, Z, Y, H, type, sc = 1, 
                vrs="C", custfun=NULL, 
                resids.in = rep(1,length(Y)),
                toler=1e-7, imax=1000){
  # Generalized cross-validation
  ncv = 6
  if(!is.null(custfun)) ncv = 7
  vrs = match.arg(vrs, c("C", "R"))
  fit.r <- IRLS(Z, Y, lambda, H, type, sc = sc, vrs=vrs, 
                resids.in = resids.in,
                toler=toler, imax=imax)
  GCV.scores <- GCV_crit(fit.r$resids,fit.r$hat_values,
                         custfun=custfun)
  return(c(GCV.scores, fit.r$converged, fit.r$ic))
}

#' Cross-validation for location estimation
#'
#' Provides the cross-validation indices from \link{GCV_crit} in
#' conjunction with the \link{IRLS} function and \link{ridge} function directly 
#' as an argument of the parameter \code{lambda}. Applicable for location
#' estimation and function \link{ts_location}.
#'
#' @param lambda A candidate parameter value; non-negative real number.
#'
#' @param Z Data matrix of dimension \code{n}-times-\code{p}, where \code{n} is
#' the number of observations, \code{p} is the dimension.
#'
#' @param Y Vector of responses of length \code{n}.
#'
#' @param H Penalty matrix of size \code{p}-times-\code{p} that
#' is used inside the quadratic term for penalizing estimated parameters.
#' 
#' @param type The type of the loss function used in the minimization problem.
#' Accepted are \code{type="absolute"} for the absolute loss \code{rho(t)=|t|}; 
#' \code{type="square"} for the square loss \code{rho(t)=t^2}; 
#' \code{type="Huber"} for the Huber loss \code{rho(t)=t^2/2} if 
#' \code{|t|<tuning} and \code{rho(t)=tuning*(|t|-tuning/2)} otherwise; and 
#' \code{type="logistic"} for the logistic loss 
#' \code{rho(t)=2*t + 4*log(1+exp(-t))-4*log(2)}.
#' 
#' @param w Vector of length \code{n} of weights attached to the elements of 
#' \code{Y}. If \code{w=NULL} (default), a constant vector with values 
#' \code{1/n} is used. 
#'
#' @param vrs Version of the algorhitm to be used in function \link{IRLS}; 
#' either \code{vrs="C"} for the \code{C++} version, or \code{vrs="R"} for the 
#' \code{R} version. Both should give (nearly) identical results, see 
#' \link{IRLS} and \link{ridge}.
#' 
#' @param method A method for estimating the fit. Possible options are 
#' \code{"IRLS"} for the IRLS algorithm, or \code{"ridge"} for ridge regression.
#' Ridge is applicable only if \code{type="square"}; this method is much faster,
#' but provides only a non-robust fit.
#'
#' @param custfun A custom function combining the residuals \code{resids} and
#' the hat values \code{hats}. The result of the function must be numeric, 
#' see \link{GCV_crit}.
#'
#' @details Function \code{custfun} has two arguments 
#' corresponding to \code{resids} and \code{hats}. The output of the function 
#' must be numeric. 
#' 
#' @return A named numerical vector of values. The length of the vector depends
#' on the input. The vector contains the values:
#' \itemize{
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
#' }
#'
#' @examples
#' n = 50      # sample size
#' p = 10      # dimension of predictors
#' Z = matrix(rnorm(n*p),ncol=p) # design matrix
#' Y = Z[,1]   # response vector
#' lambda = 1  # tuning parameter for penalization
#' H = diag(p) # penalty matrix
#' type = "absolute" # absolute loss
#' 
#' # Run with the IRLS procedure
#' res = IRLS(Z, Y, lambda, H, type)
#' with(res,GCV_crit(resids,hat_values))
#' with(res,GCV_crit(resids,hat_values,custfun = function(r,h) 
#'     sum(r^2)))
#'     
#' w = rep(1/n,n)
#' GCV_location(lambda,Z,Y,H,type,w,custfun = function(r,h) sum(r^2))

GCV_location <- function(lambda, Z, Y, H, type, w, vrs="C", 
                         method="IRLS",
                         custfun=NULL, 
                         resids.in = rep(1,length(Y)),
                         toler=1e-7, imax=1000){
  
  method = match.arg(method,c("IRLS", "ridge", "HuberQp"))
  type = match.arg(type,c("square","absolute","Huber","logistic"))
  if(method=="ridge" & type!="square") 
    stop("method 'ridge' available only for type 'square'.")
  if(method=="HuberQp" & type!="huber") 
    stop("method 'HuberQp' available only for type 'huber'.")
  
  # Generalized cross-validation
  ncv = 6
  if(!is.null(custfun)) ncv = 7
  vrs = match.arg(vrs, c("C", "R"))
  
  if(method=="IRLS"){
    fit.r <- IRLS(Z, Y, lambda, H, type=type, w=w, vrs=vrs, 
                  resids.in = resids.in,
                  toler=toler, imax=imax)
    GCV.scores <- GCV_crit(fit.r$resids,fit.r$hat_values,
                           custfun=custfun)
    return(c(GCV.scores, fit.r$converged, fit.r$ic))
  }
  if(method=="ridge"){
    fit.r <- ridge(Z, Y, lambda, H, w=w, vrs=vrs)
    GCV.scores <- GCV_crit(fit.r$resids,fit.r$hat_values,custfun=custfun)
    return(c(GCV.scores, 1, 0))
  }
  if(method=="HuberQp"){
    fit.r <- HuberQp(Z, Y, lambda, H, w=w, vrs=vrs)
    GCV.scores <- GCV_crit(fit.r$resids,fit.r$hat_values,custfun=custfun)
    return(c(GCV.scores, 1, 0))
  }
}

#' Cross-validation for Ridge Regression
#'
#' Provides the cross-validation indices from \link{GCV_crit} in
#' conjunction with the \link{ridge} function directly as an argument of the
#' parameter \code{lambda}. 
#'
#' @param lambda A candidate parameter value; non-negative real number.
#'
#' @param Z Data matrix of dimension \code{n}-times-\code{p}, where \code{n} is
#' the number of observations, \code{p} is the dimension.
#'
#' @param Y Vector of responses of length \code{n}.
#'
#' @param H Penalty matrix of size \code{p}-times-\code{p} that
#' is used inside the quadratic term for penalizing estimated parameters.
#' 
#' @param vrs Version of the algorhitm to be used in function \link{ridge}; 
#' either \code{vrs="C"} for the \code{C++} version, or \code{vrs="R"} for the 
#' \code{R} version. Both should give (nearly) identical results, see 
#' \link{ridge}.
#'
#' @param custfun A custom function combining the residuals \code{resids}, the 
#' hat values \code{hats}. The result of the function must be numeric, see
#' \link{GCV_crit}.
#'
#' @details Function \code{custfun} has two arguments 
#' corresponding to \code{resids} and \code{hats}. The output of the function 
#' must be numeric. Additional parameters passed to \code{ridge} may be passed
#' to the function.
#'
#' @return A named numerical vector of values. The length of the vector depends
#' on the input. The vector contains the values:
#' \itemize{
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
#' }
#'
#' @examples
#' n = 50      # sample size
#' p = 10      # dimension of predictors
#' Z = matrix(rnorm(n*p),ncol=p) # design matrix
#' Y = Z[,1]   # response vector
#' lambda = 1  # tuning parameter for penalization
#' H = diag(p) # penalty matrix
#'     
#' # Run with the ridge function
#' res = ridge(Z, Y, lambda, H)
#' with(res,GCV_crit(resids,hat_values))
#' with(res,GCV_crit(resids,hat_values,custfun = function(r,h) 
#'     sum((r/(1-h))^2)))
#'     
#' GCV_ridge(lambda,Z,Y,H,custfun = function(r,h) sum((r/(1-h))^2))

GCV_ridge <- function(lambda,Z,Y,H,vrs="C",custfun=NULL){
  # Generalized cross-validation for ridge
  vrs = match.arg(vrs, c("C", "R"))
  fit.r <- ridge(Z,Y,lambda,H,vrs=vrs)
  GCV.scores <- GCV_crit(fit.r$resids, fit.r$hat_values, custfun=custfun)
  return(GCV.scores)
}

#' Cross-validation for Huber Regression
#'
#' Provides the cross-validation indices from \link{GCV_crit} in
#' conjunction with the \link{HuberQp} function directly as an argument of the
#' parameter \code{lambda}. 
#'
#' @param lambda A candidate parameter value; non-negative real number.
#'
#' @param Z Data matrix of dimension \code{n}-times-\code{p}, where \code{n} is
#' the number of observations, \code{p} is the dimension.
#'
#' @param Y Vector of responses of length \code{n}.
#'
#' @param H Penalty matrix of size \code{p}-times-\code{p} that
#' is used inside the quadratic term for penalizing estimated parameters.
#' 
#' @param vrs Version of the algorhitm to be used in function \link{HuberQp}; 
#' either \code{vrs="C"} for the \code{C++} version, or \code{vrs="R"} for the 
#' \code{R} version. Both should give (nearly) identical results, see 
#' \link{HuberQp}.
#'
#' @param custfun A custom function combining the residuals \code{resids}, the 
#' hat values \code{hats}. The result of the function must be numeric, see
#' \link{GCV_crit}.
#'
#' @details Function \code{custfun} has two arguments 
#' corresponding to \code{resids} and \code{hats}. The output of the function 
#' must be numeric. Additional parameters passed to \code{HuberQp} may be passed
#' to the function.
#'
#' @return A named numerical vector of values. The length of the vector depends
#' on the input. The vector contains the values:
#' \itemize{
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
#' }
#'
#' @examples
#' n = 50      # sample size
#' p = 10      # dimension of predictors
#' Z = matrix(rnorm(n*p),ncol=p) # design matrix
#' Y = Z[,1]   # response vector
#' lambda = 1  # tuning parameter for penalization
#' H = diag(p) # penalty matrix
#'     
#' # Run with the ridge function
#' res = ridge(Z, Y, lambda, H)
#' with(res,GCV_crit(resids,hat_values))
#' with(res,GCV_crit(resids,hat_values,custfun = function(r,h) 
#'     sum((r/(1-h))^2)))
#'     
#' GCV_ridge(lambda,Z,Y,H,custfun = function(r,h) sum((r/(1-h))^2))

GCV_HuberQp <- function(lambda, Z, Y, H, vrs="C", custfun=NULL){
  # Generalized cross-validation for ridge
  vrs = match.arg(vrs, c("C", "R"))
  fit.r <- HuberQp(Z,Y,lambda,H,vrs=vrs)
  GCV.scores <- GCV_crit(fit.r$resids, fit.r$hat_values, custfun=custfun)
  return(GCV.scores)
}

#' Criteria used for cross-validation and for tuning parameter lambda
#'
#' Several criteria commonly used for selection of the tuning parameter 
#' \code{lambda} in functions \link{IRLS} and \link{ridge}. 
#'
#' @param resids A vector of residuals of length \code{n}.
#'
#' @param hats A vector of hat values of length \code{n}.
#' 
#' @param custfun A custom function combining the residuals \code{resids} and 
#' the hat values \code{hats}. The 
#' result of the function must be numeric.
#'
#' @details Function \code{custfun} has two arguments 
#' corresponding to \code{resids} and \code{hats}. The output of the function 
#' must be numeric.
#'
#' @return A named numerical vector of values. The length of the vector depends
#' on the input. The vector contains (some of) the values:
#' \itemize{
#'  \item{"AIC"}{ Akaike's information criterion given by 
#'  \code{mean(resids^2)+2*mean(hats)}, where \code{n} is the length of
#'  both \code{resids} and \code{hats}.}
#'  \item{"GCV"}{ Leave-one-out cross-validation criterion given by
#'  \code{mean((resids^2)/((1-hats)^2))}.}
#'  \item{"GCV(tr)"}{ Modified leave-one-out cross-validation criterion 
#'  given by \code{mean((resids^2)/((1-mean(hats))^2))}.}
#'  \item{"BIC"}{ Bayes information criterion given by 
#'  \code{mean(resids^2)+log(n)*mean(hats)}.}
#'  \item{"rGCV"}{ A robust version of \code{GCV} where mean is replaced
#'  by a robust M-estimator of scale of \code{resids/(1-hats)}, see 
#'  \link[robustbase]{scaleTau2} for details.}
#'  \item{"rGCV(tr)"}{ Modified version of a \code{rGCV} given by 
#'  a robust M-estimator of scale of \code{resids/(1-mean(hats))}.}
#'  \item{"custom"}{ The custom criterion given by function \code{custfun}. 
#'  Works only if \code{custfun} is part of the input.}
#' }
#'
#' @examples
#' n = 50      # sample size
#' p = 10      # dimension of predictors
#' Z = matrix(rnorm(n*p),ncol=p) # design matrix
#' Y = Z[,1]   # response vector
#' lambda = 1  # tuning parameter for penalization
#' H = diag(p) # penalty matrix
#' type = "absolute" # absolute loss
#' 
#' # Run with the IRLS procedure
#' res = IRLS(Z, Y, lambda, H, type)
#' with(res,GCV_crit(resids,hat_values))
#' with(res,GCV_crit(resids,hat_values,custfun = function(r,h) 
#'     sum(r^2)))
#'     
#' # Run with the ridge function
#' res = ridge(Z, Y, lambda, H)
#' with(res,GCV_crit(resids,hat_values))
#' with(res,GCV_crit(resids,hat_values,custfun = function(r,h) 
#'     sum((r/(1-h))^2)))

GCV_crit = function(resids, hats, custfun=NULL){
  n = length(resids)
  GCVs = rep(NA,7)
  names(GCVs) = c("AIC","GCV","GCV(tr)","BIC","rGCV","rGCV(tr)",
                  "custom")
  GCVs[1:6] = c(
    mean(resids^2)+2*mean(hats),
    mean((resids^2)/((1-hats)^2)),
    mean((resids^2)/((1-mean(hats))^2)),
    mean(resids^2)+log(n)*mean(hats),
    robustbase::scaleTau2(resids/(1-hats), c2 = 5),
    robustbase::scaleTau2(resids/(1-mean(hats)), c2 = 5)
  )
  if(is.null(custfun)){
    return(GCVs[1:6])
  } else {
    GCVs[7] = custfun(resids,hats)
    return(GCVs) 
  }
}

#' k-fold cross-validation for the IRLS procedure
#'
#' Provides the k-fold version of the cross-validation procedure in
#' conjunction with the \link{IRLS} function directly as an argument of the
#' parameter \code{lambda}. 
#'
#' @param lambda A candidate parameter value; non-negative real number.
#'
#' @param Z Data matrix of dimension \code{n}-times-\code{p}, where \code{n} is
#' the number of observations, \code{p} is the dimension.
#'
#' @param Y Vector of responses of length \code{n}.
#'
#' @param H Penalty matrix of size \code{p}-times-\code{p} that
#' is used inside the quadratic term for penalizing estimated parameters.
#' 
#' @param type The type of the loss function used in the minimization problem.
#' Accepted are \code{type="absolute"} for the absolute loss \code{rho(t)=|t|}; 
#' \code{type="square"} for the square loss \code{rho(t)=t^2}; 
#' \code{type="Huber"} for the Huber loss \code{rho(t)=t^2/2} if 
#' \code{|t|<tuning} and \code{rho(t)=tuning*(|t|-tuning/2)} otherwise; and 
#' \code{type="logistic"} for the logistic loss 
#' \code{rho(t)=2*t + 4*log(1+exp(-t))-4*log(2)}.
#'
#' @param k Number of folds to consider, positive integer. By default
#' set to 5.
#' 
#' @param vrs Version of the algorhitm to be used in function \link{IRLS}; 
#' either \code{vrs="C"} for the \code{C++} version, or \code{vrs="R"} for the 
#' \code{R} version. Both should give (nearly) identical results, see 
#' \link{IRLS}.
#'
#' @details Cross-validation based on the median of squared residuals obtained
#' in each fold of the data.
#'
#' @return A numerical of medians of squared residuals for each fold of the 
#' data. A numerical vector of \code{k} non-negative values.
#'
#' @examples
#' n = 50      # sample size
#' p = 10      # dimension of predictors
#' Z = matrix(rnorm(n*p),ncol=p) # design matrix
#' Y = Z[,1]   # response vector
#' lambda = 1  # tuning parameter for penalization
#' H = diag(p) # penalty matrix
#' type = "absolute" # absolute loss
#' 
#' kCV(lambda,Z,Y,H,type,k=5)

kCV = function(lambda,Z,Y,H,type,k=5,vrs="C"){
  # k-fold cross-validation
  vrs = match.arg(vrs, c("C", "R"))
  n = length(Y)
  if(n<k) stop("For k-fold cross-validation we need n>=k.")
  #
  nk = floor(n/k) # number of functions in each batch but the last
  batchi = c(rep(1:k,each=nk),rep(k,n-nk*k)) 
  # batch index for each function
  crit = rep(NA,k)
  for(ki in 1:k){
    bi = (batchi==ki)
    Ze = Z[!bi,] # Z for estimation
    Ye = Y[!bi]  # Y for estimation
    Zt = Z[bi,]  # Z for testing
    Yt = Y[bi]   # Y for testing
    fit.r <- IRLS(Ze,Ye,lambda,H,type,vrs=vrs)
    rest = Yt - Zt%*%fit.r$theta_hat # residuals for the testing part
    # if(fit.r$converged==0) rest = rep(Inf,length(Yt)) # If IRLS did not converge
    crit[ki] = median(rest^2) # robustbase::scaleTau2(rest^2, c2 = 5)
    # criterion for ki-th batch
  }
  return(crit) # criterion for each batch
} 