#' Iteratively Reweighted Least Squares for robust functional regression
#'
#' Iteratively Reweighted Least Squares (IRLS) algorithm that is used to 
#' estimate a vector of regression parameters in a (possibly robust and 
#' penalized) linear regression model. Weights can be supplied as well.
#'
#' @param Z Data matrix of dimension \code{n}-times-\code{p}, where \code{n} is
#' the number of observations, \code{p} is the dimension.
#'
#' @param Y Vector of responses of length \code{n}.
#'
#' @param lambda Tuning parameter, a non-negative real number.
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
#' @param sc Scale parameter to be used in the IRLS. By default \code{sc=1}, 
#' that is no scaling is performed.
#'
#' @param resids.in Initialization of the vector of residuals used to launch 
#' the IRLS algorithm.
#' 
#' @param tuning A non-negative tuning constant for the absolute loss function 
#' (that is, \code{type="absolute"}). For \code{tuning = 0} the standard 
#' absolute loss \code{rho(t) = |t|} is used. For \code{tuning > 0}, the Huber 
#' loss is used, that is \code{rho(t)} is quadratic for \code{|t|<tuning} and 
#' linear for \code{|t|>=tuning}. The function is chosen so that \code{rho} 
#' is always continuously differentiable.
#' 
#' @param toler A small positive constant specifying the tolerance level for 
#' terminating the algorithm. The prcedure stops if the maximum absolute 
#' distance between the residuals in the previous iteration and the new 
#' residuals drops below \code{toler}.
#' 
#' @param imax Maximum number of allowed iterations of IRLS. 
#'
#' @param vrs Version of the algorhitm to be used. The program is prepared in
#' two versions: i) \code{vrs="C"} calls the \code{C++} version of the 
#' algorithm, programmed within the \code{RCppArmadillo} framework for
#' manipulating matrices. This is typically the fastest version. 
#' ii) \code{vrs="R"} calls the \code{R} version. The two versions may 
#' give slightly different results due to the differences in evaluating inverse
#' matrices. With \code{vrs="C"} one uses the function \code{solve} directly
#' from \code{Armadillo} library in \code{C++}; with \code{vrs="R"} the 
#' standard function \code{solve} from \code{R} package \code{base} is used 
#' with the option \code{tol = toler_solve}.
#'
#' @param toler_solve A small positive constant to be passed to function
#' \link[base]{solve} as argument \code{tol}. Used to handle numerically 
#' singular matrices whose inverses need to be approximated. By default set to
#' 1e-35.
#'
#' @details Especially for extremely small values of \code{lambda}, numerically
#' singular matrices must be inverted in the procedure. This may cause numerical
#' instabilites, and is the main cause for differences in results when using
#' \code{vrs="C"} and \code{vrs="R"}. In case when IRLS does not converge within
#' \code{imax} iterations, a warning is given.
#'
#' @return A list composed of:
#' \itemize{
#'  \item{"theta_hat"}{ A numerical matrix of size \code{p}-times-\code{1} of 
#'  estimated regression coefficients.}
#'  \item{"converged"}{ Indicator whether the IRLS procedure successfully 
#'  converged. Takes value 1 if IRLS converged, 0 otherwise.}
#'  \item{"ic"}{ Number of iterations needed to reach convergence. If 
#'  \code{converged=0}, always \code{ic=imax}.}
#'  \item{"resids"}{ A numerical vector of length \code{n} containing the final
#'  set of residuals in the fit of \code{Y} on \code{Z}.}
#'  \item{"hat_values"}{ Diagonal terms of the (possibly penalized) hat matrix of
#'  the form \code{Z*solve(t(Z)*W*Z+n*lambda*H)*t(Z)*W}, where \code{W} 
#'  is the diagonal weight matrix in the final iteration of IRLS.}
#'  \item{"last_check"}{ The final maximum absolute difference between 
#'  \code{resids} and the residuals from the previous iteration. We have 
#'  \code{resids < toler} if and only if the IRLS converged (that is, 
#'  \code{converged=1}).}
#'  \item{"weights"}{ The vector of weights given to the observations in the 
#'  final iteration of IRLS. For squared loss (\code{type="square"}) this gives 
#'  a vector whose all elements are 2.}
#'  \item{"fitted"}{ Fitted values in the model. A vector of length \code{n} 
#'  corresponding to the fits of \code{Y}.}
#' }
#'
#' @references
#' Ioannis Kalogridis and Stanislav Nagy. (2023). Robust functional regression 
#' with discretely sampled predictors. 
#' \emph{Under review}.
#'
#' Peter. J. Huber. (1981). Robust Statistics, \emph{New York: John Wiley.}
#'
#' @seealso \link{ridge} for a faster (non-robust) version of this
#' function with \code{type="square"}.
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
#' # Run the two versions of the IRLS procedure
#' res_C = IRLS(Z, Y, lambda, H, type, vrs="C")
#' res_R = IRLS(Z, Y, lambda, H, type, vrs="R")
#' # Check whether both versions converged after the same number of iterations
#' res_C$ic
#' res_R$ic
#' # Check the maximum absolute difference between the results
#' max(abs(res_C$theta_hat-res_R$theta_hat))
#' # Visualise the difference between the results
#' plot(res_C$theta_hat ~ res_R$theta_hat)

IRLS = function(Z, Y, lambda, H, type, w=NULL, sc = 1, 
                resids.in = rep(1,length(Y)), 
                tuning=NULL, toler=1e-7, imax=1000, vrs="C", 
                toler_solve=1e-35){
  
  IRLS_R <- function(Z, Y, lambda, H, type, w=NULL, sc, resids.in, 
                      tuning, toler, imax, toler_solve){
    n = length(Y)
    if(is.null(w)) w = rep(1/n,n)
    if(length(w)!=n) stop("Weights w must have the same length as Y.")
    ic = 0
    istop = 0
    clim = c(.5, 1, .5*tuning, 1)[type] # tail behavior constant for psiw
    while(istop == 0 & ic < imax){
      ic = ic + 1
      #
      Wdiag = c(w*psiw(resids.in/sc,type,tuning)/sc^2)
      naind = (is.na(Wdiag)) & (abs(resids.in)<tuning)
      Wdiag[naind] = w[naind]*1 # division 0/0
      naind = (is.na(Wdiag))
      Wdiag[naind] = w[naind]*clim/abs(resids.in*sc) 
      # if resids are too small exp(-resids) ~ Inf but in the limit always
      # psiw(t) ~ clim/abs(t)
      Z.s = scale(t(Z), center = FALSE, scale = 1/Wdiag);
      Z1 = Z.s%*%Z + lambda*H
      theta_new = solve(Z1, Z.s%*%Y, tol=toler_solve)
      resids1 <- c(Y - Z%*%theta_new)
      check = max(abs(resids1-resids.in))/sc 
      if(check < toler){istop=1}
      resids.in <- resids1
    }
    if(istop==0) warning(
      paste0("log(lambda) ",sprintf("%.3f",log(lambda)),
             ": Estimator did not converge."))
    hat.values = diag(hatm<-Z%*%solve(Z1, Z.s, tol=toler_solve))
    #
    return(list(theta_hat = theta_new,
                converged=istop,
                ic=ic,
                resids = resids.in, 
                hat_values = hat.values,
                last_check = check,
                weights = 2*Wdiag,
                fitted = hatm%*%Y))
  }
  
  n = length(Y)
  if(is.null(w)) w = rep(1/n,n)
  vrs = match.arg(vrs,c("C","R"))
  type = match.arg(type,c("square","absolute","Huber","logistic"))
  type = switch(type, absolute = 1, square = 2, Huber = 3, logistic = 4)
  if(type==3 & is.null(tuning)){
    tuning = 1.345
    # warning("Huber loss, setting constant to default 1.345.")
  }
  if(type==1 & is.null(tuning)){
    tuning = 1/100
    # warning("absolute loss, setting tuning to default 1/100.")    
  }
  if(is.null(tuning)) tuning = 1/100 # tuning for logistic regression
  #
  if(nrow(Z)!=length(Y)) 
    stop("Number of rows of Z must equal the lenght of Y.")
  if(nrow(H)!=ncol(Z))
    stop("H must be a square matrix with the same number of columns as Z.")
  if(ncol(H)!=ncol(Z))
    stop("H must be a square matrix with the same number of columns as Z.")
  if(tuning<0) stop("tuning must be a non-negative number.")
  if(lambda<0) stop("lambda must be a non-negative number.")
  if(sc<=0) stop("Scale estimator must be strictly positive.")
  if(vrs=="C"){
    rs = tryCatch(
      error = function(cnd){
        warning(paste0("Solve in C++ crashed, switching to R version, ",cnd))
        IRLS_R(Z, Y, lambda, H, type, w, sc, resids.in, 
               tuning, toler, imax, toler_solve)
      }, {
        IRLSC(Z, Y, lambda, H, type, w, sc, resids.in, 
              tuning, toler, imax)
      })
    return(rs)
  }
  
  if(vrs=="R") return(IRLS_R(Z, Y, lambda, H, type, w, sc, resids.in, 
                             tuning, toler, imax, toler_solve))
}

#' Fast Ridge Regression with given penalty matrix
#'
#' A (weighted) ridge regression estimator with a specified penalty matrix
#' in a linear regression model. The solution corresponds to the 
#' result of function \link{IRLS} with \code{type="square"}.
#'
#' @param Z Data matrix of dimension \code{n}-times-\code{p}, where \code{n} is
#' the number of observations, \code{p} is the dimension.
#'
#' @param Y Vector of responses of length \code{n}.
#'
#' @param lambda Tuning parameter, a non-negative real number.
#'
#' @param H Penalty matrix of size \code{p}-times-\code{p} that
#' is used inside the quadratic term for penalizing estimated parameters.
#' 
#' @param w Vector of length \code{n} of weights attached to the elements of 
#' \code{Y}. If \code{w=NULL} (default), a constant vector with values 
#' \code{1/n} is used.
#' 
#' @param vrs Version of the algorhitm to be used. The program is prepared in
#' two versions: i) \code{vrs="C"} calls the \code{C++} version of the 
#' algorithm, programmed within the \code{RCppArmadillo} framework for
#' manipulating matrices. This is typically the fastest version. 
#' ii) \code{vrs="R"} calls the \code{R} version. The two versions may 
#' give slightly different results due to the differences in evaluating inverse
#' matrices. With \code{vrs="C"} one uses the function \code{solve} directly
#' from \code{Armadillo} library in \code{C++}; with \code{vrs="R"} the 
#' standard function \code{solve} from \code{R} package \code{base} is used 
#' with the option \code{tol = toler_solve}.
#'
#' @param toler_solve A small positive constant to be passed to function
#' \link[base]{solve} as argument \code{tol}. Used to handle numerically 
#' singular matrices whose inverses need to be approximated. By default set to
#' 1e-35.
#'
#' @details Especially for extremely small values of \code{lambda}, numerically
#' singular matrices must be inverted in the procedure. This may cause numerical
#' instabilites, and is the main cause for differences in results when using
#' \code{vrs="C"} and \code{vrs="R"}. This function is equivalent with 
#' \link{IRLS} when used with the square loss \code{type="square"}, but faster
#' and more stable as it does not perform the iterative algorithm. Instead, it 
#' computes the estimator directly.
#'
#' @return A list composed of:
#' \itemize{
#'  \item{"theta_hat"}{ A numerical matrix of size \code{p}-times-\code{1} of 
#'  estimated regression coefficients.}
#'  \item{"resids"}{ A numerical vecotor of length \code{n} containing the final
#'  set of residuals in the fit of \code{Y} on \code{Z}.}
#'  \item{"hat_values"}{ Diagonal terms of the (penalized) hat matrix of
#'  the form \code{Z*solve(t(Z)*Z + n*lambda*H)*t(Z)}.}
#'  \item{"fitted"}{ Fitted values in the model. A vector of length \code{n} 
#'  correponding to the fits of \code{Y}.}
#' }
#' 
#' @seealso \link{IRLS} for a robust version of this
#' function.
#'
#'
#' @examples
#' n = 50      # sample size
#' p = 10      # dimension of predictors
#' Z = matrix(rnorm(n*p),ncol=p) # design matrix
#' Y = Z[,1]   # response vector
#' lambda = 1  # tuning parameter for penalization
#' H = diag(p) # penalty matrix
#' 
#' res_C = ridge(Z, Y, lambda, H, vrs="C")
#' res_R = ridge(Z, Y, lambda, H, vrs="R")
#' # Check the maximum absolute difference between the results
#' max(abs(res_C$theta_hat-res_R$theta_hat))
#' # Visualise the difference between the results
#' plot(res_C$theta_hat ~ res_R$theta_hat)
#' 
#' # Compare the output with function IRLS
#' res_IRLS = IRLS(Z, Y, lambda, H, type="square")
#' max(abs(res_C$theta_hat-res_IRLS$theta_hat))

ridge = function(Z, Y, lambda, H, w=NULL, vrs="C", toler_solve=1e-35){
  
  n = length(Y)
  if(is.null(w)) w = rep(1/n,n)
  vrs = match.arg(vrs,c("C","R"))
  if(nrow(Z)!=length(Y)) 
    stop("Number of rows of Z must equal the lenght of Y.")
  if(nrow(H)!=ncol(Z))
    stop("H must be a square matrix with the same number of columns as Z.")
  if(ncol(H)!=ncol(Z))
    stop("H must be a square matrix with the same number of columns as Z.")
  if(lambda<0) stop("lambda must be a non-negative number.")
  # if(sc<=0) stop("Scale estimator must be strictly positive.")
  if(vrs=="C") return(ridgeC(Z, Y, lambda, H, w))
  
  ridge_R <- function(Z, Y, lambda, H, w=NULL, toler_solve){
    n = length(Y)
    th = solve(t(Z)%*%diag(w)%*%Z+lambda*H,t(Z)%*%diag(w)%*%Y,
               tol=toler_solve)
    fitted =  Z%*%th
    hat = Z%*%solve(t(Z)%*%diag(w)%*%Z+lambda*H,t(Z)%*%diag(w),tol=toler_solve)
    resid = Y - fitted
    return(list(theta_hat = th,
                resids = resid, 
                hat_values = diag(hat),
                fitted = fitted))
  }
  
  if(vrs=="R") return(ridge_R(Z, Y, lambda, H, w, toler_solve))
}

#' Fast Regression with Huber Penalty - Quadratic programming solution
#'
#' A Huber regression estimator with a specified penalty matrix
#' in a linear functional regression model. The solution corresponds to the 
#' result of function \link{IRLS} with \code{type="Huber"}.
#'
#' @param Z Data matrix of dimension \code{n}-times-\code{p}, where \code{n} is
#' the number of observations, \code{p} is the dimension.
#'
#' @param Y Vector of responses of length \code{n}.
#'
#' @param lambda Tuning parameter, a non-negative real number.
#'
#' @param H Penalty matrix of size \code{p}-times-\code{p} that
#' is used inside the quadratic term for penalizing estimated parameters.
#' 
#' @param w Vector of length \code{n} of weights attached to the elements of 
#' \code{Y}. If \code{w=NULL} (default), a constant vector with values 
#' \code{1/n} is used.
#' 
#' @param vrs Version of the algorhitm to be used. The program is prepared in
#' two versions: i) \code{vrs="C"} calls the \code{C++} version of the 
#' algorithm, programmed within the \code{RCppArmadillo} framework for
#' manipulating matrices. This is typically the fastest version, especially if p >> n. 
#' ii) \code{vrs="R"} calls the \code{R} version. The two versions may 
#' give slightly different results due to the different tolerances used for the solution
#' of the qudratic problem.
#'
#' @param toler_solve Unused at the moment
#'
#' @details This function is equivalent with 
#' \link{IRLS} when used with the square loss \code{type="Huber"}, but faster 
#' and more stable as it does not perform the iterative algorithm. Note that this function
#' is faster than \link{IRLS} only if the number of predictors p is larger than the number of functional observations n.
#' On the contrary, if n > p, \link{IRLS} is faster, even though less accurate. When p > n, the  \code{vrs="C"} version is
#' faster than the  \code{vrs="R"} one.
#'
#' @return A list composed of:
#' \itemize{
#'  \item{"theta_hat"}{ A numerical matrix of size \code{p}-times-\code{1} of 
#'  estimated regression coefficients.}
#'  \item{"resids"}{ A numerical vecotor of length \code{n} containing the final
#'  set of residuals in the fit of \code{Y} on \code{Z}.}
#'  \item{"hat_values"}{ Diagonal terms of the (penalized) hat matrix of
#'  the form \code{Z*solve(t(Z)*Z + n*lambda*H)*t(Z)}.}
#'  \item{"fitted"}{ Fitted values in the model. A vector of length \code{n} 
#'  correponding to the fits of \code{Y}.}
#' }
#' 
#' @seealso \link{IRLS} for an iterative version of this
#' function.
#'
#'
#' @examples
#' n = 50      # sample size
#' p = 10      # dimension of predictors
#' Z = matrix(rnorm(n*p),ncol=p) # design matrix
#' Y = Z[,1]   # response vector
#' lambda = 1  # tuning parameter for penalization
#' H = diag(p) # penalty matrix
#' 
#' res_C = HuberQp(Z, Y, lambda, H, vrs="C")
#' res_R = HuberQp(Z, Y, lambda, H, vrs="R")
#' # Check the maximum absolute difference between the results
#' max(abs(res_C$theta_hat-res_R$theta_hat))
#' # Visualize the difference between the results
#' plot(res_C$theta_hat ~ res_R$theta_hat)
#' 
#' # Compare the output with function IRLS
#' res_IRLS = IRLS(Z, Y, lambda, H, type="Huber")
#' max(abs(res_C$theta_hat-res_IRLS$theta_hat))

HuberQp = function(Z, Y, lambda, H, w=NULL, vrs="C", toler_solve=1e-35){
  n = length(Y)
  if(is.null(w)) w = rep(1/n,n)
  # Scale the weights i order to obtain the exact result of IRLS (the entire equation is scaled)
  w = w*n 
  vrs = match.arg(vrs,c("C","R"))
  if(nrow(Z)!=length(Y)) 
    stop("Number of rows of Z must equal the lenght of Y.")
  if(nrow(H)!=ncol(Z))
    stop("H must be a square matrix with the same number of columns as Z.")
  if(ncol(H)!=ncol(Z))
    stop("H must be a square matrix with the same number of columns as Z.")
  if(lambda<0) stop("lambda must be a non-negative number.")
  # if(sc<=0) stop("Scale estimator must be strictly positive.")
  
  ### Solve the problem (C++ version)
  if(vrs=="C"){
    return(HuberQpC(Z, Y, 2*n*lambda*H, w))
  }
  
  huber_qp_osqp_penalized <- function(X, y, H, w, delta = 1.345) {
    # Input validation
    if (!is.matrix(X)) stop("X must be a matrix")
    if (!is.vector(y)) stop("y must be a vector")
    if (any(is.na(X)) || any(is.na(y))) stop("X and y must not contain NA values")
    n <- nrow(X)
    p <- ncol(X)
    if (length(y) != n) stop("Length of y must equal number of rows in X")
    if (!is.matrix(H)) stop("H must be a matrix")
    if (nrow(H) != p || ncol(H) != p) stop("H must be p x p, where p is ncol(X)")
    if (!isSymmetric(H)) stop("H must be symmetric")
    
    ### Applying the weights (Check that we obtain the same results with IRLS) #### CHECK
    ### AAA
    X1 = diag(w)%*%X
    y1 = w*y # diag(w) %*% Y if Y is a columns vector
    
    nvars <- p + 2 * n  # variables: beta (p), a (n), t (n)
    ncons <- 3 * n      # constraints: 2n for t_i >= |y_i - X_i^T beta - a_i|, n for t_i >= 0
    
    # Quadratic term: P as sparse block-diagonal with H, I_n, and 0
    P <- bdiag(H + 1e-6 * diag(p), Diagonal(n, 1 + 1e-6), Diagonal(n, 1e-6))
    
    # Linear term: qvec = delta for t_i's
    qvec <- rep(0, nvars)
    qvec[(p + n + 1):nvars] <- delta
    
    # Constraints: A %*% x >= lvec
    Amat <- matrix(0, ncons, nvars)  # ncons rows x nvars cols
    bvec <- rep(0, ncons)
    
    for (i in 1:n) {
      # Constraint 1: t_i + X_i^T beta + a_i >= y_i
      row1 <- (i - 1) * 2 + 1
      Amat[row1, 1:p] <- X1[i, ]
      Amat[row1, p + i] <- 1
      Amat[row1, p + n + i] <- 1
      bvec[row1] <- y1[i]
      
      # Constraint 2: t_i - X_i^T beta - a_i >= -y_i
      row2 <- (i - 1) * 2 + 2
      Amat[row2, 1:p] <- -X1[i, ]
      Amat[row2, p + i] <- -1
      Amat[row2, p + n + i] <- 1
      bvec[row2] <- -y1[i]
      
      # Constraint 3: t_i >= 0
      row3 <- 2 * n + i
      Amat[row3, p + n + i] <- 1
      bvec[row3] <- 0
    }
    
    # Convert Amat to sparse matrix
    A <- as(Amat, "dgCMatrix")
    
    # Set constraint bounds
    lvec <- bvec
    uvec <- rep(Inf, ncons)
    
    # Solve QP using osqp (This is already a C++ routine)
    settings <- osqpSettings(verbose = FALSE)
    model <- osqp(P = P, q = qvec, A = A, l = lvec, u = uvec, pars = settings)
    sol <- model$Solve()
    
    # Extract beta
    beta <- sol$x[1:p]
    
    # Compute Hat matrix 
    hat = X%*%solve(t(X)%*%diag(w)%*%X+lambda*H,t(X)%*%diag(w),tol=toler_solve)
    
    # Compute residuals
    fitted = X%*%beta
    resid = y - fitted
    
    return(list(theta_hat = beta,
                resids = resid,
                hat_values = diag(hat),
                fitted = fitted))
  }
  
  ### Solve using the R version (relies on the C++ oqsp routine)
  if(vrs=="R") return(huber_qp_osqp_penalized(Z, Y, 2*n*lambda*H, w))
}

#' Weight function for the IRLS algorithm
#'
#' Returns a vector of weights given by \code{psi(t)/(2*t)}, where 
#' \code{psi} is the derivative of the loss function \code{rho}. 
#'
#' @param t Vector of input values of length \code{n}.
#'
#' @param type Integer code for the type of loss function. Accepted are
#' \code{type=1} for the absolute loss \code{rho(t)=|t|}; \code{type=2} for
#' the square loss \code{rho(t)=t^2}; \code{type=3} for the Huber loss
#' \code{rho(t)=t^2/2} if \code{|t|<tuning} and 
#' \code{rho(t)=tuning*(|t|-tuning/2)} otherwise; and \code{type=4} for the
#' logistic loss \code{rho(t)=2*t + 4*log(1+exp(-t))-4*log(2)}.
#'
#' @param tuning Tuning parameter, a non-negative real number. For the absolute
#' this should be a small number that 'smooths' out the numerical effects of
#' the kink of \code{rho} near the origin (by default \code{tuning = 1/100}. 
#' For the Huber loss \code{tuning} is the constant to be used in the function
#' (by default \code{tuning = 1.345}. For \code{type=2} or \code{type=4} this 
#' constant is not used. 
#'
#'
#' @return A numerical vector of values.
#'
#' @examples
#' curve(psiw(x,type=1,tuning=0.1),-5,5) # absolute loss with tuning
#' curve(psiw(x,type=2),-5,5) # square loss
#' curve(psiw(x,type=3),-5,5) # Huber loss
#' curve(psiw(x,type=4),-5,5) # logistic loss

psiw = function(t,type,tuning=NULL){
  # Type is now only the code 1-4
  if(type==3 & is.null(tuning)){
    tuning = 1.345
    # warning("Huber loss, setting constant to default 1.345.")
  }
  if(type==1 & is.null(tuning)){
    tuning = 1/100
    # warning("absolute loss, setting tuning to default 1/100.")    
  }
  if(is.null(tuning)) tuning = 0
  return(psiwC(t,type,tuning))
}

#' Loss functions
#'
#' Several typically used loss functions \code{rho}. 
#'
#' @param t Vector of input values of length \code{n}.
#'
#' @param type Integer code for the type of loss function. Accepted are
#' \code{type="absolute"} for the absolute loss \code{rho(t)=|t|}; 
#' \code{type="square"} for the square loss \code{rho(t)=t^2}; 
#' \code{type="Huber"} for the Huber loss \code{rho(t)=t^2/2} if 
#' \code{|t|<Hk} and \code{rho(t)=Hk*(|t|-Hk/2)} otherwise; and 
#' \code{type="logistic"} for the logistic loss 
#' \code{rho(t)=2*t + 4*log(1+exp(-t))-4*log(2)}.
#'
#' @param Hk Tuning parameter, a non-negative real number, affects only
#' the Huber loss (by default \code{Hk = 1.345}.
#'
#' @return A numerical vector of values.
#'
#' @examples
#' curve(rho(x,type="absolute"),-5,5) # absolute loss with tuning
#' curve(rho(x,type="square"),-5,5) # square loss
#' curve(rho(x,type="Huber"),-5,5) # Huber loss
#' curve(rho(x,type="logistic"),-5,5) # logistic loss

rho = function(t,type,Hk=1.345){          # loss function
  type = match.arg(type,c("square","absolute","Huber","logistic"))
  if(Hk<0) stop("Hk must be a non-negative number.")
  if(type=="absolute") return(abs(t))
  if(type=="square") return(abs(t)^2)
  if(type=="Huber"){
    return((abs(t)<=Hk)*(t^2/2) + (abs(t)>Hk)*(Hk*(abs(t)-Hk/2)))
  }
  if(type=="logistic") return(2*t + 4*log(1+exp(-t))-4*log(2))
}

#' Radial functions for thin-plate splines
#'
#' The function \code{eta_{m,d}} involved in the definition of a 
#' \code{d}-dimensional thin-plate spline of order \code{m}.
#'
#' @param x Vector of input values of length \code{n}.
#' 
#' @param d Dimension of the domain, positive integer.
#' 
#' @param m Order of the spline, positive integer.
#'
#' @return A numerical vector of values.
#'
#' @examples
#' curve(eta(x,d=1,m=3),-5,5) 

eta = function(x,d,m){
  if(d%%2==0){
    x[x==0] = 1 # log(t) will be changed from -Infty to 0
    return((-1)^{m+1+d/2}/(2^{2*m-1}*pi^{d/2}*
                             factorial(m-1)*factorial(m-d/2))*
             x^{2*m-d}*log(x))
  }
  if(d%%2==1) return(gamma(d/2-m)/(2^{2*m}*pi^{d/2}*factorial(m-1))*
                       x^{2*m-d})
}

#' Area of the cells in Voronoi tesselation
#'
#' For one-dimensional or two-dimensional sets of points, 
#' finds the Voronoi tesselation of the data and returns
#' a vector of volumes (lengths or areas) of all the Voronoi 
#' cells.
#'
#' @param x The dataset whose Voronoi tesselation should be considered.
#' \code{p}-times-\code{d} matrix, one row per point in 
#' \code{d}-dimensional space. If \code{d==1}, the elements of \code{x} must
#' be sorted.
#'
#' @param I.method Input method for the complete domain \code{I}
#' where the Voronoi tesselation is evaluated. Takes a value \code{"box"}
#' for \code{I} the smallest axes-aligned box, or \code{"chull"}
#' for \code{I} being a convex hull of points. By default set to 
#' \code{"chull"}. In dimension \code{d=1}, the two methods \code{"box"}
#' and \code{"chull"} are equivalent.
#' 
#' @param I A set of points specifying the complete domain \code{I}.
#' In general, can be a \code{q}-times-\code{d} matrix, where \code{q}
#' is the number of points and \code{d} is the dimension. The matrix \code{I}
#' then specifies the point from which to compute the domain \code{I}: 
#' (a) If \code{method.I=="chull"}, the domain is the convex hull of \code{I}; 
#' (b) If \code{method.I=="box"}, the domain is the smallest axes-aligned box
#' that contains \code{I}. If \code{I} is \code{NULL} (by default), then 
#' \code{I} is taken to be the same as \code{x}. If \code{I.method=="box"}, 
#' \code{I} can be specified also by a pair of real values \code{a<b}, 
#' in which case we take \code{I} to be the axis-aligned sqare \code{[a,b]^d}.
#' 
#' @param scale Should the resulting volumes be scaled to sum to \code{1}? By
#' default set to \code{TRUE}.
#' 
#' @param plot A logical indicator of whether the resulting Voronoi cells should
#' be plotted. By default \code{FALSE}.
#'
#' @return A numerical vector of length \code{p} of volumes of the Voronoi cells
#' corresponding to the rows of \code{x}. If \code{scale==TRUE} (default), this
#' vector is scaled so that the sum of its elements is \code{1}. If some of the
#' elements of this vector is numerically zero (that is, some of the Voronoi 
#' cells are empty), a warning is given.
#'
#' @examples
#' p = 50      # number of observed points in domain
#' x = runif(p) # x-coordinates of the points
#' y = runif(p) # y-coordinates of the points
#' 
#' # Voronoi areas in the plane:
#' # (a) I is the smallest box containing the data
#' vorArea(cbind(x,y), I.method="box",plot=TRUE)
#' 
#' # (b) I is the box [0,1]^2
#' vorArea(cbind(x,y),I.method="box",I=c(0,1),plot=TRUE)
#' 
#' # (c) I is the convex hull of the data
#' vorArea(cbind(x,y),I.method="chull",plot=TRUE)
#' 
#' # Voronoi areas in the real line:
#' # (a) I is the smallest interval containing the data
#' vorArea(x = matrix(sort(x),ncol=1),I.method="chull",plot=TRUE)
#' # (b) I = [0,1]
#' vorArea(x = matrix(sort(x),ncol=1),I.method="box",I=c(0,1),plot=TRUE)
#' # (c) I = [-1,1]
#' vorArea(x = matrix(sort(x),ncol=1),I.method="chull",I=matrix(c(-1,1),ncol=1),plot=TRUE)

vorArea = function(x, I.method = "chull", I=NULL, scale = TRUE, plot = FALSE){
  # input:
  # x: p-times-d matrix, one row per time point in R^d
  # I.method: "box" for I being the smallest axes-aligned box
  #           "chull" for I being the convex hull of I
  # I: q-times-d matrix, specifying the observation window I
  #    If I is NULL, then I = x.
  #    If I.method="box", I is a pair of points a<b, and I=[a,b]^d
  # output:
  # vector of length p of areas of the Voronoi cells
  # scaled to have total sum 1
  
  if(!is.matrix(x)) stop("Points in the domain x must be a matrix of dimension p-times-d.")
  d = ncol(x)
  p = nrow(x) 
  if(d>2){
    warning("vorArea currently implemented only for dimensions d<=2, output is only a vector of the same value 1/p.")
    return(rep(1/p,length=p))
  }
  I.method = match.arg(I.method, c("box","chull"))
  if(is.null(I)) I = x
  if(I.method=="box"){
    if(!is.matrix(I)) I = matrix(I,nrow=2,ncol=d)
    rng = apply(I,2,range)
    bnds = lapply(seq_len(ncol(rng)), function(i) rng[,i])
    I = expand.grid(bnds)
    # I = expand.grid(replicate(d, 0:1, simplify=FALSE))  
  }
  #
  if(ncol(I)!=d) stop("The dimensions of x and I in vorArea do not match.")
  #
  if(d==1){
    I = apply(I,2,range)
    if(I[1,1]>=I[2,1]) stop("The points defining the domain I must be ordered.")
    tobs = x[,1]
    if(any(order(tobs)!=1:p))
      stop("For numerical integration, values of tobs must be ordered.")
    mids = c(I[1,1],(tobs[-1]+tobs[-p])/2,I[2,1]) # midpoints between adjacent intervals
    mids[mids<=I[1,1]] = I[1,1]
    mids[mids>=I[2,1]] = I[2,1]
    ret = diff(mids)
    #
    if(plot){
      plot(rep(0,p) ~ tobs, pch = 19, cex=.5)
      abline(h=0,lty=3)
      points(rep(0,length(mids)) ~ mids, pch=3, col=2)
    }
    # if(sum(integral_weights)!=1)
    #  stop("Problem with computing integrating weights.")
    # weights given in the one-dimensional integration to the points tobs
    if(any(ret<=1e-25)) warning("Some Voronoi cells are of zero area, possibly conv(tobs) is not a subset of I.")
    if(scale) ret = ret/sum(ret) # scaling for total sum 1
    return(unname(ret))
  }
  #
  if(d==2){
    # deldir package
    I = I[chull(I),]
    I = list(x=I[,1],y=I[,2])
    #
    tesselation = tryCatch(
      error = function(cnd){
        warning(paste0("vorArea failed, jittering data, ",cnd))
        deldir(x[,1], x[,2])}, {
        deldir(jitter(x[,1]), jitter(x[,2]))
      })
      # deldir(x[,1], x[,2])
    tiles <- tile.list(tesselation, clipp=I)
    # indices of points whose cells are non-empty
    inds = sapply(tiles,function(x) x$ptNum)
    areas = sapply(tiles,function(x) x$area)
    #
    if(plot) plot(tiles, pch = 19)
    #
    ret = rep(0,p)
    ret[inds] = areas
    if(any(ret<=1e-25)) warning("Some Voronoi cells are of zero area, possibly conv(tobs) is not a subset of I.")
    if(scale) ret = ret/sum(ret) # scaling for total sum 1
    return(unname(ret))
  }
}
