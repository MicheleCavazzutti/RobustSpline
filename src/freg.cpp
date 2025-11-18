# define ARMA_WARN_LEVEL 1  // Turns off warnings about inverses of nearly singular matrices
# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <osqp/osqp.h>
#include <vector>
#include <limits>
#include <cstring>


// Stanislav Nagy
// nagy@karlin.mff.cuni.cz
// 31.08.2025

// [[Rcpp::export()]]
arma::vec psiwC (arma::vec t, int type, double tuning){
// function that gives the weights psi(t)/(2*t) (with 2 in the denominator)
// for reweighting in the IRLS algorithm

// type: the type of the loss function used. Encoding:
//       1 for the absolute loss
//       2 for the square loss
//       3 for the Huber loss, in this case tuning is the constant in the loss
//       4 for the logistic loss
  int n = t.n_elem;
  arma::vec res(n);
  // absolute loss, with slight Huberization at constant tuning
  // for numerical reasons
  if(type == 1){
    if(tuning == 0){  // if tuning == 0 => absolute loss without tuning
      for(int i=0; i<n; i++){
        if(t(i)>0) res(i) = 0.5/t(i);
        if(t(i)<0) res(i) = -0.5/t(i);
        if(t(i)==0) res(i) = R_PosInf; 
      }
    } else {
      for(int i=0; i<n; i++){
        if(abs(t(i))<=tuning){
          res(i) = 0.5/tuning;
        } else {
          if(t(i)>0) res(i) = 0.5/t(i);
          if(t(i)<0) res(i) = -0.5/t(i);
        }
      }
    }
    return(res);
  }
  // square loss
  if(type == 2){
    res.ones();
    return(res);
  }
  // Huber loss
  if(type == 3){
    if(tuning == 0){ // tuning == 0 => absolute loss
      for(int i=0; i<n; i++){
        if(t(i)>0) res(i) = 0.5/t(i);
        if(t(i)<0) res(i) = -0.5/t(i);
        if(t(i)==0) res(i) = R_PosInf;
      } 
    } else {
      for(int i=0; i<n; i++){
        if(abs(t(i))<=tuning){
          res(i) = 0.5;
        } else { 
          res(i) = 0.5*tuning/abs(t(i)); 
        }      
      }      
    }
  return(res);      
  }
  // logistic loss
  if(type == 4){
    for(int i=0; i<n; i++){
      if(t(i)==0){
        res(i) = 0.5; 
      } else {
        if(abs(t(i))>=500){
        res(i) = 1/abs(t(i)); // if abs(t) is too large, exp(t) might give
        // numerically Inf, but in limit we have that psiw converges to 1/abs(t)
        } else {
        res(i)=(-2*exp(-t(i))/(1+exp(-t(i)))+1)/t(i);
        }
      }
    }
  return(res);
  }
  // If type is unknown returns zeros
  res.zeros();
  return(res);
}

// [[Rcpp::export()]]
Rcpp::List IRLSC (arma::mat Z, arma::vec Y, double lambda, arma::mat H,
  int type, arma::vec W, double sc, arma::vec residsin, double tuning, double toler, int imax){
  
  int ic = 0, istop = 0;
  int n = Y.n_elem, p1 = Z.n_cols;
  arma::vec w(n);
  arma::mat Zt = trans(Z), ZW(p1,n), Z1(p1,p1); // Zt is of dimension p1-times-n
  arma::mat hat_mat(n,n); // hat matrix for hat values
  arma::vec theta_new(p1), resids1(n), hat_diag(n), fitted(n);
  double check = 0;
  
  while(istop == 0 && ic < imax){
    ic++;
    w = W % psiwC(residsin/sc, type, tuning); // element-wise product
    w = w/(sc*sc); // scaling affects only the weights
    for(int j1=0; j1<n; j1++){
      ZW.col(j1) = Zt.col(j1)*w(j1);        
    } // scaled matrix t(Z)%*%diag(W)
    Z1 = ZW*Z + lambda*H;
    theta_new = solve(Z1, ZW*Y); // need to solve the numerical inverse issue
    resids1 = Y - Z*theta_new;
    check = max(abs(resids1-residsin)/sc); 
    if(check < toler){
      istop = 1;
    }
    // Rcpp::Rcout << "It. " << ic << " Check: " << check << " Theta : \n" << theta_new << "\n";
    residsin = resids1;
  }
  if(istop==0) Rcpp::warning("Estimator did not converge.");
  hat_mat = Z*solve(Z1, ZW);
  fitted = hat_mat*Y;
  hat_diag = hat_mat.diag();
  return Rcpp::List::create(Rcpp::Named("theta_hat") = theta_new,
                            Rcpp::Named("converged") = istop,
                            Rcpp::Named("ic") = ic,
                            Rcpp::Named("resids") = residsin,
                            Rcpp::Named("hat_values") = hat_diag,
                            Rcpp::Named("last_check") = check,
                            Rcpp::Named("weights") = 2*w,
                            Rcpp::Named("fitted") = fitted);
}

// [[Rcpp::export()]]
Rcpp::List IRLSCmult (arma::mat Z, arma::vec Y, arma::vec lambda, arma::mat H,
  int type, arma::vec W, double sc, double tuning, double toler, int imax){

  int n = Y.n_elem, p1 = Z.n_cols, nlambda = lambda.n_elem;  
  // Rcpp::Rcout << nlambda << std::endl;
  arma::vec istop(nlambda, arma::fill::zeros), ic(nlambda, arma::fill::zeros);
  arma::vec w(n), resids1(n);
  arma::mat Zt = trans(Z), ZW(p1,n), Z1(p1,p1); // Zt is of dimension p1-times-n
  arma::mat hat_mat(n,n); // hat matrix for hat values
  arma::mat resids(n,nlambda, arma::fill::ones), hat_diag(n,nlambda, arma::fill::zeros), fitted(n,nlambda, arma::fill::zeros);
  arma::mat theta_new(p1,nlambda, arma::fill::zeros); 
  double check = 0;
  
  for(int k=0; k<nlambda; k++){
    while(istop(k) == 0 && ic(k) < imax){
      ic(k)++;
      w = W % psiwC(resids.col(k)/sc, type, tuning); // element-wise product
      w = w/(sc*sc); // scaling affects only the weights
      for(int j1=0; j1<n; j1++){
        ZW.col(j1) = Zt.col(j1)*w(j1);        
      } // scaled matrix t(Z)%*%diag(W)
      Z1 = ZW*Z + lambda(k)*H;
      theta_new.col(k) = solve(Z1, ZW*Y, arma::solve_opts::fast); // need to solve the numerical inverse issue
      resids1 = Y - Z*theta_new.col(k);
      check = max(abs(resids1-resids.col(k)))/sc; 
      if(check < toler){
        istop(k) = 1;
      }
      resids.col(k) = resids1;
    }
    // if(istop(k)==0) Rcpp::warning("Estimator did not converge.");
    hat_mat = Z*solve(Z1, ZW, arma::solve_opts::fast);
    fitted.col(k) = hat_mat*Y;
    hat_diag.col(k) = hat_mat.diag();
  }
  return Rcpp::List::create(Rcpp::Named("theta_hat") = theta_new,
                            Rcpp::Named("converged") = istop,
                            Rcpp::Named("ic") = ic,
                            Rcpp::Named("resids") = resids,
                            Rcpp::Named("hat_values") = hat_diag,
                            Rcpp::Named("fitted") = fitted);
}

// [[Rcpp::export()]]
Rcpp::List ridgeC (arma::mat Z, arma::vec Y, double lambda, arma::mat H, arma::vec W){
  
  int n = Y.n_elem, p1 = Z.n_cols;
  arma::mat Zt = trans(Z); // Zt is of dimension p1-times-n
  arma::mat hat_mat(n,n); // hat matrix for hat values
  arma::vec theta_new(p1), residsin(n), hat_diag(n), fitted(n);
  
  theta_new = solve(Zt*arma::diagmat(W)*Z + lambda*H, Zt*arma::diagmat(W)*Y);
  fitted = Z*theta_new;
  residsin = Y - fitted;
  hat_mat = Z*solve(Zt*arma::diagmat(W)*Z + lambda*H, Zt*arma::diagmat(W));
  hat_diag = hat_mat.diag();
  return Rcpp::List::create(Rcpp::Named("theta_hat") = theta_new,
                            Rcpp::Named("resids") = residsin,
                            Rcpp::Named("hat_values") = hat_diag,
                            Rcpp::Named("fitted") = fitted);
}


// [[Rcpp::export()]]
Rcpp::List HuberQpC(
    const arma::mat Z,
    const arma::vec Y,
    const arma::mat H,
    double delta = 1.345
) {
    int n = Z.n_rows;
    int p = Z.n_cols;

    int n_vars = p + 2 * n;   // beta (p), a (n), t (n)
    int n_cons = 3 * n;       // 2n + n

    // -------------------------------------------------------------------------
    // Build dense P (block diagonal)
    // -------------------------------------------------------------------------
    arma::mat Hp = H + 1e-6 * arma::eye(p, p);
    arma::mat In1 = (1 + 1e-6) * arma::eye(n, n);
    arma::mat In2 = 1e-6 * arma::eye(n, n);

    arma::mat P = arma::zeros(n_vars, n_vars);
    P.submat(0, 0, p - 1, p - 1) = Hp;
    P.submat(p, p, p + n - 1, p + n - 1) = In1;
    P.submat(p + n, p + n, n_vars - 1, n_vars - 1) = In2;

    // -------------------------------------------------------------------------
    // q vector
    // -------------------------------------------------------------------------
    arma::vec qvec = arma::zeros(n_vars);
    qvec.subvec(p + n, n_vars - 1).fill(delta);

    // -------------------------------------------------------------------------
    // Build A matrix (dense then to sparse)
    // -------------------------------------------------------------------------
    arma::mat Amat = arma::zeros(n_cons, n_vars);
    arma::vec bvec = arma::zeros(n_cons);

    for (int i = 0; i < n; ++i) {
        int row1 = 2 * i;
        int row2 = 2 * i + 1;
        int row3 = 2 * n + i;

        Amat.submat(row1, 0, row1, p - 1) = Z.row(i);
        Amat(row1, p + i) = 1.0;
        Amat(row1, p + n + i) = 1.0;
        bvec(row1) = Y(i);

        Amat.submat(row2, 0, row2, p - 1) = -Z.row(i);
        Amat(row2, p + i) = -1.0;
        Amat(row2, p + n + i) = 1.0;
        bvec(row2) = -Y(i);

        Amat(row3, p + n + i) = 1.0;
        bvec(row3) = 0.0;
    }

    arma::vec lvec = bvec;
    arma::vec uvec(n_cons);
    uvec.fill(std::numeric_limits<double>::infinity());

    // -------------------------------------------------------------------------
    // Convert P and A to CSC arrays suitable for OSQP (new API)
    // We'll build column-by-column (CSC) arrays.
    // For P: OSQP requires the upper-triangular part only. We'll iterate j >= i.
    // -------------------------------------------------------------------------

    // ---- P: extract upper triangular CSC representation ----
    std::vector<OSQPFloat> Px;
    std::vector<OSQPInt>   Pi;
    std::vector<OSQPInt>   Pp; Pp.reserve(n_vars + 1);
    Pp.push_back(0);

    for (int j = 0; j < n_vars; ++j) {
        for (int i = 0; i <= j; ++i) {              // upper triangular (i <= j)
            double v = P(i, j);
            if (v != 0.0) {
                Px.push_back(static_cast<OSQPFloat>(v));
                Pi.push_back(static_cast<OSQPInt>(i));
            }
        }
        Pp.push_back(static_cast<OSQPInt>(Px.size()));
    }
    OSQPInt P_nnz = static_cast<OSQPInt>(Px.size());

    // ---- A: full matrix CSC representation ----
    std::vector<OSQPFloat> Ax;
    std::vector<OSQPInt>   Ai;
    std::vector<OSQPInt>   Ap; Ap.reserve(n_vars + 1);
    Ap.push_back(0);

    for (int j = 0; j < n_vars; ++j) {
        for (int i = 0; i < n_cons; ++i) {
            double v = Amat(i, j);
            if (v != 0.0) {
                Ax.push_back(static_cast<OSQPFloat>(v));
                Ai.push_back(static_cast<OSQPInt>(i));
            }
        }
        Ap.push_back(static_cast<OSQPInt>(Ax.size()));
    }
    OSQPInt A_nnz = static_cast<OSQPInt>(Ax.size());

    // ---- q, l, u as OSQPFloat arrays ----
    std::vector<OSQPFloat> q_osqp(n_vars);
    for (int i = 0; i < n_vars; ++i) q_osqp[i] = static_cast<OSQPFloat>(qvec(i));

    std::vector<OSQPFloat> l_osqp(n_cons), u_osqp(n_cons);
    for (int i = 0; i < n_cons; ++i) {
        l_osqp[i] = static_cast<OSQPFloat>(lvec(i));
        double uval = uvec(i);
        if (std::isinf(uval) && uval > 0) {
            u_osqp[i] = OSQP_INFTY;
        } else if (std::isinf(uval) && uval < 0) {
            u_osqp[i] = -OSQP_INFTY;
        } else {
            u_osqp[i] = static_cast<OSQPFloat>(uval);
        }
    }

    // -------------------------------------------------------------------------
    // Create OSQP CSC matrices (they keep pointers to our arrays)
    // Use OSQPCscMatrix_new (new API)
    // -------------------------------------------------------------------------
    OSQPCscMatrix* P_csc = OSQPCscMatrix_new(static_cast<OSQPInt>(n_vars),
                                             static_cast<OSQPInt>(n_vars),
                                             P_nnz,
                                             Px.data(),
                                             Pi.data(),
                                             Pp.data());

    OSQPCscMatrix* A_csc = OSQPCscMatrix_new(static_cast<OSQPInt>(n_cons),
                                             static_cast<OSQPInt>(n_vars),
                                             A_nnz,
                                             Ax.data(),
                                             Ai.data(),
                                             Ap.data());

    // -------------------------------------------------------------------------
    // Settings and setup
    // -------------------------------------------------------------------------
    OSQPSettings* settings = OSQPSettings_new();
    osqp_set_default_settings(settings);
    settings->alpha = 1.6;
    settings->verbose = 0;
    settings->eps_abs = 1e-5;
    settings->eps_rel = 1e-5;

    // Solver
    OSQPSolver* solver = nullptr;
    OSQPInt exitflag = osqp_setup(&solver,
                                  P_csc,
                                  q_osqp.data(),
                                  A_csc,
                                  l_osqp.data(),
                                  u_osqp.data(),
                                  static_cast<OSQPInt>(n_cons),
                                  static_cast<OSQPInt>(n_vars),
                                  settings);
    if (exitflag) {
        // free created objects
        if (solver) osqp_cleanup(solver);
        OSQPCscMatrix_free(A_csc);
        OSQPCscMatrix_free(P_csc);
        OSQPSettings_free(settings);
        Rcpp::stop("OSQP setup failed (exitflag != 0).");
    }

    // -------------------------------------------------------------------------
    // Solve
    // -------------------------------------------------------------------------
    exitflag = osqp_solve(solver);
    if (exitflag) {
        // cleanup and error
        osqp_cleanup(solver);
        OSQPCscMatrix_free(A_csc);
        OSQPCscMatrix_free(P_csc);
        OSQPSettings_free(settings);
        Rcpp::stop("OSQP solve failed (exitflag != 0).");
    }

    // -------------------------------------------------------------------------
    // Extract solution: solver->solution->x has length n_vars
    // -------------------------------------------------------------------------
    arma::vec theta_new = arma::zeros(p);
    if (solver && solver->solution && solver->info) {
        // status is a char[32]
        if (std::strncmp(solver->info->status, "solved", 6) == 0) {
            OSQPFloat* solx = solver->solution->x; // pointer to primal solution
            if (solx != nullptr) {
                for (int i = 0; i < p; ++i) theta_new(i) = static_cast<double>(solx[i]);
            }
        } else {
            Rcpp::Rcerr << "OSQP did not solve the problem. Status: " << solver->info->status << std::endl;
            // still try to extract x if available
            if (solver->solution && solver->solution->x) {
                for (int i = 0; i < p; ++i) theta_new(i) = static_cast<double>(solver->solution->x[i]);
            }
        }
    } else {
        Rcpp::Rcerr << "OSQP solver or solution/info pointer is null." << std::endl;
    }

    // -------------------------------------------------------------------------
    // Cleanup: free OSQP structures (this does not free our std::vectors memory,
    // which will be released automatically at scope exit)
    // -------------------------------------------------------------------------
    osqp_cleanup(solver);
    OSQPCscMatrix_free(A_csc);
    OSQPCscMatrix_free(P_csc);
    OSQPSettings_free(settings);

    return Rcpp::List::create(Rcpp::Named("theta_hat") = theta_new);
}
