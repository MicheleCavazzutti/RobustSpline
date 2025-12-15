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
arma::vec psiwC (arma::vec t, int type, double alpha, double tuning){
// function that gives the weights psi(t)/(2*t) (with 2 in the denominator)
// for reweighting in the IRLS algorithm

// type: the type of the loss function used. Encoding:
//       1 for the absolute (quantile) loss
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
        if(t(i)>0) res(i) = alpha/(2*t(i)); // 0.5/t(i);
        if(t(i)<0) res(i) = (alpha-1)/(2*t(i)); // -0.5/t(i);
        if(t(i)==0) res(i) = R_PosInf; 
      }
    } else {
      for(int i=0; i<n; i++){
        if(abs(t(i))<=tuning){
          res(i) = alpha/(2*tuning); // 0.5/tuning;
        } else {
          if(t(i)>0) res(i) = alpha/(2*t(i)); // 0.5/t(i);
          if(t(i)<0) res(i) = (alpha-1)/(2*t(i)); // -0.5/t(i);
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
  int type, double alpha, arma::vec W, double sc, arma::vec residsin, double tuning, double toler, int imax){
  
  int ic = 0, istop = 0;
  int n = Y.n_elem, p1 = Z.n_cols;
  arma::vec w(n);
  arma::mat Zt = trans(Z), ZW(p1,n), Z1(p1,p1); // Zt is of dimension p1-times-n
  arma::mat hat_mat(n,n); // hat matrix for hat values
  arma::vec theta_new(p1), resids1(n), hat_diag(n), fitted(n);
  double check = 0;
  
  while(istop == 0 && ic < imax){
    ic++;
    w = W % psiwC(residsin/sc, type, alpha, tuning); // element-wise product
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
  int type, double alpha, arma::vec W, double sc, double tuning, double toler, int imax){

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
      w = W % psiwC(resids.col(k)/sc, type, alpha, tuning); // element-wise product
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


using namespace Rcpp;

// ----------------------------------------------------------------------------
// Helper: convert arma::sp_mat -> CSC vectors (column-major, rows sorted)
// ----------------------------------------------------------------------------
static void spmat_to_csc_vectors(
    const arma::sp_mat& S,
    std::vector<OSQPFloat>& x_out,
    std::vector<OSQPInt>& i_out,
    std::vector<OSQPInt>& p_out
) {
    const arma::uword ncol = S.n_cols;
    x_out.clear();
    i_out.clear();
    p_out.clear();
    p_out.reserve(ncol + 1);
    p_out.push_back(0);

    for (arma::uword col = 0; col < ncol; ++col) {
        // iterate over nonzeros in this column (Armadillo gives row-ordered iteration)
        for (arma::sp_mat::const_iterator it = S.begin_col(col); it != S.end_col(col); ++it) {
            double val = static_cast<double>(*it);
            if (val != 0.0) {
                x_out.push_back(static_cast<OSQPFloat>(val));
                i_out.push_back(static_cast<OSQPInt>(it.row()));
            }
        }
        p_out.push_back(static_cast<OSQPInt>(x_out.size()));
    }
}

// ----------------------------------------------------------------------------
// Helper: build block-diagonal sparse P like bdiag(H+epsI, (1+eps)I_n, eps I_n)
// and return lower-triangular part (equivalent to Matrix::bdiag behavior)
// ----------------------------------------------------------------------------
static arma::sp_mat make_block_P_sp(const arma::mat& H, int p, int n, double eps = 1e-6) {
  arma::mat Hp = 0.5 * (H + H.t());
  
  for (int i = 0; i < p; ++i) {
    if (std::fabs(Hp(i,i)) < 1e-12)
      Hp(i,i) = 1e-12;
  }
  
  //----------------------------------------------------------------------
  // 3. Add eps * I (guaranteed > 0)
  //----------------------------------------------------------------------
  Hp = Hp + eps * arma::eye<arma::mat>(p, p);
  
  arma::sp_mat Hs(Hp);   // convert to sparse
  
  //----------------------------------------------------------------------
  // 4. Build (1+eps) I_n and eps I_n with guaranteed positive diagonal
  //----------------------------------------------------------------------
  arma::sp_mat I1 = arma::speye<arma::sp_mat>(n, n);
  arma::sp_mat I2 = arma::speye<arma::sp_mat>(n, n);
  
  for (int i = 0; i < n; ++i) {
    I1(i, i) = 1.0 + eps;
    I2(i, i) = eps > 1e-12 ? eps : 1e-12;
  }
    
  //----------------------------------------------------------------------
  // 5. Assemble block diagonal sparse matrix
  //----------------------------------------------------------------------
  int N = p + n + n;
  arma::sp_mat Psp(N, N);
    
  Psp.submat(0,      0,      p - 1,     p - 1)     = Hs;
  Psp.submat(p,      p,      p + n - 1, p + n - 1) = I1;
  Psp.submat(p + n,  p + n,  p + 2*n - 1, p + 2*n - 1) = I2;

  //----------------------------------------------------------------------
  // 6. Return lower triangular part (OSQP expects triangular)
  //----------------------------------------------------------------------
  //return arma::trimatl(Psp);
  return arma::trimatu(Psp);
}

// ----------------------------------------------------------------------------
// Helper: convert dense Amat (arma::mat) to sparse arma::sp_mat (removes zeros)
// ----------------------------------------------------------------------------
static arma::sp_mat make_A_sp(const arma::mat& Adense) {
    arma::sp_mat Asp(Adense); // Armadillo will drop exact zeros
    return Asp;
}

// [[Rcpp::export()]]
Rcpp::List HuberQpC(
    const arma::mat Z,
    const arma::vec Y,
    const arma::mat H,
    const arma::vec w,
    double delta = 1.345
) {
    int n = (int) Z.n_rows;
    int p = (int) Z.n_cols;

    int n_vars = p + 2 * n;   // beta (p), a (n), t (n)
    int n_cons = 3 * n;       // 2n + n

    // Apply weight: diag(w) %*% Z  (Z each_col % w)
    arma::mat Zw = Z.each_col() % w;

    // w * Y (element-wise)
    arma::vec Yw = w % Y;

    // -------------------------------------------------------------------------
    // Build P as sparse, following bdiag(H + eps I, (1+eps)I_n, eps I_n)
    // -------------------------------------------------------------------------
    arma::sp_mat Psp = make_block_P_sp(H, p, n, 1e-6);

    // -------------------------------------------------------------------------
    // q vector
    // -------------------------------------------------------------------------
    arma::vec qvec = arma::zeros<arma::vec>(n_vars);
    for (int i = p + n; i < n_vars; ++i) qvec(i) = delta;

    // -------------------------------------------------------------------------
    // Build A matrix (dense then convert to sparse)
    // -------------------------------------------------------------------------
    arma::mat Amat = arma::zeros<arma::mat>(n_cons, n_vars);
    arma::vec bvec = arma::zeros<arma::vec>(n_cons);

    for (int ii = 0; ii < n; ++ii) {
        int row1 = 2 * ii;
        int row2 = 2 * ii + 1;
        int row3 = 2 * n + ii;

        // Zw.row(ii) is 1 x p
        Amat.submat(row1, 0, row1, p - 1) = Zw.row(ii);
        Amat(row1, p + ii) = 1.0;
        Amat(row1, p + n + ii) = 1.0;
        bvec(row1) = Yw(ii);

        Amat.submat(row2, 0, row2, p - 1) = -Zw.row(ii);
        Amat(row2, p + ii) = -1.0;
        Amat(row2, p + n + ii) = 1.0;
        bvec(row2) = -Yw(ii);

        Amat(row3, p + n + ii) = 1.0;
        bvec(row3) = 0.0;
    }

    arma::vec lvec = bvec;
    std::vector<double> uvec(n_cons);
    for (int i = 0; i < n_cons; ++i) uvec[i] = std::numeric_limits<double>::infinity();

    // Convert Amat to sparse (drops zeros, sorts indices)
    arma::sp_mat Asp = make_A_sp(Amat);

    // -------------------------------------------------------------------------
    // Convert sparse arma::sp_mat -> CSC vectors for OSQP (P and A)
    // -------------------------------------------------------------------------
    std::vector<OSQPFloat> Px;
    std::vector<OSQPInt>   Pi;
    std::vector<OSQPInt>   Pp;
    spmat_to_csc_vectors(Psp, Px, Pi, Pp);
    OSQPInt P_nnz = static_cast<OSQPInt>(Px.size());

    std::vector<OSQPFloat> Ax;
    std::vector<OSQPInt>   Ai;
    std::vector<OSQPInt>   Ap;
    spmat_to_csc_vectors(Asp, Ax, Ai, Ap);
    OSQPInt A_nnz = static_cast<OSQPInt>(Ax.size());

    // ---- q, l, u as OSQPFloat arrays ----
    std::vector<OSQPFloat> q_osqp(n_vars);
    for (int i = 0; i < n_vars; ++i) q_osqp[i] = static_cast<OSQPFloat>(qvec(i));

    std::vector<OSQPFloat> l_osqp(n_cons), u_osqp(n_cons);
    for (int i = 0; i < n_cons; ++i) {
        l_osqp[i] = static_cast<OSQPFloat>(lvec(i));
        double uval = uvec[i];
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
    //settings->alpha = 1.6;
    //settings->verbose = 0;
    settings->eps_abs = 1e-1;
    settings->eps_rel = 1e-1;
    settings->max_iter = 10000;

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
        if (solver) osqp_cleanup(solver);
        OSQPCscMatrix_free(A_csc);
        OSQPCscMatrix_free(P_csc);
        OSQPSettings_free(settings);
	Rcpp::Rcerr << "OSQP error code: " << exitflag << "\n";
	if (solver && solver->info) {
	  Rcpp::Rcerr << "OSQP info: " << solver->info->status << "\n";
	}
        Rcpp::stop("OSQP setup failed (exitflag != 0).");
    }

    // -------------------------------------------------------------------------
    // Solve
    // -------------------------------------------------------------------------
    exitflag = osqp_solve(solver);
    if (exitflag) {
        osqp_cleanup(solver);
        OSQPCscMatrix_free(A_csc);
        OSQPCscMatrix_free(P_csc);
        OSQPSettings_free(settings);
        Rcpp::stop("OSQP solve failed (exitflag != 0).");
    }

    // -------------------------------------------------------------------------
    // Extract solution: solver->solution->x has length n_vars
    // -------------------------------------------------------------------------
    arma::vec theta_new = arma::zeros<arma::vec>(p);
    if (solver && solver->solution && solver->info) {
        if (std::strncmp(solver->info->status, "solved", 6) == 0) {
            OSQPFloat* solx = solver->solution->x;
            if (solx != nullptr) {
                for (int i = 0; i < p; ++i) theta_new(i) = static_cast<double>(solx[i]);
            }
        } else {
            Rcpp::Rcerr << "OSQP did not solve the problem. Status: " << solver->info->status << std::endl;
            if (solver->solution && solver->solution->x) {
                for (int i = 0; i < p; ++i) theta_new(i) = static_cast<double>(solver->solution->x[i]);
            }
        }
    } else {
        Rcpp::Rcerr << "OSQP solver or solution/info pointer is null." << std::endl;
    }

    // -------------------------------------------------------------------------
    // Cleanup: free OSQP structures (this does not free our std::vectors memory)
    // -------------------------------------------------------------------------
    osqp_cleanup(solver);
    OSQPCscMatrix_free(A_csc);
    OSQPCscMatrix_free(P_csc);
    OSQPSettings_free(settings);

    // Prepare the output (same as prima)
    arma::mat ZtW = Z.t() * (Z.each_col() % w);   // p x p

    arma::mat ZtW2 = Z.t();
    ZtW2.each_row() %= w.t();                     // p x n

    // Matrix to invert: Xt W X + lambda H  (qui H è già passato correttamente)
    arma::mat Ainv = ZtW + H;

    arma::mat hat = Z * arma::solve(Ainv, ZtW2, arma::solve_opts::likely_sympd);

    arma::vec fitted = Z * theta_new;
    arma::vec resid = Y - fitted;

    arma::vec hat_diag = hat.diag();

    return Rcpp::List::create(Rcpp::Named("theta_hat") = theta_new,
                              Rcpp::Named("resids") = resid,
                              Rcpp::Named("hat_values") = hat_diag,
                              Rcpp::Named("fitted") = fitted);
}
