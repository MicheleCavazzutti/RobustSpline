# define ARMA_WARN_LEVEL 1  // Turns off warnings about inverses of nearly singular matrices
# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//#include <osqp/osqp.h>
#include "./osqp/include/public/osqp.h"
#include <vector>
#include <limits>
#include <cstring>
#include <algorithm>
//#include "osqp.h"


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
// Helper: Fast CSC extraction
// Directly maps Armadillo sparse structure to OSQP vectors to save time/memory
// ----------------------------------------------------------------------------
static void spmat_to_csc_pointers(
    const arma::sp_mat& S,
    std::vector<OSQPFloat>& x_out,
    std::vector<OSQPInt>& i_out,
    std::vector<OSQPInt>& p_out
) {
    x_out.assign(S.values, S.values + S.n_nonzero);
    i_out.resize(S.n_nonzero);
    for(arma::uword k=0; k < S.n_nonzero; ++k) i_out[k] = static_cast<OSQPInt>(S.row_indices[k]);
    p_out.resize(S.n_cols + 1);
    for(arma::uword k=0; k <= S.n_cols; ++k) p_out[k] = static_cast<OSQPInt>(S.col_ptrs[k]);
}

// ----------------------------------------------------------------------------
// Helper: Build Sparse P for Huber
// P = bdiag(H + eps*I, (1+eps)*I_n, eps*I_n)
// ----------------------------------------------------------------------------
static arma::sp_mat make_P_huber(const arma::mat& H, int p, int n, double eps) {
    arma::sp_mat Hs = arma::sp_mat(0.5 * (H + H.t()));
    Hs.diag() += eps;

    arma::sp_mat I1 = arma::speye<arma::sp_mat>(n, n) * (1.0 + eps);
    arma::sp_mat I2 = arma::speye<arma::sp_mat>(n, n) * std::max(eps, 1e-12);

    // Efficient sparse block assembly
    arma::sp_mat P = arma::join_cols(
        arma::join_rows(Hs, arma::sp_mat(p, 2 * n)),
        arma::join_cols(
            arma::join_rows(arma::sp_mat(n, p), arma::join_rows(I1, arma::sp_mat(n, n))),
            arma::join_rows(arma::sp_mat(n, p + n), I2)
        )
    );
    return arma::trimatu(P);
}

// ----------------------------------------------------------------------------
// Helper: Build Sparse P for Quantile
// P = bdiag(2*lambda*H + eps*I, eps*I_n, eps*I_n)
// ----------------------------------------------------------------------------
static arma::sp_mat make_P_quantile(const arma::mat& H, double lambda, int p, int n, double eps) {
    arma::sp_mat Hs = arma::sp_mat(lambda * (H + H.t())); // 2 * 0.5 * lambda
    Hs.diag() += eps;

    arma::sp_mat In_eps = arma::speye<arma::sp_mat>(n, n) * std::max(eps, 1e-12);

    arma::sp_mat P = arma::join_cols(
        arma::join_rows(Hs, arma::sp_mat(p, 2 * n)),
        arma::join_cols(
            arma::join_rows(arma::sp_mat(n, p), arma::join_rows(In_eps, arma::sp_mat(n, n))),
            arma::join_rows(arma::sp_mat(n, p + n), In_eps)
        )
    );
    return arma::trimatu(P);
}

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// ----------------------------------------------------------------------------
// Helper: Fast CSC extraction per OSQP
// ----------------------------------------------------------------------------
static void spmat_to_osqp_csc(
    const arma::sp_mat& S,
    std::vector<OSQPFloat>& x,
    std::vector<OSQPInt>& i,
    std::vector<OSQPInt>& p
) {
    x.assign(S.values, S.values + S.n_nonzero);
    i.resize(S.n_nonzero);
    for(arma::uword k = 0; k < S.n_nonzero; ++k) {
        i[k] = static_cast<OSQPInt>(S.row_indices[k]);
    }
    p.resize(S.n_cols + 1);
    for(arma::uword k = 0; k <= S.n_cols; ++k) {
        p[k] = static_cast<OSQPInt>(S.col_ptrs[k]);
    }
}

// [[Rcpp::export]]
Rcpp::List HuberQpC(const arma::sp_mat& Z, 
                    const arma::vec& Y, 
                    const arma::mat& H, 
                    const arma::vec& w, 
                    double delta = 1.345) {
    int n = static_cast<int>(Z.n_rows);
    int p = static_cast<int>(Z.n_cols);
    int n_vars = p + 2 * n;
    int n_cons = 3 * n;

    // 1. Matrice P
    arma::sp_mat Psp(n_vars, n_vars);
    Psp.submat(0, 0, p-1, p-1) = arma::sp_mat(H); // H regolarizzazione (già 2*lambda*H)
    for(int i = 0; i < n; ++i) {
      Psp(p + i, p + i) = w(i) + 1e-12;        // Pesi sulla parte quadratica 'a'
      Psp(p + n + i, p + n + i) = 1e-12;       // Regolarizzazione per 't'
    }
    
    Psp = arma::trimatu(Psp);

    // 2. Vettore q
    arma::vec qvec = arma::zeros<arma::vec>(n_vars);
    for(int i = 0; i < n; ++i) {
      qvec(p + n + i) = w(i) * delta;          // Peso * delta solo su t
    }

    // 3. Matrice A
    // FORZIAMO la valutazione dell'espressione per evitare l'errore di conversione
    arma::sp_mat sp_In = arma::speye<arma::sp_mat>(n, n);
    arma::sp_mat Asp(n_cons, n_vars);
    Asp.submat(0, 0, n-1, p-1) = Z;
    Asp.submat(0, p, n-1, p+n-1) = sp_In;
    Asp.submat(0, p+n, n-1, n_vars-1) = sp_In;
    Asp.submat(n, 0, 2*n-1, p-1) = -Z;
    Asp.submat(n, p, 2*n-1, p+n-1) = -sp_In;
    Asp.submat(n, p+n, 2*n-1, n_vars-1) = sp_In;
    Asp.submat(2*n, p+n, 3*n-1, n_vars-1) = sp_In;

    arma::vec lvec(n_cons);
    lvec.subvec(0, n-1)        = Y;
    lvec.subvec(n, 2*n-1)      = -Y;
    lvec.subvec(2*n, 3*n-1).zeros();
    
    arma::vec uvec(n_cons);
    uvec.fill(OSQP_INFTY);

    // 4. OSQP Setup
    std::vector<OSQPFloat> Px, Ax;
    std::vector<OSQPInt> Pi, Pp, Ai, Ap;
    spmat_to_osqp_csc(Psp, Px, Pi, Pp);
    spmat_to_osqp_csc(Asp, Ax, Ai, Ap);

    OSQPCscMatrix* P_osqp = OSQPCscMatrix_new(n_vars, n_vars, (OSQPInt)Px.size(), Px.data(), Pi.data(), Pp.data());
    OSQPCscMatrix* A_osqp = OSQPCscMatrix_new(n_cons, n_vars, (OSQPInt)Ax.size(), Ax.data(), Ai.data(), Ap.data());
    OSQPSettings* settings = (OSQPSettings*)malloc(sizeof(OSQPSettings));
    if(settings) osqp_set_default_settings(settings);
    settings->verbose = 0;
    settings->eps_abs = 1e-6; 
    settings->eps_rel = 1e-6;

    OSQPSolver* solver = nullptr;
    osqp_setup(&solver, P_osqp, qvec.memptr(), A_osqp, lvec.memptr(), uvec.memptr(), n_cons, n_vars, settings);
    osqp_solve(solver);

    arma::vec theta = arma::zeros<arma::vec>(p);
    if (solver && solver->solution) {
        for (int i = 0; i < p; ++i) theta(i) = solver->solution->x[i];
    }

    osqp_cleanup(solver);
    free(P_osqp); free(A_osqp); free(settings);

    // Variabili intermedie per Rcpp
    arma::vec fitted_vals = Z * theta;
    arma::vec residuals = Y - fitted_vals;

    return Rcpp::List::create(
        Rcpp::Named("theta_hat") = theta,
        Rcpp::Named("fitted") = fitted_vals,
        Rcpp::Named("resids") = residuals
    );
}

// [[Rcpp::export]]
Rcpp::List QuantileQpC(const arma::sp_mat& Z, 
                       const arma::vec& Y, 
                       double lambda, 
                       const arma::mat& H, 
                       const arma::vec& w, 
                       double alpha = 0.5) {
    int n = static_cast<int>(Z.n_rows);
    int p = static_cast<int>(Z.n_cols);
    int n_vars = p + 2 * n;
    int n_cons = 3 * n;

    arma::sp_mat Psp(n_vars, n_vars);
    arma::mat H_p = 2.0 * lambda * H;
    H_p.diag() += 1e-6; 
    Psp.submat(0, 0, p-1, p-1) = arma::sp_mat(H_p);
    for(int j = p; j < n_vars; ++j) Psp(j, j) = 1e-6;

    Psp = arma::trimatu(Psp);

    arma::vec qvec = arma::zeros<arma::vec>(n_vars);
    qvec.subvec(p, p + n - 1) = alpha * w;
    qvec.subvec(p + n, n_vars - 1) = (1.0 - alpha) * w;

    arma::sp_mat sp_In = arma::speye<arma::sp_mat>(n, n);
    arma::sp_mat Asp(n_cons, n_vars);
    Asp.submat(0, 0, n-1, p-1) = Z;
    Asp.submat(0, p, n-1, p+n-1) = sp_In;
    Asp.submat(0, p+n, n-1, n_vars-1) = -sp_In;
    Asp.submat(n, p, 2*n-1, p+n-1) = sp_In;
    Asp.submat(2*n, p+n, 3*n-1, n_vars-1) = sp_In;

    arma::vec lvec = arma::join_cols(Y, arma::zeros<arma::vec>(2 * n));
    arma::vec uvec = lvec;
    uvec.subvec(n, 3 * n - 1).fill(OSQP_INFTY);

    std::vector<OSQPFloat> Px, Ax;
    std::vector<OSQPInt> Pi, Pp, Ai, Ap;
    spmat_to_osqp_csc(Psp, Px, Pi, Pp);
    spmat_to_osqp_csc(Asp, Ax, Ai, Ap);

    OSQPCscMatrix* P_osqp = OSQPCscMatrix_new(n_vars, n_vars, (OSQPInt)Px.size(), Px.data(), Pi.data(), Pp.data());
    OSQPCscMatrix* A_osqp = OSQPCscMatrix_new(n_cons, n_vars, (OSQPInt)Ax.size(), Ax.data(), Ai.data(), Ap.data());
    OSQPSettings* settings = (OSQPSettings*)malloc(sizeof(OSQPSettings));
    if(settings) osqp_set_default_settings(settings);
    settings->verbose = 0;
    settings->eps_abs = 1e-4; 
    settings->eps_rel = 1e-4; 

    OSQPSolver* solver = nullptr;
    osqp_setup(&solver, P_osqp, qvec.memptr(), A_osqp, lvec.memptr(), uvec.memptr(), n_cons, n_vars, settings);
    osqp_solve(solver);

    arma::vec theta = arma::zeros<arma::vec>(p);
    if (solver && solver->solution) {
        for (int i = 0; i < p; ++i) theta(i) = solver->solution->x[i];
    }

    osqp_cleanup(solver);
    free(P_osqp); free(A_osqp); free(settings);

    arma::vec fitted_vals = Z * theta;
    arma::vec residuals = Y - fitted_vals;

    return Rcpp::List::create(
        Rcpp::Named("theta_hat") = theta,
        Rcpp::Named("fitted") = fitted_vals,
        Rcpp::Named("resids") = residuals
    );
}
