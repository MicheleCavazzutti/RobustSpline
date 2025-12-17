// we include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

#include <Rcpp.h>
using namespace Rcpp;

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat distn(arma::mat X, int n, int d)
{
    // matrix of distances of all rows of X from each other
    arma::mat res(n,n, arma::fill::zeros);
    for(int i = 0; i < n; i++){
        for(int j = i+1; j < n; j++){
            res(i,j) = res(j,i) = sqrt(sum(pow(X.row(i) - X.row(j),2)));
        }
    }
    return(res);
}

// [[Rcpp::export]]
arma::mat distnAB(arma::mat A, arma::mat B, int m, int n, int d)
{
    // matrix of distances of all rows of A from all rows of B
    arma::mat res(m,n, arma::fill::zeros);
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            res(i,j) = sqrt(sum(pow(A.row(i) - B.row(j),2)));
        }
    }
    return(res);
}

// [[Rcpp::export]]
arma::mat solveC(arma::mat A, arma::mat B)
{
    // solve function implemented in armadillo
    return(solve(A, B));
}
