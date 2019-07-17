#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// sign function
double sign(double x){
  return (x > 0) - (0 > x);
}

// Soft-thresholding function
double soft(double a, double lambda){
  return sign(a)*std::max(0.0, std::abs(a)-lambda);
}

// BIC criterion for w1, given w2.
double myBIC(int n, const arma::mat& R1, const arma::colvec& d,
           const arma::colvec& w1, int BICtype){
  // d = R12*w2. w1 should not be normalized.
  double myBIC = 0.0;
  double rss = as_scalar(trans(w1)*R1*w1 - 2*trans(w1)*d + 1);
  if (BICtype == 1){ // BIC1
    myBIC = n*rss + sum(w1!=0)*log(n);
  } else if (BICtype == 2){ // BIC2
    myBIC = n*log(rss*n/(n-sum(w1!=0))) + sum(w1!=0)*log(n);
  } else if (BICtype == 3){ // BIC as in Wilms and Croux paper
    myBIC = 2*rss + sum(w1!=0)*log(n);
  } else {
    warning("Choose BIC type either 1 or 2.\n");
  }
  return myBIC;
}


// find lasso solution for w1 given w2.
// [[Rcpp::export]]
Rcpp::List lassobic(int n, const arma::mat& R1, const arma::mat& R2, const arma::mat& R12,
                    arma::colvec w1init, arma::colvec w2, const arma::colvec& lamseq,
                    int BICtype, int maxiter = 1000, double tol = 0.01,
                    bool convcheck = true){
  // basically same as solveLasso for fixed w2 and all lambda values.
  // find w1 with smallest bic

  // normalize the initial value vectors
  w1init = w1init/sqrt(as_scalar(trans(w1init)*R1*w1init));
  w2 = w2/sqrt(as_scalar(trans(w2)*R2*w2));

  // initialize variables
  arma::colvec d = R12*w2;
  int p = d.n_elem;
  int nlam = lamseq.n_elem;
  int ind = 0;
  arma::colvec bicvec(nlam);
  arma::colvec w1(p);
  w1.fill(0);
  arma::mat wmat(p, nlam);

  if (arma::min(lamseq) < arma::max(arma::abs(d))) {
    arma::colvec r = d - R1*w1init;
    // for each lambda value, estimate w1
    for (int j=0; j<nlam; j++){
      double error = 1000.0;
      int iter = 0;
      // optimization
      while( iter <= maxiter && error > tol ){
        iter += 1;
        error = 0.0;
        for (int i=0; i<p; i++){
          double wold = w1init[i];
          w1init[i] = soft(r[i] + wold, lamseq[j]);
          r += R1.col(i)*(wold - w1init[i]);
          error += std::abs(w1init[i] - wold);
        }
      }// done with lasso optimization for one lambda value.
      wmat.col(j) = w1init; // save
      bicvec[j] = myBIC(n, R1, d, w1init, BICtype);
      // print error if iter reached maxiter
      if (iter == (maxiter + 1) && convcheck){
        warning("Failed to converge: lasso part at %i-th lambda = %f with error = %f\n", j, lamseq[j], error);
      }
    }// finished checking all lambda values.
    ind = bicvec.index_min();
    w1 = wmat.col(ind);
    if (as_scalar(trans(w1)*R1*w1) == 0){
      w1.fill(0);
    } else {
      w1 = w1/sqrt(as_scalar(trans(w1)*R1*w1));
    }

  }
  return List::create(Rcpp::Named("finalcoef") = w1,
                      Rcpp::Named("wmat") = wmat,
                      Rcpp::Named("bicInd") = ind+1, // vector index in cpp starts from 0, so had to revise "bicInd" variable by adding 1.
                      Rcpp::Named("lamseq") = lamseq,
                      Rcpp::Named("bicvalues") = bicvec);
}


