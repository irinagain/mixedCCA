//[[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include "stdio.h"
using namespace Rcpp;
using namespace arma;

double sign(double x){ // sign function
  return (x > 0) - (0 > x);
}

double soft(double a, double lambda){ // Soft-thresholding function
  return sign(a)*std::max(0.0, std::abs(a)-lambda);
}

// Using BIC criterion,
double BIC(int n, const arma::mat& R1, const arma::mat& R2, const arma::mat& R12, const arma::colvec& w1, const arma::colvec& w2, int BICtype){
  double BIC = datum::nan;
  double rss = as_scalar(trans(w1)*R1*w1 - 2*trans(w1)*R12*w2 + trans(w2)*R2*w2);

  if (BICtype == 1){
    BIC = rss + sum(w1!=0)*log(n)/n;
  } else if (BICtype == 2){
    BIC = log(rss*n/(n-sum(w1!=0))) + sum(w1!=0)*log(n)/n;
  } else { printf("Choose BIC type either 1 or 2.\n"); }

  return BIC;
}


arma::colvec lassobic(int n, const arma::mat& R1, const arma::mat& R2, const arma::mat& R12, arma::colvec winit1, const arma::colvec& w2, const arma::colvec& lambda,
                      int BICtype,
                      int maxiter = 100, double tol = 0.001, bool verbose = true){
  // basically same as solveLasso for fixed w2 and all lambda values.
  // find w1 with smallest bic
  arma::colvec d = R12*w2;
  int p = d.n_elem;
  int nlam = lambda.n_elem;
  arma::colvec bicvec(nlam);

  if (arma::min(lambda) >= arma::max(arma::abs(d))) {
    std::fill(winit1.begin(), winit1.end(), 0.0);
    return winit1;
  } else {

    arma::colvec r = d - R1*winit1;
    arma::mat wmat(p, nlam);

    for (int j=0; j<nlam; j++){
      double error = 1000.0;
      int iter = 0;

      while( iter <= maxiter && error > tol ){
        iter += 1;
        error = 0.0;
        for (int i=0; i<p; i++){
          double wold = winit1[i];
          winit1[i] = soft(r[i] + wold, lambda[j]);
          r += R1.col(i)*(wold - winit1[i]);
          error += std::abs(winit1[i] - wold);
        }
      }
      wmat.col(j) = winit1;
      if (iter == (maxiter + 1) && verbose){
        printf("Failed to converge. Try increasing the number of iterations. (lasso part)\n error = %f\n", error);
      }
      bicvec[j] = BIC(n, R1, R2, R12, winit1, w2, BICtype);
    }
    int ind = bicvec.index_min();
    return wmat.col(ind);

  }
}


// [[Rcpp::export]]
arma::colvec find_w12bic(int n, const arma::mat& R1, const arma::mat& R2, const arma::mat& R12,
                         const arma::vec& lambda1, const arma::vec& lambda2,
                         arma::colvec w1init, arma::colvec w2init,
                         int BICtype,
                         int maxiter = 1000, double tol = 0.001, bool verbose = true){
  // Same as find_w12 function. SolveLasso part is replaced with lassobic.
  int p1 = R1.n_cols;
  int p2 = R2.n_cols;
  w1init = w1init/sqrt(as_scalar(trans(w1init)*R1*w1init));
  w2init = w2init/sqrt(as_scalar(trans(w2init)*R2*w2init));

  double error = 1000.0;
  int iter = 0;
  arma::colvec w1(p1);
  arma::colvec w2(p2);

  while( iter <= maxiter && error > tol ){
    iter += 1;
    error = 0.0;

    w1 = lassobic(n, R1, R2, R12, w1init, w2init, lambda1, BICtype);
    if (sum(arma::abs(w1))==0){
      std::fill(w2.begin(), w2.end(), 0.0);
      w1.insert_rows(w1.n_elem, w2);
      return w1;
    }
    w1 = w1/sqrt(as_scalar(trans(w1)*R1*w1));
    error += sum(arma::abs(w1-w1init));
    w1init = w1;

    w2 = lassobic(n, R2, R1, trans(R12), w2init, w1init, lambda2, BICtype);
    if (sum(arma::abs(w2))==0){
      std::fill(w1.begin(), w1.end(), 0.0);
      w1.insert_rows(w1.n_elem, w2);
      return w1;
    }
    w2 = w2/sqrt(as_scalar(trans(w2)*R2*w2));
    error += sum(arma::abs(w2-w2init));
    w2init = w2;
  }
  if (iter == (maxiter + 1) && verbose){
    printf("Failed to converge. Try increasing the number of iterations. (findw12 part)\n error = %f\n", error);
  }
  w1init.insert_rows(w1init.n_elem, w2init);
  return w1init;
}
