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
double BICw1(double n, const arma::mat& R1, const arma::colvec& d,
           const arma::colvec& w1, int BICtype){
  // d = R12*w2. w1 should not be normalized.
  double myBIC = 0.0;
  double rss = as_scalar(trans(w1)*R1*w1 - 2*trans(w1)*d + 1); // trans(w2)*R2*w2 should be 1, so it added as + 1 at the end.
  if (BICtype == 1){ // BIC1
    myBIC = n*rss + sum(w1!=0)*log(n);
  } else if (BICtype == 2){ // BIC2
    myBIC = n*log(rss*n/(n-sum(w1!=0))) + sum(w1!=0)*log(n);
  } else if (BICtype == 3){ // BIC in Wilms and Croux Rcode
    myBIC = 2*rss + sum(w1!=0)*log(n);
  } else if (BICtype == 4){ // BIC in Wilms and Croux paper
    myBIC = n*log(rss) + sum(w1!=0)*log(n);
  } else {
    warning("Choose BIC type either 1 or 2.\n");
  }
  return myBIC;
}


// find lasso solution for w1 given w2.
// [[Rcpp::export]]
Rcpp::List lassobic(double n, const arma::mat& R1, const arma::colvec& d, //d = R12*w2
                    arma::colvec w1init, const arma::colvec& lamseq,
                    int BICtype, int maxiter = 1000, double tol = 0.0001, int lassoverbose = 1){
  // basically same as solveLasso for fixed w2 and all lambda values.
  // find w1 with smallest bic

  // normalize the initial value vectors
  w1init = w1init/sqrt(as_scalar(trans(w1init)*R1*w1init));

  // initialize variables
  int p = d.n_elem;
  int nlam = lamseq.n_elem;
  int ind = 0;
  arma::colvec bicvec(nlam);
  arma::colvec w1(p); w1.fill(0); // This w1 is for final output. Not a starting value. Starting value is w1init.
  // And in case minimum lambda value is larger than max(abs(R12*w2)), this function result in zero solution.
  arma::mat wmat(p, nlam); //Needed to initialize to track all solutions (not normalized). BUT in case the minimal lambda is below lambda max, it will return zero vector in line 94.

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
          bicvec[j] = BICw1(n, R1, d, w1init, BICtype);
            // print error if iter reach maxiter
            if (iter >= maxiter && lassoverbose == 1){
              warning("Failed to converge: lasso part at %i-th lambda = %f with error = %f\n", j, lamseq[j], error);
            }
        }// finished checking all lambda values.

    ind = bicvec.index_min(); // Based on bic, choose the tuning parameter with the smallest bic value.
    w1 = wmat.col(ind); // This is the lasso solution. Unless it has all zeros, it will be normalized below.
    if (as_scalar(trans(w1)*R1*w1) == 0){
      w1.fill(0);
    } else {
      w1 = w1/sqrt(as_scalar(trans(w1)*R1*w1));
    }

  } else { //in case the minimal lambda is below lambda max,
    //ind will be a zero scalar and w1 is a zero vector of length p.
    wmat.zeros(p, 1); // it will return zero for wmat. Originally initiated as p by nlam matrix, but here become p by 1 matrix.
    bicvec.fill(0); // a zero vector of length nlam.
  }
  return List::create(Rcpp::Named("finalcoef") = w1,
                      Rcpp::Named("wmat") = wmat, // NOTE THAT wmat is "NOT" normalized lasso solution.
                      Rcpp::Named("bicInd") = ind+1, // vector index in cpp starts from 0, so had to revise "bicInd" variable by adding 1.
                      Rcpp::Named("lamseq") = lamseq,
                      Rcpp::Named("bicvalues") = bicvec);
}


