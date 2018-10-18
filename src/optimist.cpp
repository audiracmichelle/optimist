#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double fn_contribParity(NumericVector &x, NumericMatrix &Sigma) {
  int N = x.size();
  NumericVector g(N);
  NumericVector R(N);
  double f = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      R(i) += Sigma(i,j) * x(j);
    }
  }
  for (int i = 0; i < N - 1; i++) {
    for (int j = i + 1; j < N; j++) {
      double C = x(i) * R(i) - x(j) * R(j);
      f += 2 * C * C;
    }
  }
  return 1e6 * f;
}
                    
// [[Rcpp::export]]              
NumericVector gr_contribParity(NumericVector &x, NumericMatrix &Sigma) {
  int N = x.size();
  NumericVector g(N);
  NumericVector R(N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      R(i) += Sigma(i,j) * x(j);
    }
  }
  for (int i = 0; i < N - 1; i++) {
    for (int j = i + 1; j < N; j++) {
      double C = 4.0 * (x(i) * R(i) - x(j) * R(j));
      for (int k = 0; k < N; k++) {
        double V = x(i) * Sigma(k,i) - x(j) * Sigma(k,j);
        g(k) += C * V;
        if ((i == k) && (j != k)) {
          g(k) += C * R(i);
        }
        if ((i != k) && (j == k)) {
          g(k) += - C * R(j);
        }
      }
    }
  }
  return 1e6 * g;
}