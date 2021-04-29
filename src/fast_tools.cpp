#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double fsd(SEXP x)
{
  double M = 0.0;
  double S = 0.0;
  int k = 1;
  int n = LENGTH(x);

  if (Rf_isReal(x))
  {
    double* px = REAL(x);

    for (int i = 0; i < n; ++i)
    {
      double tmpM = M;
      M += (px[i] - tmpM) / k;
      S += (px[i] - tmpM) * (px[i] - M);
      k++;
    }
  }
  else if (Rf_isInteger(x))
  {
    int* px = INTEGER(x);

    for (int i = 0; i < n; ++i)
    {
      double tmpM = M;
      M += (px[i] - tmpM) / k;
      S += (px[i] - tmpM) * (px[i] - M);
      k++;
    }
  }
  else
  {
    Rf_error("x is not INTEGER or REAL");
  }

  return std::sqrt(S / (k-2));
}

// [[Rcpp::export]]
double ffsd(SEXP x)
{
  double stdDev = 0;
  double sumAll = 0;
  double sumAllQ = 0;
  int n = LENGTH(x);

  if (Rf_isReal(x))
  {
    double* px = REAL(x);

    for (int i = 0; i < n; ++i)
    {
      double x = px[i];
      sumAll += x;
      sumAllQ += x * x;
    }
  }
  else if (Rf_isInteger(x))
  {
    int* px = INTEGER(x);

    for (int i = 0; i < n; ++i)
    {
      double x = px[i];
      sumAll += x;
      sumAllQ += x * x;
    }
  }
  else
  {
    Rf_error("x is not INTEGER or REAL");
  }

  double var = (sumAllQ - (sumAll * sumAll) / n) * (1.0 / (n - 1));
  return std::sqrt(var);
}
