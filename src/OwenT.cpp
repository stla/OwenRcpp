#include <Rcpp.h>
#include <cmath>
#include <cfloat>

// [[Rcpp::export]]
double OwenT01(double h, double a, int jmax, double cutpoint){
  const double twopi = 6.2831853071795862;
  if(h > cutpoint){
    return atan(a) * exp(-0.5 * (h*h) * a / atan(a)) *
      (1 + 0.00868 * pow(h*a,4)) / twopi;
  }
  double cumsum = 1.0;
  const double K = exp(-h*h/2);
  // double matrjk[jmax+1];
  // matrjk[0] = 1-K;
  // double jk[jmax+1];
  // jk[0] = 1.0;
  double crossprod = (1-K)*a;
  int i;
  for(i=1; i<=jmax; i++){
    cumsum += exp(2*i*log(h) - i*log(2.0) - lgamma(i+1));
    // jk[i] = 1.0/(2.0*i+1);
    // if(i%2==1){
    //   jk[i] = -jk[i];
    // }
    // matrjk[i] = (1.0 - K*cumsum)*jk[i];
    // crossprod += matrjk[i]*pow(a,2*i+1);
    double jk = 1.0/(2.0*i+1);
    if(i%2==1){
      jk = -jk;
    }
    crossprod += (1.0 - K*cumsum)*jk*pow(a,2*i+1);
  }
  return (atan(a)-crossprod)/twopi;
}
double OwenT01(double h, double a, int jmax, double cutpoint);

//****************************************************************************80
double pNorm(double q){
  return R::pnorm(q,0.0,1.0,1,0);
}
double pNorm(double q);

//****************************************************************************80
// [[Rcpp::export]]
double OwenT(double h, double a, int jmax, double cutpoint){
  double absa = fabs(a);
  if(absa == 0.0){
    return 999.0;
  }
  double absh = fabs(h);
  if(absa >= DBL_MAX){
    return R::sign(a) * pNorm(-absh) / 2;
  }
  if(absa <= 1.0){
    return OwenT01(absh, absa, jmax, cutpoint);
  }
  return R::sign(a) * (pNorm(absh)/2 + pNorm(absa*absh) * (0.5 - pNorm(absh)) -
               OwenT01(absa*absh, 1/absa, jmax, cutpoint));
}
