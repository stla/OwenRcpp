# include <Rcpp.h>
using namespace Rcpp;
# include <cstdlib>
# include <cmath>
using namespace std;

//****************************************************************************80
// [[Rcpp::export]]
double tfn ( double x, double fx )
{
# define NG 5

  double fxs;
  // double r[NG] = {
  //   0.1477621,
  //   0.1346334,
  //   0.1095432,
  //   0.0747257,
  //   0.0333357 };
  // double r1;
  // double r2;
  double rt;
  double tp = 0.159155;
  double tv1 = 1.0E-35;
  double tv2 = 15.0;
  double tv3 = 15.0;
  // double tv4 = 1.0E-05;
  // double u[NG] = {
  //   0.0744372,
  //   0.2166977,
  //   0.3397048,
  //   0.4325317,
  //   0.4869533 };
  double value;
  // double x1;
  double x2;
  double xs;
//
//  Test for X near zero.
//
  if ( fabs ( x ) < tv1 )
  {
    value = tp * atan ( fx );
    return value;
  }
//
//  Test for large values of abs(X).
//
  if ( tv2 < fabs ( x ) )
  {
    value = 0.0;
    return value;
  }
//
//  Test for FX near zero.
//
  if ( fabs ( fx ) < tv1 )
  {
    value = 0.0;
    return value;
  }
//
//  Test whether abs ( FX ) is so large that it must be truncated.
//
  xs = - 0.5 * x * x;
  x2 = fx;
  fxs = fx * fx;
//
//  Computation of truncation point by Newton iteration.
//
  if ( tv3 <= log1p ( fxs ) - xs * fxs )
  {
    double tv4 = 1.0E-05;
    double x1;
    x1 = 0.5 * fx;
    fxs = 0.25 * fxs;

    for ( ; ; )
    {
      rt = fxs + 1.0;

      x2 = x1 + ( xs * fxs + tv3 - log ( rt ) )
      / ( 2.0 * x1 * ( 1.0 / rt - xs ) );

      fxs = x2 * x2;

      if ( fabs ( x2 - x1 ) < tv4 )
      {
        break;
      }
      x1 = x2;
    }
  }
//
//  Gaussian quadrature.
//
  rt = 0.0;
  {
    // double r1;
    // double r2;
    double r[NG] = {
      0.1477621,
      0.1346334,
      0.1095432,
      0.0747257,
      0.0333357 };
    double u[NG] = {
      0.0744372,
      0.2166977,
      0.3397048,
      0.4325317,
      0.4869533 };
    int i;
    for ( i = 0; i < NG; i++ )
    {
      double r1 = 1.0 + fxs * pow ( 0.5 + u[i], 2 );
      double r2 = 1.0 + fxs * pow ( 0.5 - u[i], 2 );
      rt = rt + r[i] * ( exp ( xs * r1 ) / r1 + exp ( xs * r2 ) / r2 );
    }
  }

  value = rt * x2 * tp;

  return value;
# undef NG
}

double tfn ( double x, double fx );

//****************************************************************************80
// [[Rcpp::export]]
double tha ( double h1, double h2, double a1, double a2 )
{
  double a;
  // double absa;
  double ah;
  // double c1;
  // double c2;
  // double ex;
  double g;
  // double gah;
  // double gh;
  double h;
  // double lam;
  double twopi = 6.2831853071795864769;
  double value;

  if ( h2 == 0.0 )
  {
    value = 0.0;
    return value;
  }

  h = h1 / h2;

  if ( a2 == 0.0 )
  {
    g = R::pnorm(h, 0.0, 1.0, 1, 0);

    if ( h < 0.0 )
    {
      value = g / 2.0;
    }
    else
    {
      value = ( 1.0 - g ) / 2.0;
    }

    if ( a1 < 0.0 )
    {
      value = - value;
    }
    return value;
  }

  a = a1 / a2;

  if ( fabs ( h ) < 0.3 && 7.0 < fabs ( a ) )
  {
    double lam = fabs ( a * h );
    double ex = exp ( - lam * lam / 2.0 );
    g = R::pnorm ( lam, 0.0, 1.0, 1, 0 );
    double c1 = ( ex / lam + sqrt ( twopi ) * ( g - 0.5 ) ) / twopi;
    double c2 = ( ( lam * lam + 2.0 ) * ex / lam / lam / lam
      + sqrt ( twopi ) * ( g - 0.5 ) ) / ( 6.0 * twopi );
    ah = fabs ( h );
    value = 0.25 - c1 * ah + c2 * ah * ah * ah;
    if ( a < 0.0 )
    {
      value = - fabs ( value );
    }
    else
    {
      value = fabs ( value );
    }
  }
  else
  {
    double absa = fabs ( a );

    if ( absa <= 1.0 )
    {
      value = tfn ( h, a );
      return value;
    }

    ah = absa * h;
    double gh = R::pnorm ( h, 0.0, 1.0, 1, 0 );
    double gah = R::pnorm ( ah, 0.0, 1.0, 1, 0 );
    value = 0.5 * ( gh + gah ) - gh * gah - tfn ( ah, 1.0 / absa );

    if ( a < 0.0 )
    {
      value = - value;
    }
  }
  return value;
}

double tha ( double h1, double h2, double a1, double a2 );

//****************************************************************************80
NumericVector sSequence(int n, double a, double b, double d){
  NumericVector A(n);
  NumericVector M(n);
  double sB = sqrt(b);
  M[0] = a * sB * R::dnorm(d*sB, 0.0, 1.0, 0) * R::pnorm(d*a*sB, 0.0, 1.0, 1, 0);
  if(n>1){
    double sqrt2pi = 2.506628274631000502415765284811;
    A[1] = 1.0;
    M[1] = b * (d * a * M[0] + a * R::dnorm(d, 0.0, 1.0, 0) / sqrt2pi);
    if(n>2){
      int k;
      for ( k = 2; k < n; k++ ){
        A[k] = 1.0 / (k-1.0) / A[k-1];
        M[k] = (k-1.0)/k * b * (A[k-1] * d * a * M[k-1] + M[k-2]);
      }
    }
  }
  return M;
}

NumericVector sSequence(int n, double a, double b, double d);

//****************************************************************************80
// [[Rcpp::export]]
double pStudent(double q, int nu, double delta){
  double a = q/sqrt(nu);
  double b = nu/(nu+q*q);
  double sB = sqrt(b);
  if(nu % 2 == 1){
    double C = R::pnorm(-delta*sB, 0.0, 1.0, 1, 0) + 2.0 * tfn(delta*sB,a);
    if(nu == 1){
      return C;
    }else{
      int i;
      NumericVector M = sSequence(nu, a, b, delta);
      double sum = 0.0;
      for(i=1; i<nu-1; i+=2){
        sum += M[i];
      }
      return C + 2.0*sum;
    }
  }else{
    int i;
    NumericVector M = sSequence(nu, a, b, delta);
    double sum = 0.0;
    for(i=0; i<nu-1; i+=2){
      sum += M[i];
    }
    double sqrt2pi = 2.506628274631000502415765284811;
    return R::pnorm(-delta, 0.0, 1.0, 1, 0) + sqrt2pi * sum;
  }
}

double dNorm(double x){
  return R::dnorm(x,0.0,1.0,0);
}
double dNorm(double x);

double pNorm(double q){
  return R::pnorm(q,0.0,1.0,1,0);
}
double pNorm(double q);

//****************************************************************************80
// [[Rcpp::export]]
double owenQ1(int nu, double t, double delta, double R){
  double a = t/sqrt(nu);
  double b = nu/(nu+t*t);
  double sB = sqrt(b);
  if(nu==1){
    double C = pNorm(R) - 2*tha(R, 1, a*R-delta, R) -
      2*tha(delta*sB, 1, delta*a*b-R, b*delta) + 2*tfn(delta*sB, a) -
      (delta>=0);
    return C;
  }
  int n = nu-1;
  NumericVector H(n);
  NumericVector M(n);
  H[0] = -dNorm(R) * pNorm(a*R-delta);
  M[0] = a*sB*dNorm(delta*sB)*(pNorm(delta*a*sB)-pNorm((delta*a*b-R)/sB));
  if(nu >= 3){
    H[1] = R * H[0];
    M[1] = b*(delta*a*M[0] + a*dNorm(delta*sB)*(dNorm(delta*a*sB)-dNorm((delta*a*b-R)/sB)));
    if(nu >= 4){
      NumericVector A(n);
      NumericVector L(n-2);
      A[0] = 1;
      A[1] = 1;
      L[0] = a * b * R * dNorm(R) * dNorm(a*R-delta) / 2;
      int k;
      for(k=2; k<n; k++){
        A[k] = 1.0/k/A[k-1];
      }
      if(nu >= 5){
        for(k=1; k<n-2; k++){
          L[k] = A[k+2] * R * L[k-1];
        }
      }
      for(k=2; k<n; k++){
        H[k] = A[k] * R * H[k-1];
        M[k] = (k-1.0)/k * b * (A[k-2] * delta * a * M[k-1] + M[k-2]) - L[k-2];
      }
    }
  }
  if(nu % 2 == 0){
    double sqrt2pi = 2.506628274631000502415765284811;
    double sum = 0.0;
    int i;
    for(i=0; i<nu-1; i+=2){
      sum += M[i]+H[i];
    }
    return pNorm(-delta) + sqrt2pi * sum;
  }else{
    double sum = 0.0;
    int i;
    for(i=1; i<nu-1; i+=2){
      sum += M[i]+H[i];
    }
    double C = pNorm(R) - 2*tha(R, 1, a*R-delta, R) -
      2*tha(delta*sB, 1, delta*a*b-R, b*delta) + 2*tfn(delta*sB, a) -
      (delta>=0);
    return C+2*sum;
  }
}
