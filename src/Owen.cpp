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
  int i;
  double r[NG] = {
    0.1477621, 
    0.1346334, 
    0.1095432, 
    0.0747257, 
    0.0333357 };
  double r1;
  double r2;
  double rt;
  double tp = 0.159155;
  double tv1 = 1.0E-35;
  double tv2 = 15.0;
  double tv3 = 15.0;
  double tv4 = 1.0E-05;
  double u[NG] = {
    0.0744372, 
    0.2166977, 
    0.3397048,
    0.4325317, 
    0.4869533 };
  double value;
  double x1;
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
  if ( tv3 <= log ( 1.0 + fxs ) - xs * fxs )
  {
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
  for ( i = 0; i < NG; i++ )
  {
    r1 = 1.0 + fxs * pow ( 0.5 + u[i], 2 );
    r2 = 1.0 + fxs * pow ( 0.5 - u[i], 2 );

    rt = rt + r[i] * ( exp ( xs * r1 ) / r1 + exp ( xs * r2 ) / r2 );
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
  double absa;
  double ah;
  double c1;
  double c2;
  double ex;
  double g;
  double gah;
  double gh;
  double h;
  double lam;
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
    lam = fabs ( a * h );
    ex = exp ( - lam * lam / 2.0 );
    g = R::pnorm ( lam, 0.0, 1.0, 1, 0 );
    c1 = ( ex / lam + sqrt ( twopi ) * ( g - 0.5 ) ) / twopi;
    c2 = ( ( lam * lam + 2.0 ) * ex / lam / lam / lam
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
    absa = fabs ( a );

    if ( absa <= 1.0 )
    {
      value = tfn ( h, a );
      return value;
    }

    ah = absa * h;
    gh = R::pnorm ( h, 0.0, 1.0, 1, 0 );
    gah = R::pnorm ( ah, 0.0, 1.0, 1, 0 );
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
    int k;
    if(n>2){
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
