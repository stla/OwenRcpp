#include <Rcpp.h>
using namespace Rcpp;
// # include <cstdlib>
#include <cmath>
// using namespace std;
#include "OwenIntegrator.h"

//****************************************************************************80
double pNorm(double q){
  return R::pnorm(q,0.0,1.0,1,0);
}
double pNorm(double q);

//****************************************************************************80
// [[Rcpp::export]]
double Tha(double h, double a, double error=1e-16)
{
  double result;
  if(a < 0){
    result = -Tha(h,-a);
  }
  if(a <= 1){
	  result = OwenIntegrator::Integrate(h, a, error) / 6.2831853071795862;
  }else{
    result = (pNorm(h)+pNorm(a*h))/2-pNorm(h)*pNorm(a*h) - Tha(a*h, 1/a, error);
  }
  return result;
}

double Tha ( double h, double a, double error );

//****************************************************************************80
double dNorm(double x){
  return R::dnorm(x,0.0,1.0,0);
}
double dNorm(double x);

//****************************************************************************80
NumericVector sSequence(int n, double a, double b, double d){
  NumericVector A(n);
  NumericVector M(n);
  double sB = sqrt(b);
  M[0] = a * sB * dNorm(d*sB) * pNorm(d*a*sB);
  if(n>1){
    double sqrt2pi = 2.506628274631000502415765284811;
    A[1] = 1.0;
    M[1] = b * (d * a * M[0] + a * dNorm(d) / sqrt2pi);
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
    double C = pNorm(-delta*sB) + 2.0 * Tha(delta*sB,a);
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
    return pNorm(-delta) + sqrt2pi * sum;
  }
}


//****************************************************************************80
// [[Rcpp::export]]
double owenQ1(int nu, double t, double delta, double R){
  double a = t/sqrt(nu);
  double b = nu/(nu+t*t);
  double sB = sqrt(b);
  if(nu==1){
    double C = pNorm(R) - (delta>=0) - 2*Tha(R, (a*R-delta)/R);
    if(delta != 0){
      C += -2*Tha(delta*sB, (delta*a*b-R)/(b*delta)) + 2*Tha(delta*sB, a);
    }else{
      C += 0.5 + atan(a)/3.1415926535897931;
    }
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
    double C = pNorm(R) - (delta>=0) - 2*Tha(R, (a*R-delta)/R);
    if(delta != 0){
      C += -2*Tha(delta*sB, (delta*a*b-R)/(b*delta)) + 2*Tha(delta*sB, a);
    }else{
      C += 0.5 + atan(a)/3.1415926535897931;
    }
    return C+2*sum;
  }
}
