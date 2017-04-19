#' @title Owen T-function
#' @description Evaluates the Owen T-function.
#' @param h numeric scalar
#' @param a numeric scalar
#' @return A number between 0 and 1.
#' @export
#' @useDynLib OwenRcpp
#' @importFrom Rcpp evalCpp
#' @importFrom stats pnorm
#' @examples
#' # OwenT(h,a) = OwenT(-h,a)
#' OwenT(2,1) == OwenT(-2,1)
#' # OwenT(0,a) = atan(a)/2pi
#' a <- runif(1, -1000, 1000)
#' OwenT(0,a) - atan(a)/(2*pi)
#' # OwenT(h,1) = Phi(h)(1-Phi(h))/2
#' h <- runif(1, -3, 3)
#' OwenT(h,1) - pnorm(h)*(1-pnorm(h))/2
#' # OwenT(h,Inf) = (1-Phi(|h|))/2 :
#' OwenT(1,10000) - (1-pnorm(abs(1)))/2
#' OwenT(1,Inf) == (1-pnorm(abs(1)))/2
OwenT <- function(h, a){
  if(is.infinite(a)){
    return(sign(a)*(1-pnorm(abs(h)))/2)
  }else if(is.infinite(h)){
    return(0)
  }else{
    return(tfn(h, a))
  }
}

#' @title Owen T-function for ratio values of the parameters
#' @description Evaluates the Owen T-function at \code{h=h1/h2} and \code{a=a1/a2}.
#' @param h1,h2,a1,a2 numeric scalars
#' @return A number between 0 and 1, the value of \eqn{T(h1/h2,a1/a2)}
#' @export
#' @useDynLib OwenRcpp
#' @importFrom Rcpp evalCpp
#' @examples
#' OwenT(2,2) - OwenTratios(2, 1, 2, 1)
#' # T(h,Inf) = (1-Phi(|h|))/2 :
#' OwenTratios(3, 1, 1, 0) - (1-pnorm(abs(3)))/2
#' # T(Inf,a) = 0:
#' OwenTratios(1, 0, 3, 1)
OwenTratios <- function(h1, h2, a1, a2){
  tha(h1, h2, a1, a2)
}
