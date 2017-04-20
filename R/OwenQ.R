#' @title First Owen Q-function
#' @description Evaluates the first Owen Q-function (integral from \eqn{0} to \eqn{R})
#' for an integer value of the degrees of freedom.
#' @param nu integer greater than \eqn{1}, the number of degrees of freedom
#' @param t finite number, positive or negative
#' @param delta finite number, positive or negative
#' @param R finite positive number, the upper bound of the integral
#' @useDynLib OwenRcpp
#' @return A number between \eqn{0} and \eqn{1}, the value of the integral from \eqn{0} to \eqn{R}.
#' @export
#' @examples
#' # OwenQ1(nu, t, delta, Inf) = ptOwen(t, nu, delta)
#' OwenQ1(4, 3, 2, 100)
#' ptOwen(3, 4, 2)
OwenQ1 <- function(nu, t, delta, R){
  if(R<0){
    stop("R must be positive.")
  }
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  if(is.infinite(t) || is.infinite(delta) || is.infinite(R)){
    stop("Parameters must be finite.")
  }
  owenQ1(nu=nu, t=t, delta=delta, R=R)
}

