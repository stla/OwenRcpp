#' @title Student CDF with integer number of degrees of freedom
#' @description Cumulative distribution function of the noncentrel Student
#' distribution with an integer number of degrees of freedom.
#' @param q quantile
#' @param nu integer greater than \eqn{1}, the number of degrees of freedom
#' @param delta noncentrality parameter
#' @return Numeric value, the CDF evaluated at \code{q}.
#' @export
#' @useDynLib OwenRcpp
#' @examples
#' ptOwen(2, 3) - pt(2, 3)
#' ptOwen(2, 3, delta=1) - pt(2, 3, ncp=1)
ptOwen <- function(q, nu, delta=0){
  if(isNotPositiveInteger(nu)){
    stop("`nu` must be an integer >=1.")
  }
  if(is.infinite(q) || is.infinite(delta)){
    stop("Parameters must be finite.")
  }
  pStudent(as.double(q), as.integer(nu), as.double(delta))
}