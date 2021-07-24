#' Probability Density Function of Stomped Normal (StoN) Distribution
#'
#' @param x Input random variable or vector.
#' @param mu Parameter \eqn{\mu} of the distribution.
#' @param sig Parameter \eqn{\sigma} of the distribution.
#' @param k Parameter \emph{k} of the distribution.
#' @return The PDF of StoN with parameters(mu,sigma,k) at X=x.
#' @examples
#' dStoN(1,0,5,1)
#' dStoN(c(1,5,4),0,1,0.5)

#' @export
dStoN = function(x, mu, sig, k) {

  xs = (x-mu)/sig
  D = 2*(1-pnorm(k)+k*dnorm(k))
  xs[xs>-k & xs<k] = k
  pdf = dnorm(xs)/(D*sig)

  return(pdf)
}
