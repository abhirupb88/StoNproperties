#' Cumulative Density Function of Stomped Normal (StoN) Distribution
#'
#' @param x Input random variable or vector.
#' @param mu Parameter \eqn{\mu} of the distribution.
#' @param sig Parameter \eqn{\sigma} of the distribution.
#' @param k Parameter \emph{k} of the distribution.
#' @return The CDF of StoN with parameters(mu,sigma,k) at X=x.
#' @examples
#' pStoN(1,0,5,1)
#' pStoN(c(1,5,4),0,1,0.5)

#' @export
pStoN = function(x, mu, sig, k) {
  if( is.vector(x) ){
    d = matrix(x, length(x), 1)
    apply(d, 1, pStN, mu=mu, sig=sig, k=k)
  }else pStN(x,mu,sig,k)
}

pStN = function(x, mu, sig, k) {
  if( missing(mu)  ) mu  = 0.0
  if( missing(sig) ) sig = 1.0

  xs = (x-mu)/sig
  D = 2*(1-pnorm(k)+k*dnorm(k))
  if( xs<=-k ) {
    cdf = pnorm(xs)/D
  } else if( xs>-k && xs<k ){
    cdf = 0.5+(xs*dnorm(k))/D
  } else cdf = 1-pnorm(-xs)/D

  return(cdf)
}
