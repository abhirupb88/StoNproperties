#' Random Sample Generation from Stomped Normal (StoN) Distribution
#'
#' @param n The number of observations.
#' @param mu Parameter \eqn{\mu} of the distribution.
#' @param sig Parameter \eqn{\sigma} of the distribution.
#' @param k Parameter \emph{k} of the distribution.
#' @return The random sample of size n from StoN with parameters(mu,sigma,k).
#' @examples
#' rStoN(1000,0,5,1)

#' @export
rStoN = function(n, mu, sig, k) {
  u = runif(n,0,1)
  D = 2*(1-pnorm(k)+k*dnorm(k))
  y = vector("numeric",n)
  for(i in 1:n){
    if( u[i]<=pnorm(-k)/D ){
      y[i] = mu+sig*qnorm(D*u[i])
    } else if( u[i]>pnorm(-k)/D && u[i]<1-pnorm(-k)/D ){
      y[i] = mu+sig*D*(u[i]-0.5)/dnorm(k)
    } else y[i] = mu-sig*qnorm(D*(1-u[i]))
  }
  return(y)
}
