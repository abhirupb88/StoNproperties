#' Method of Moments Estimators of the Parameters of Stomped Normal (StN) Distribution
#'
#' @param XS Input dataset.
#' @return MoM Estimators of Parameters (mu, sig, k).
#' @examples
#' estimate.mom(rStoN(1000,0,5,1))

#' @export
estimate.mom = function(XS) {

  mu.mom = mean(XS); s2 = var(XS)
  m4 = mean((XS - mu.mom)^4)
  solve.k = function(k) s2^2*(3+(dnorm(k)*(k^3+k^5/5))/(1-pnorm(k)+k*dnorm(k)))-m4*(1+(dnorm(k)*k^3)/(3*(1-pnorm(k)+k*dnorm(k))))^2
  k.mom = try(uniroot(solve.k, c(0, 10))$root,silent = TRUE)
  if(is(k.mom,'numeric')) {
    sig.mom = sqrt(s2/(1+(dnorm(k.mom)*k.mom^3)/(3*(1-pnorm(k.mom)+k.mom*dnorm(k.mom)))))
  } else {
    sig.mom = NA; k.mom = NA
  }
  rm(solve.k,s2,m4)

  return( c(mu.mom,sig.mom,k.mom) )
}
