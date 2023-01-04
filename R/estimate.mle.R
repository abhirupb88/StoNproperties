#' Likelihood Iterative Estimators (LIEs) of the Parameters of Stomped Normal (StN) Distribution
#'
#' @param XS Input dataset.
#' @return LI Estimators of Parameters (mu, sig, k).
#' @examples
#' estimate.mle(rStoN(1000,0,5,1))

#' @export
estimate.mle = function(XS) {

  mom = estimate.mom(XS); mu0 = mom[1]
  if (is.na(mom[2])) sig0 = sqrt(var(XS)) else sig0 = mom[2]
  if (is.na(mom[3])) k0 = 0 else k0 = mom[3]

  max_iter = 500; diff = 1.0; iter = 1
  while(diff > 0.09 & iter < max_iter) {
    k1 = MLestimate.k(XS, mu0, sig0)
    sig1 = optim(par=sig0, function(sig) sum(log(dStoN(XS,mu0,sig,k1))), method = "Brent",
                 lower = 0.0001, upper = 50, control=list(trace=2,fnscale=-1))$par
    mu1  = mean(XS[abs((XS-mu0)/sig1)<=k1])

    if (sum(log(dStoN(XS,mu1,sig1,k1)))-sum(log(dStoN(XS,mu0,sig1,k1))) > -0.1 & !is.na(mu1)) mu0 = mu1
    diff = sum(log(dStoN(XS,mu0,sig1,k1)))-sum(log(dStoN(XS,mu0,sig0,k0)))
    sig0 = sig1; k0 = k1; iter = iter+1
    rm(mu1,sig1,k1)
  }
  rm(mom,max_iter,diff,iter)

  return( c(mu0, sig0, k0) )
}
