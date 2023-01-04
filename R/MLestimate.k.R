#' MLE of Width Parameter \emph{k} of Stomped Normal (StN) Distribution
#'
#' @param data Input random vector.
#' @param mu Parameter \eqn{\mu} of the distribution.
#' @param sig Parameter \eqn{\sigma} of the distribution.
#' @return MLE of Parameter \emph{k} of the distribution.
#' @examples
#' MLestimate.k(rStoN(1000,0,5,1),0,5)

#' @export
MLestimate.k = function(data, mu, sig) {
  dummy = unique(sort(abs(data-mu)/sig))
  dummy = unique(round(dummy,2))
  k.est = 0; diff = 1
  while(diff > 0.000001){
    dumLik = sapply(dummy, function(x) sum(log(dStoN(data, mu, sig, x))))
    idx = max(which(dumLik==max(dumLik)))
    if (idx == 1 | idx == length(dummy)) {
      k.est = dummy[idx]; diff = 0
    } else {
      k.tmp = dummy[idx]; k.low = dummy[idx-1]; k.top = dummy[idx+1]
      dummy = seq(k.low, k.top, length.out = 100)
      diff  = abs(sum(log(dStoN(data, mu, sig, k.tmp))) - sum(log(dStoN(data, mu, sig, k.est))))
      k.est = k.tmp; rm(k.tmp, k.low, k.top)
    }
    rm(dumLik)
  }
  rm(dummy, diff)
  return(k.est)
}
