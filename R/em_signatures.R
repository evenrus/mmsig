#' Expectation maximization algorithm to estimate the signature contribution
#'
#' @param sigts.defn mutational signature definitions
#' @param mut.freqs mutational profile of sample
#' @param max.iter maximum number of iterations
#' @param dbg boolean whether to print or not
#'
#' @return estimated mutational signature contributions
#'
em_signatures = function(sigts.defn, mut.freqs, max.iter, dbg) {
  nof.sigts.defn <- ncol(sigts.defn)
  alpha = stats::runif(nof.sigts.defn); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
  for (i in 1:max.iter) {
    contr = t(array(alpha, dim=c(nof.sigts.defn,96))) * sigts.defn
    probs = contr/array(rowSums(contr), dim=dim(contr))
    probs[is.na(probs)] = 0
    probs = probs * mut.freqs
    old_alpha = alpha
    alpha = colSums(probs)/sum(probs)
    if (sum(abs(alpha-old_alpha))<1e-5) {
      break
    }
  }
  spit(dbg, "em: exit iteration: %d", i)

  return( alpha )
}
