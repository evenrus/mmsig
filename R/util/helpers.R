## Helper functions for mmsig

spit = function (dbg, mess, ...) {
    if (dbg) { print( sprintf(mess,...) ) }
}

spat = function (dbg, name, var) {
    if (dbg) {print(name); print(var)}
}

oppstrand = function(x) {   # comp() complements dna sequences
    trin1 = paste(comp(rev(strsplit(substr(x,1,3),"")[[1]]), forceToLower=F), collapse="")
    trin2 = paste(comp(rev(strsplit(substr(x,5,7),"")[[1]]), forceToLower=F), collapse="")
    return(sprintf("%s>%s",trin1,trin2))
}

# EM algorithm to estimate the signature contribution
em_signatures = function(sigts.defn,mut.freqs,max.iter,dbg) {
    nof.sigts.defn <- ncol(sigts.defn)
    alpha = runif(nof.sigts.defn); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
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