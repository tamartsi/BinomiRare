BRtest <-
function(d, probs, G){
	if (length(probs) != nrow(G)) stop("Number of individuals with disease probablity different than number of rows in genotyping matrix")
	if (sum(is.na(probs)) > 0){
		message(paste0("Removing ", sum(is.na(probs)) , " individuals with missing disease probabilityes"))
		inds.rm <- which(is.na(probs))
		probs <- probs[-inds.rm]
		G <- G[-inds.rm,, drop = F]
	}
	
	message("Testing variants...")
	
	res <- data.frame(variant = colnames(G), n.carrier = NA, n.D.carrier = NA, expected.n.D.carrier = NA, pval = NA, stringsAsFactors = F)
	
	for (i in 1:ncol(G)){
		carrier.inds <- which(G[,i] > 0)
		res$n.carrier[i] <- length(carrier.inds)
		cur.prob.vec <- probs[carrier.inds]
		res$expected.n.D.carrier[i] <- sum(cur.prob.vec)
		res$n.D.carrier[i] <- sum(d[carrier.inds])
		
		res$pval[i] <- poibin.midp(n.carrier = res$n.carrier[i], n.D.carrier = res$n.D.carrier[i], prob.vec = cur.prob.vec)		 
	}
	
	return(res)	
}
