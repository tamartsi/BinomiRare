prepareForMetaBRtest <-
function(d, probs, G, variant.annot, output.file = "carriers_prob_dat.csv", test = FALSE, return.result.object = FALSE){
	if (length(probs) != nrow(G)) stop("Number of individuals with disease probablity different than number of rows in genotyping matrix")
	if (sum(is.na(probs)) > 0){
		message(paste0("Removing ", sum(is.na(probs)) , " individuals with missing disease probabilityes"))
		inds.rm <- which(is.na(probs))
		probs <- probs[-inds.rm]
		G <- G[-inds.rm,, drop = F]
	}
	
	required.annot.cols <- c("variant", "chromosome", "position", "alleleA", "other.alleles")
	## check annotation columns
	for (colname.ind in 1:length(required.annot.cols)){
		if (!is.element(required.annot.cols[colname.ind], colnames(variant.annot))) stop(paste0("Variant annotation do not have column named ", required.annot.cols[colname.ind]))
	}
	
	vars.no.annot <- setdiff(colnames(G), variant.annot$variant)
	if (length(vars.no.annot) > 0){
		message(paste0("Removing ", length(vars.no.annot), " variants without annotation..."))
		inds.no.annot <- match(vars.no.annot, colnames(G))
		G <- G[,-inds.no.annot, drop = F]
		if (ncol(G) == 0) stop("No more variants to test")
	} else message("All vairants have annotation info")
	variant.inds.in.annot <- match(colnames(G), variant.annot$variant)
	
	message("Preparing information...")
	
	
	
	if (!test){
		res <- data.frame(variant = colnames(G),
							 chromosome = variant.annot$chromosome[variant.inds.in.annot], 
							 position = variant.annot$position[variant.inds.in.annot], 
							 alleleA = variant.annot$alleleA[variant.inds.in.annot], 
							 other.alleles = variant.annot$other.alleles[variant.inds.in.annot], 
							 n.carrier = NA, n.D.carrier = NA, carrier.probs = NA, stringsAsFactors = F)
	} else{ ## also test
		res <- data.frame(variant = colnames(G), 
							chromosome = variant.annot$chromosome[variant.inds.in.annot], 
							 position = variant.annot$position[variant.inds.in.annot], 
							 alleleA = variant.annot$alleleA[variant.inds.in.annot], 
							 other.alleles = variant.annot$other.alleles[variant.inds.in.annot], 
							n.carrier = NA, n.D.carrier = NA, expected.n.D.carrier = NA, pval = NA, 
							carrier.probs = NA, stringsAsFactors = F)
		}
	
	
	for (i in 1:ncol(G)){
		carrier.inds <- which(G[,i] > 0)
		res$n.carrier[i] <- length(carrier.inds)
		cur.prob.vec <- probs[carrier.inds]
		res$n.D.carrier[i] <- sum(d[carrier.inds])
		res$carrier.probs[i] <- paste(cur.prob.vec, collapse = ";")
		
		if (test){
			res$expected.n.D.carrier[i] <- sum(cur.prob.vec)			
			res$pval[i] <- poibin.midp(res$n.carrier[i], res$n.D.carrier[i], cur.prob.vec)		
		}
	}

	write.csv(res, file = output.file, row.names = FALSE)
	if (return.result.object) return(res)
}
