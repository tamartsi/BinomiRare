BRtestMeta <-
function(folder, recursive = TRUE, output.file = "BR_meta_results.csv", error.log.file = "BR_meta_inconsistencies.txt" , return.result.object = FALSE){
	files <- list.files(folder, full.names = TRUE, recursive = recursive)
	if (length(files) == 0) stop("No files in folder")
	inds.no.csv <- setdiff(1:length(files), grep(".csv", files, fixed = TRUE))
	if (length(inds.no.csv) > 0){
		message(paste0("Removing ", length(inds.no.csv), "files that are not in csv format"))
		files <- files[-inds.no.csv]
		if (length(files) == 0) stop("No files left after removing non-csv format files")
	}
	
	message("\t Reading information from files...")
	required.cols <- c("variant", "n.carrier", "n.D.carrier", "carrier.probs", "chromosome", "position", "alleleA", "other.alleles")
	### initialize with the first file:
	i <- 1
	dat.all <- read.csv(files[i], stringsAsFactors = F)
	
	## check columns
	for (colname.ind in 1:length(required.cols)){
		if (!is.element(required.cols[colname.ind], colnames(dat.all))) stop(paste0("File ", files[i], " do not have column named ", required.cols[colname.ind]))
	}
	dat.all$alleleA <- toupper(dat.all$alleleA)
	dat.all$other.alleles <- toupper(dat.all$other.alleles)
	dat.all <- dat.all[,required.cols]
	
	inds.no.carriers <- which(dat.all$n.carrier < 1)
	if (length(inds.no.carriers) > 0) 	dat.all <- dat.all[-inds.no.carriers,]
	
	dat.all$n.studies <- 1
	
	if (length(files) == 1) message("Only one file available, using information from one study...")

	if (length(files) > 1){
		
		for (i in 2:length(files)){
			dat.i <- read.csv(files[i], stringsAsFactors = F)
			
			## check columns
			for (colname.ind in 1:length(required.cols)){
				if (!is.element(required.cols[colname.ind], colnames(dat.i))) stop(paste0("File ", files[i], " do not have column named ", required.cols[colname.ind]))
			}
			dat.i$alleleA <- toupper(dat.i$alleleA)
			dat.i$other.alleles <- toupper(dat.i$other.alleles)
			
			## remove variants with no carriers:
			inds.no.carriers <- which(dat.i$n.carrier < 1)
			if (length(inds.no.carriers) > 0) 	dat.i <- dat.i[-inds.no.carriers,]
	
			#### find SNPs that are common to dat.all and dat.i, check match in chromosome, position, and alleleA, and update dat.all:
			common.vars <- intersect(dat.all$variant, dat.i$variant)
			common.vars.inds.in.dat.all <- match(common.vars, dat.all$variant)
			common.vars.inds.in.dat.i <- match(common.vars, dat.i$variant)
			
			## check that annotation is the same:
			inds.no.annot.match <- c(which(dat.all$chromosome[common.vars.inds.in.dat.all] != dat.i$chromosome[common.vars.inds.in.dat.i]),
									which(dat.all$position[common.vars.inds.in.dat.all] != dat.i$position[common.vars.inds.in.dat.i]),
									which(dat.all$alleleA[common.vars.inds.in.dat.all] != dat.i$alleleA[common.vars.inds.in.dat.i]))
			if (length(inds.no.annot.match) > 0){
				inds.no.annot.match <- unique(inds.no.annot.match)	
				message("Some variants have different annotations, and are removed from some studies. Please check log file.")
				write.table(paste0("Variants ", dat.all$variant[common.vars.inds.in.dat.all[inds.no.annot.match]], 
						" have chromosomes ", dat.all$chromosome[common.vars.inds.in.dat.all[inds.no.annot.match]], 
							"in some files, but chromosomes ", dat.i$chromosome[common.vars.inds.in.dat.i[inds.no.annot.match]], 
							"in ", files[i]), file = error.log.file, append = TRUE, col.names = FALSE, row.names = FALSE)
				
				write.table(paste0("Variants ", dat.all$variant[common.vars.inds.in.dat.all[inds.no.annot.match]], 
						" have positions ", dat.all$position[common.vars.inds.in.dat.all[inds.no.annot.match]], 
							"in some files, but positions ", dat.i$position[common.vars.inds.in.dat.i[inds.no.annot.match]], 
							"in ", files[i]), file = error.log.file, append = TRUE, col.names = FALSE, row.names = FALSE)
							
				write.table(paste0("Variants ", dat.all$variant[common.vars.inds.in.dat.all[inds.no.annot.match]], 
						" have alleleAs ", dat.all$alleleA[common.vars.inds.in.dat.all[inds.no.annot.match]], 
							"in some files, but alleleAs ", dat.i$alleleA[common.vars.inds.in.dat.i[inds.no.annot.match]], 
							"in ", files[i]), file = error.log.file, append = TRUE, col.names = FALSE, row.names = FALSE)
							
									
				common.vars.inds.in.dat.all <- common.vars.inds.in.dat.all[-inds.no.annot.match]
				common.vars.inds.in.dat.i <- common.vars.inds.in.dat.i[-inds.no.annot.match]
			}
			

			
			dat.all[common.vars.inds.in.dat.all ,]$n.carrier <- dat.all[common.vars.inds.in.dat.all ,]$n.carrier + dat.i[common.vars.inds.in.dat.i,]$n.carrier
			dat.all[common.vars.inds.in.dat.all ,]$n.D.carrier <- dat.all[common.vars.inds.in.dat.all ,]$n.D.carrier + dat.i[common.vars.inds.in.dat.i,]$n.D.carrier
			dat.all[common.vars.inds.in.dat.all ,]$carrier.probs <- paste0(dat.all[common.vars.inds.in.dat.all ,]$carrier.probs, ";", dat.i[common.vars.inds.in.dat.i,]$carrier.probs)
			dat.all[common.vars.inds.in.dat.all ,]$n.studies <- dat.all[common.vars.inds.in.dat.all ,]$n.studies + 1
			
			#### find SNPs that are new in dat.i, and add to dat.all:
			dat.i.vars <- setdiff(dat.i$variant, common.vars)
			
			if (length(dat.i.vars) > 0){
				no.match.inds.in.dat.i <- match(dat.i.vars, dat.i$variant)
				dat.to.add <- dat.i[no.match.inds.in.dat.i, required.cols]
				dat.to.add$n.studies <- 1
				dat.all <- rbind(dat.all, dat.to.add)
			}
			
		}
	}	

	message("\t Finished combining information, testing association...")

	dat.all$pval <- NA
	dat.all$expected.n.D.carrier <- NA
	for (var.ind in 1:nrow(dat.all)){
		cur.prob.vec <- strsplit(dat.all$carrier.probs[var.ind], split = ";")[[1]]
		cur.prob.vec <- as.numeric(cur.prob.vec)
		dat.all$expected.n.D.carrier[var.ind] <- sum(cur.prob.vec)
		
		dat.all$pval[var.ind] <- poibin.midp(dat.all$n.carrier[var.ind], dat.all$n.D.carrier[var.ind], cur.prob.vec)		
	}
	
	write.csv(dat.all, file = output.file, row.names = FALSE)
	if (return.result.object) {
		dat.all <- dat.all[,-grep("carrier.probs", colnames(dat.all))]
		return(dat.all)
	}
	
}
