\name{BRtestMeta}
\alias{BRtestMeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Applies BinomiRare (BR) test in meta-analysis based on previously 
prepared files from multiple studies.
}
\description{
Given a folder with files of data prepared for BR meta-analysis, 
merge the files together and applies the BR test.
}
\usage{
BRtestMeta(folder, recursive = TRUE, output.file = "BR_meta_results.csv", 
	error.log.file = "BR_meta_inconsistencies.txt", return.result.object = FALSE)
}
\arguments{
  \item{folder}{
A folder in which the files from the various studies are at.
Only .csv files are considered.  
}
  \item{recursive}{
Should folders with the main folders be looked at?
}
  \item{output.file}{
File with meta-analysis test results.
}
  \item{error.log.file}{
File for which the meta-analysis function prints identified errors
(e.g. non-matching alleles, etc).
}
  \item{return.result.object}{
Should an R data.frame with results be returned, in addition the printing
of the results to file?
}
}

\value{
prints results to file, and if return.result.object == TRUE, also 
returns a data.frame with columns: variant (variant name), 
n.carrier (total across studies), n.D.carrier (total across studies),
chromosome, position, alleleA, other.alleles, n.studies 
(the number of studied that had individuals carrying the variant), 
pval (BR p-value), and expected.n.D.carrier (the expected number of 
diseased individuals among the carriers).
}

\author{
Tamar Sofer
}
\note{
other.allele is still not used. It is assume that the same 
effect allele (alleleA) is tested by all studies that privde it. If there's no
match in annotation for a given vairant, information about the variant 
is thrown away. 
}


\examples{
require(poibin)	
dir.create(file.path(getwd(), "Files_for_meta"))
setwd(file.path(getwd(), "Files_for_meta"))


variant.annot <- data.frame(variant = paste0("var_", 1:500), chromosome = 1, 
	position = 1:500, alleleA = sample(c("A", "C", "G", "T"), 
	size = 500, replace = TRUE), other.alleles = NA, stringsAsFactors = FALSE)

for (i in 1:nrow(variant.annot)){
	variant.annot$other.alleles[i] <- sample(setdiff(c("A", "C", "G", "T"), 
		variant.annot$alleleA[i]), size = 1)
}


n <- 10000
effect.size <- 0
n.in.pop <- n/5
pop.risk <- c(rep(-2.6, n.in.pop), rep(-2.2, n.in.pop), rep(-2, n.in.pop), 
			rep(-1.8, n.in.pop) , rep(-1.6, n.in.pop) )
pop.inds <- c(rep(1, n.in.pop), rep(2, n.in.pop), rep(3, n.in.pop), 
			rep(4, n.in.pop) , rep(5, n.in.pop) )

### for each population, prepare information for meta-analysis and write to file.
x <- rnorm(n, sd = 0.02)
x <- pmax(x, 0)

g <- rbinom(n, size = 2, prob = x)  ##
G <- matrix(rbinom(n*100, size = 1, prob = 0.001), nrow = n)  ### another 100 null variants
G <- cbind(g, G)
colnames(G) <- paste0("var_", 1:ncol(G))
rownames(G) <- paste0("person_", 1:nrow(G))
p <- expit(pop.risk +   g*effect.size + 20*x)
d <- rbinom(n, 1, p)

for (pop in 1:5){
	d.pop <- d[which(pop.inds == pop)]
	G.pop <- G[which(pop.inds == pop),]
	x.pop <- x[which(pop.inds == pop)]
	prob.mod <- glm(d.pop ~ x.pop, family = binomial)
	prob.d <- expit(predict(prob.mod))
	prepareForMetaBRtest(d.pop, prob.d, G.pop, variant.annot, output.file = 
		file.path(getwd(), paste0("carriers_prob_study_", pop, ".csv")), 
		test = FALSE, return.result.object = FALSE)	
}


############# Now meta-analyze results from files:
results.folder <- getwd() 

res <- BRtestMeta(folder = results.folder, return.result.object = TRUE)


}
