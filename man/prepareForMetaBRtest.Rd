\name{prepareForMetaBRtest}
\alias{prepareForMetaBRtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepares a file with information required for meta-analysis of BinomiRare
test.
}
\description{
Based on vectors of disease status, pre-computed disease probabilities, 
a matrix of genetic variants, and variant annotation, writes a file
with information for meta-anlaysis. 
}
\usage{
prepareForMetaBRtest(d, probs, G, variant.annot, output.file = "carriers_prob_dat.csv", 
	test = FALSE, return.result.object = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
	An n vector of disease status
}
  \item{probs}{
	An n vector of (estimated) disease probabilities corresponding 
	to the individuals with disease status. 
}
  \item{G}{
	An nxp matrix of counts of (assumed rare) p variants genotpyed on n individuals. 
	Columns correspond to variants, rows to individuals. 
}  \item{variant.annot}{
A data-frame with annotation of the genetic variants under considerations.
It has to have column named "variant" (variant name), "chromosome", "position", 
 "alleleA" (effect allele), and "other.alleles" (the alleles that aren't tested). 
}
  \item{output.file}{
Name of file to print information to. 
}
  \item{test}{
In addition to printing information to file, should the BR test be performed?
}
  \item{return.result.object}{
If return.result.object is true, R will return a data.frame with results,
in addition to printing them to file. 
}
}

\value{
a matrix is printed to the file output.file.
In addition, if requested a data frame is returned.
Returned/printed columns are variant, chromosome, position, alleleA,
other.alleles, n.carrier, n.D.carrier,  and if
test == TRUE, also expected.n.D.carrier (expected number
of diseased carriers), and pval (the p-value from the BinomiRare test).
}
\author{
Tamar Sofer
}
\note{
other.alleles are currently not used, but in the future there will 
be added functionality that may flip alleles if needed, or return
information and alleles do not match, or when triallelic SNPs are used. 
}



\examples{
require(poibin)

##########  Example 1:  a single data set.
########## Simulate data
n <- 10000
effect.size <- 1
pop.risk <- -2.6

x <- rnorm(n, sd = 0.01)
x <- pmax(x, 0)
	
g <- rbinom(n, size = 2, prob = x)  ## one causal variant, x is a confounder
G <- matrix(rbinom(n*100, size = 1, prob = 0.001), nrow = n)  ### another 100 null variants
G <- cbind(g, G)
colnames(G) <- paste0("var_", 1:ncol(G))
rownames(G) <- paste0("person_", 1:nrow(G))
p <- expit(pop.risk +   g*effect.size + 20*x)
d <- rbinom(n, 1, p) 
names(d) <- paste0("person_", 1:nrow(G))


########### Now that we have outcome d, genotypes G and a covariate x: 
########### Estimate disease probability model

prob.mod <- glm(d ~ x, family = binomial)
prob.d <- expit(predict(prob.mod))

system.time(res <- BRtest(d, prob.d, G)) ### super quick


######### If we wanted to contribute this data to meta-analysis
######### here we need annotation as well!
### simulate annotation:
variant.annot <- data.frame(variant = paste0("var_", 1:500), chromosome = 1, 
	position = 1:500, alleleA = sample(c("A", "C", "G", "T"), size = 500, replace = TRUE), 
	other.alleles = NA, stringsAsFactors = FALSE)
	
for (i in 1:nrow(variant.annot)){
	variant.annot$other.alleles[i] <- sample(setdiff(c("A", "C", "G", "T"), 
		variant.annot$alleleA[i]), size = 1)
}

res <- prepareForMetaBRtest(d, prob.d, G, variant.annot, output.file = 
	"carriers_prob_dat.csv", test = TRUE, return.result.object = TRUE)
### CHECK: if we had missing annotation for variant 1:
# res <- prepareForMetaBRtest(d, prob.d, G, variant.annot[-1,], 
# output.file = "carriers_prob_dat.csv", test = TRUE, return.result.object = TRUE)
### CHECK: if the variant.annot data frame did not have a necessary column:
#res <- prepareForMetaBRtest(d, prob.d, G, variant.annot[,-1], output.file = 
# "carriers_prob_dat.csv", test = TRUE, return.result.object = TRUE)

}

