\name{BinomiRare-package}
\alias{BinomiRare-package}
\alias{BinomiRare}
\docType{package}
\title{
\packageTitle{BinomiRare}
}
\description{
\packageDescription{BinomiRare}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{BinomiRare}
\packageIndices{BinomiRare}
Functions to test whether the number of diseased individuals among carriers
 of rare variants is consistent with expectation, where the expectation
 is determined via any pre-specified disease probability model. It also included
 functions to export data to be used in meta-analysis, and to perform meta-analysis
 across several studies. 
 In this version, not all possible quality check are implemented (e.g. the meta-analysis
 does not flip alleles, and uses information for a given variant 
 between studies only if their effect alleles match).
}
\author{
\packageAuthor{BinomiRare}

Maintainer: Tamar Sofer \packageMaintainer{BinomiRare}
}
\references{
~~ Literature or other references for background information ~~
}

\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
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
#res <- prepareForMetaBRtest(d, prob.d, G, variant.annot[-1,], output.file =
#  "carriers_prob_dat.csv", test = TRUE, return.result.object = TRUE)
### CHECK: if the variant.annot data frame did not have a necessary column:
#res <- prepareForMetaBRtest(d, prob.d, G, variant.annot[,-1], output.file = 
# "carriers_prob_dat.csv", test = TRUE, return.result.object = TRUE)


##########  Example 2:  simulate multiple data set, and then meta-analyze them. 

########## First simulate data:
dir.create(file.path(getwd(), "Files_for_meta"))
setwd(file.path(getwd(), "Files_for_meta"))

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
		file.path(getwd(), paste0("carriers_prob_study_", pop, ".csv")), test = FALSE, 
		return.result.object = FALSE)	
}


############# Now meta-analyze results from files:
results.folder <- getwd() 

res <- BRtestMeta(folder = results.folder, return.result.object = TRUE)

}
