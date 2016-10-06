\name{BRtest}
\alias{BRtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tests the association between a set of variants and disease status.
}
\description{
Tests a set of variants represented in the numberical matrix G, 
given disease status and previously estimated probabilities of disease. 
}
\usage{
BRtest(d, probs, G)
}
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
}
}
\details{
 For each variant, calculates the expected number of disease individuals 
 among the carriers, and calculates p-value from the BinomiRare test of 
 association between the variant and the disease status.
}
\value{
The function returns a data frame with variant name, the number of 
carriers of the variant, number of diseased carriers, the expected
number of diseased carriers (according to the supplied probabilities),
and p-value. 

}
\author{
Tamar Sofer
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
}

