\name{expit}
\alias{expit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The expit function
}
\description{
returns the expit of a variable x expit(x) = exp(x)/(1 + exp(x))
}
\usage{
expit(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	a numerical vector
}
}
\details{
applies the expit transformation
}
\value{
expit(x)
}
\author{
Tamar Sofer
}
\note{
 Useful for calculating disease probabilities from logistic regression.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
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
}
