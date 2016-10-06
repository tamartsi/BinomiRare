\name{poibin.midp}
\alias{poibin.midp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculated the mid-p values based on the Poisson-Binomial distribution
}
\description{
Given a vector of disease probabilities and a number of diseased individual,
Calculated the mid-p values based on the Poisson-Binomial distribution
}
\usage{
poibin.midp(n.carrier, n.D.carrier, prob.vec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.carrier}{
 The number of carriers of a rare variant
}
  \item{n.D.carrier}{
 The number of diseased carriers of a rare variant.
 n.D.carrier cannot be larger than n.carrier.
}
  \item{prob.vec}{
vector of disease probabilities of the carriers. 
}
}

\value{
a single numeric variable - a p-value for the test that the number 
n.D.carrier is consistent with prob.vec.
}

\author{
Tamar Sofer
}
\note{
Althouth n.carrier is not strictly needed, but is useful for quality 
checks, especially when meta-analyzing. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
require(poibin)
	n <- 100

	### simulation under the null
	g <- rbinom(n, 2, 0.1)
	x <- rnorm(n)
	p <- expit(-2.3 + x)
	d <- rbinom(n, size = 1, prob = p)
	
	mod <- glm(d ~ x, family = "binomial")
	prob.d <- expit(predict(mod))
	poibin.midp(n.carrier = sum(g >0 ), n.D.carrier = sum(g*d > 0), 
		prob.vec = prob.d[which(g>0)] )
	
	##### under the alternative:	
	p <- expit(-2.3 + x + g)
	d <- rbinom(n, size = 1, prob = p)
	
	mod <- glm(d ~ x, family = "binomial")
	prob.d <- expit(predict(mod))
	poibin.midp(n.carrier = sum(g >0 ), n.D.carrier = sum(g*d > 0),
		prob.vec = prob.d[which(g>0)] )
	
}
