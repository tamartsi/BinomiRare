\name{logit}
\alias{logit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The logit function
}
\description{
Given a numerical variable x, returns it logit transfomation
(inverse of expit)
}
\usage{
logit(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A numerical vector
}
}
\details{
logit(x) = log(x/(1-x))
This is related to the logistic model of disease probability
}
\value{
a numerical vector 
}
\author{
Tamar Sofer
}

\examples{
	x <- rnorm(100)
	all(abs(logit(expit(x)) - x) < 1e-7)
}
