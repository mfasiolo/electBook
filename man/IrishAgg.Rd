\name{IrishAgg}
\alias{IrishAgg}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Irish aggregate demand data}
\description{
Data set on aggregate electricity demand from the customer in the \code{Irish} data set (see \code{?Irish}).
} %- end description

\usage{
data(IrishAgg)
}
%- maybe also `usage' for other objects documented here.
\format{\code{IrishAgg} is a \code{data.frame}, most of its variables are explained in \code{?Irish} with the exception of:
\describe{
\item{dem}{total electricity demand across the customers in the \code{Irish} data set;}
\item{demLag}{total electricity demand at lag 1;}
\item{temp95}{exponential smooth of \code{temp}, that is \code{temp95[i] = a*temp[i] + (1-a)*temp95[i]} with a = 0.05;}
}
}

\examples{
library(electBook)
data(IrishAgg)

head(IrishAgg)
  
# Aggregate demand over one week
plot(IrishAgg$dem[1:(48*7)], type = 'l')
}
