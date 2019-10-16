\name{Irish}
\alias{Irish}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Irish smart meter data}
\description{
Data set on residential electricity demand in Ireland, from the CER Smart Metering Project on Electricity Customer Behaviour Trial (\url{http://www.ucd.ie/issda/data/commissionforenergyregulationcer/}). This is not the full data set, but a subset of it.
} %- end description

\usage{
data(Irish)
}
%- maybe also `usage' for other objects documented here.
\format{\code{Irish} is a list, where \code{Irish$indCons} is matrix where each row is the demand for an individual household. \code{Irish$survey} is \code{data.frame} containing the following variables:
\describe{
\item{ID}{individual customer ID;}
\item{meanDem}{the mean demand of each customer;}
\item{SOCIALCLASS, OWNERSHIP, ...}{see \url{http://www.ucd.ie/issda/data/commissionforenergyregulationcer/} for details}
\code{Irish$extra} is \code{data.frame} containing the following variables:
\item{time}{progressive time counter;}
\item{toy}{the time of year from 0 (1st Jan) to 1 (31st Dec);}
\item{dow}{factor variable indicating the day of the week;}
\item{holy}{binary variable indicating holidays;}
\item{tod}{the time of day, ranging from 0 to 47, where 0 indicates the period from 00:00 to
00:30, 1 the period from 00:30 to 01:00 and so on;}
\item{temp}{the external temperature in degrees Celsius;}
}
}

\value{A list.}
\examples{
library(electBook)
data(Irish)
  
# Demand for 10 customers over 4 days
matplot(Irish$indCons[1:(48 * 4), 1:10], type = 'l')
  
# Survey information
head(Irish$survey)
  
# Temperature and calendar information
head(Irish$extra)
}
