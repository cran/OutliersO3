\name{Election2005}
\alias{Election2005}
\docType{data}
\title{ Election2005 data
%%   ~~ data name/kind ... ~~
}
\description{
A data set for the German 'Bundestag' election of 2005. It includes information about the elections
in 2005 and in 2002 separately for each of the 299 constituencies and also demographic and other information about the constituencies themselves.  
}
\usage{data(Election2005)}
\format{
  A data frame with 299 observations on 70 variables.
  The variables of the data set are:
\itemize{
 \item 1-4: general information about the constituencies (ID and name of district, Bundesland)
 \item 5-40: demographic and economic information
 \item 41 - 48: general information about the elections in 2009 and 2005
 \itemize{
 \item \code{WBerechE}: number of eligible voters
 \item \code{WE}: votes cast
 \item \code{UngZE}: invalid second votes
 \item \code{GulZE}: valid second votes
 }
 \item 49- 70: results of the five biggest parties in 2009 and 2005 (smaller parties are summarized in the variable \code{Rest}). Everything with a 'V' or 'v' at the end is from the election in 2005. 'ze' and 'zv'  at the end refer to second votes. 
 }
  
}
\details{
This dataset has been taken from the package \pkg{mbgraphic} of Katrin Grimm.
}

\source{
Original source was \url{http://www.bundeswahlleiter.de}.
}

\examples{
dim(Election2005)
str(Election2005)
}
\keyword{datasets}
