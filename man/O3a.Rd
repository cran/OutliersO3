\name{O3a}
\alias{O3a}

\title{
Find the outliers for each combination for one method and one level
}
\description{
Potential outliers are identified for all combinations made up of from k1 to K variables in a dataset.
}
\usage{
O3a(ouF, k1=k1, K=K, 
   mm, Alpha, Coef,
   n1=n1, n2=n2)
}
\arguments{
  \item{ouF}{
dataset to be checked for outliers
}
  \item{k1}{
lowest number of variables in a combination
}
  \item{K}{
highest number of variables in a combination
}
  \item{mm}{
method(s) used for identifying outliers
}
  \item{Alpha}{
outlier level
}
  \item{Coef}{
boxplot limit if method cannot be used for a single variable
}
  \item{n1}{
number of variables in the dataset
}
  \item{n2}{
number of cases in the dataset
}
}
\details{
To check outliers for all possible combinations of variables set k1=1 and K=number of variables in the dataset.

The optional methods are "HDo" \command{HDoutliers} (from  \pkg{HDoutliers}), "PCS" \command{FastPCS} (\pkg{FastPCS}), "BAC" \command{mvBACON} (\pkg{robustX}), "adjOut" \command{adjOutlyingness} (\pkg{robustbase}), "DDC" \command{DectectDeviatingCells} (\pkg{Cellwise}).

If only one method is specified, then up to three significance levels (Alphas) and three boxplot limits (Coefs) can be chosen.

\code{Alpha} is the parameter determining outlyingness and should be set low, as in \command{HDoutliers} and \command{mvBACON}.  For the other methods, which take a high \code{Alpha}, \code{(1-Alpha)} is used.

Methods "HDo" and "adjOut" can analyse single variables.  For the other methods boxplot limits are used for single variables and any case > (Q3 + Coef*IQR) or < (Q1 - Coef*IQR) is classed an outlier, where \code{Coef} is the limit specified.
}
\value{
outM (list of combinations and outliers found)
}

\author{
Antony Unwin unwin@math.uni-augsburg.de
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link[HDoutliers]{HDoutliers}} in \pkg{HDoutliers}, \code{\link[FastPCS]{FastPCS}} in \pkg{FastPCS}, \code{\link[robustX]{mvBACON}} in \pkg{robustX}, \code{\link[robustbase]{adjOutlyingness}} in \pkg{robustbase}, \code{\link[cellWise]{DetectDeviatingCells}} in \pkg{cellWise}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
