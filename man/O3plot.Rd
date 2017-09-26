\name{O3plot}
\alias{O3plot}

\title{
Draws the Overview of Outliers (O3) plot and a parallel coordinates plot
}
\description{
The main function for preparing and drawing the Overview of Outliers (O3) plot and a supporting parallel coordinates plot.
}
\usage{
O3plot(ouF, k1=1, K=ncol(ouF), mm="HDo", Alphas=0.05, Coefs=6,
       outColours  = c("khaki", "yellow", "red", "salmon",
       "chocolate", "red", "slategray1", "slategray2",
       "slategray3", "slategray4", "red"))
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
method(s) used for identifying outliers (up to five can be used)
}
  \item{Alphas}{
outlier level(s) (up to three can be used if only one method is specified, otherwise just one)
}
  \item{Coefs}{
boxplot limit(s) if method does not apply for a single variable (there should be the same number of Coefs as Alphas)
}
  \item{outColours}{
A set of colours for the three kinds of plot.  There are 11 in all and the defaults are khaki, yellow, red, salmon, chocolate, red, slategray1, slategray2, slategray3, slategray4, red.

The first three (1-3) are for plots with three different significance levels;
the next three (4-6) are for plots comparing results for two methods;
the final five (7-11) are for plots combining results of up to five methods.
}
}
\details{
To check outliers for all possible combinations of variables choose k1=1 and K=number of variables in the dataset.

The optional methods are "HDo" \command{HDoutliers} (from  \pkg{HDoutliers}), "PCS" \command{FastPCS} (\pkg{FastPCS}), "BAC" \command{mvBACON} (\pkg{robustX}), "adjOut" \command{adjOutlyingness} (\pkg{robustbase}), "DDC" \command{DectectDeviatingCells} (\pkg{Cellwise}).  References for all these methods can be found on their help pages, linked below.

If only one method is specified, then up to three significance levels (Alphas) and three boxplot limits (Coefs) can be specified.

\code{Alpha} is the parameter determining outlyingness and should be set low, as in \command{HDoutliers} and \command{mvBACON}.  For the other methods, which take a high \code{Alpha}, \code{(1-Alpha)} is used.  Using the same \code{Alpha} for all methods does not make them directly comparable.  The next version of the package will allow \code{Alpha} to be set differently for each method.

Methods "HDo" and "adjOut" can analyse single variables.  For the other methods boxplot limits are used for single variables and any case > (Q3 + Coef*IQR) or < (Q1 - Coef*IQR) is classed an outlier, where \code{Coef} is the limit specified.
}
\value{
  \item{nrout}{
numbers of outliers found by method(s)/significance levels
}
  \item{gpcp}{
a parallel coordinate plot of all the data with cases ever found to be outliers coloured red
}
  \item{gO3}{
an O3 plot
}

}

\author{
Antony Unwin unwin@math.uni-augsburg.de
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link[HDoutliers]{HDoutliers}} in \pkg{HDoutliers}, \code{\link[FastPCS]{FastPCS}} in \pkg{FastPCS}, \code{\link[robustX]{mvBACON}} in \pkg{robustX}, \code{\link[robustbase]{adjOutlyingness}} in \pkg{robustbase}, \code{\link[cellWise]{DetectDeviatingCells}} in \pkg{cellWise}
}
\examples{
a0 <- O3plot(stackloss, mm="PCS", Alphas=0.05, Coefs=3)
a0$nrout
a0$gpcp
a0$gO3

b0 <- O3plot(stackloss, mm="BAC", k1=2, Alphas=0.01, Coefs=6)
b0$nrout
b0$gpcp
b0$gO3

\dontrun{
a1 <- O3plot(stackloss, mm="PCS", Alphas=c(0.1, 0.05, 0.01), Coefs=c(3, 6, 10))
a1$nrout
a1$gpcp
a1$gO3

a2 <- O3plot(stackloss, mm=c("HDo", "BAC", "DDC"), Alphas=0.05, Coefs=3)
a2$gpcp
a2$gO3
a2$nrout
}

# It is advisable with large datasets to check the number of outliers identified (nrout)
# before drawing graphics.  Occasionally methods find very many outliers.
\dontrun{ 
data(diamonds, package="ggplot2")
ouF <- diamonds[1:5000, c(1, 5, 6, 8:10)]
pPa <- O3plot(ouF, mm="PCS", Alphas=0.01, Coefs=10)
pPa$nrout
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
