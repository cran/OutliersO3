\name{O3prep}
\alias{O3prep}

\title{
Identify outliers for different combinations of variables
}
\description{
Check the dataset and parameters prior to analysis.  Identify outliers for the variable combinations and methods/tolerance levels specified.  Prepare input for the two plotting functions \code{O3plotT} and \code{O3plotM}.
}
\usage{
O3prep(data, k1=1, K=ncol(data), method="HDo", tols=0.05, boxplotLimits=c(6, 10, 12),
       tolHDo=0.05, tolPCS=0.01, tolBAC=0.001, toladj=0.05, tolDDC=0.01, tolMCD=0.000001)
}
\arguments{
  \item{data}{
dataset to be checked for outliers
}
  \item{k1}{
lowest number of variables in a combination
}
  \item{K}{
highest number of variables in a combination
}
  \item{method}{
method(s) used for identifying outliers (up to six can be used)
}
  \item{tols}{
outlier tolerance level(s) when only one method is specified, up to three can be used.  For consistent use of the argument, it is transformed for some of the methods.  See details below of how the argument is applied for each approach.  
}
  \item{boxplotLimits}{
up to three boxplot limits are used, matching the number of tolerance levels, if a method does not apply for a single variable.
}
  \item{tolHDo}{
an individual outlier tolerance level for the HDoutliers method.  The default in \pkg{HDoutliers}, alpha, is 0.05.
}
  \item{tolPCS}{
an individual outlier tolerance level for the FastPCS method.  This equals (1-alpha) for the argument in \pkg{FastPCS}, where the default is 0.5.
}
  \item{tolBAC}{
an individual outlier tolerance level for the mvBACON method.  The default for alpha in \pkg{robustX} is 0.95.  This seems high, but it is divided by n, the dataset size.
}
  \item{toladj}{
an individual outlier tolerance level for the adjOutlyingness method.  This equals (1-alpha.cutoff) for the argument in \pkg{robustbase}, where the default is 0.75.
}
  \item{tolDDC}{
an individual outlier tolerance level for the DDC method.  This equals (1-tolProb) for the argument in \pkg{cellWise}, where the default is 0.99.
}
  \item{tolMCD}{
an individual outlier tolerance level for the covMcd method.  The default is 0.025 (based on the help page for plot.mcd in \pkg{robustbase}). This is NOT the alpha argument in \command{covMcd}, which is used for determining subset size and set to 0.9 in OutliersO3.
}

}

\details{
To check outliers for all possible combinations of variables choose k1=1 and K=number of variables in the dataset (the default).

The optional methods are "HDo" \command{HDoutliers} (from  \pkg{HDoutliers}), "PCS" \command{FastPCS} (\pkg{FastPCS}), "BAC" \command{mvBACON} (\pkg{robustX}), "adjOut" \command{adjOutlyingness} (\pkg{robustbase}), "DDC" \command{DDC} (\pkg{Cellwise}), "MCD" \command{covMcd} (\pkg{robustbase}).  References for all these methods can be found on their help pages, linked below.  (Note that \pkg{Cellwise} has renamed its function \command{DetectDeviatingCells}.  Since version 2.1.0 \command{DDC} is used instead.)

If only one method is specified, then up to three tolerance levels (tols) and three boxplot limits (boxplotLimits) can be specified.

\code{tol} is the argument determining outlyingness and should be set low, as in \command{HDoutliers} and \command{mvBACON}, where it is called alpha, and in \command{covMcd}.  For the other methods \code{(1-tol)} is used.  In \command{DDC} the argument is called tolProb. Using the same tolerance level for all methods does not make them directly comparable, which is why it is recommended to set them individually when drawing a comparative O3 plot.  The defaults suggested on the methods' help pages mostly found too many outliers and so other defaults have been set.  Users need to decide for themselves, possibly dependent on the dataset they are analysing.

Methods "HDo", "mvBACON", "adjOut", and "MCD" can analyse single variables.  For the other methods boxplot limits are used for single variables and any case > (Q3 + boxplotLimit*IQR) or < (Q1 - boxplotLimit*IQR) is classed an outlier, where \code{boxplotLimit} is the limit specified.
}
\value{
  \item{data}{
the dataset analysed
}
  \item{nw}{
the number of variable combinations analysed
}
  \item{mm}{
the outlier methods used
}
  \item{tols}{
the individual tolerance levels for the outlier methods used (if more than one), otherwise up to 3 tolerance levels used for one method
}
  \item{outList}{
a list for each method/tolerance level, and within that for each variable combination, of the variables used, the indices of cases identified as outliers, and the outlier distances for all cases in the dataset.
}

}

\author{
Antony Unwin unwin@math.uni-augsburg.de
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link[HDoutliers]{HDoutliers}} in \pkg{HDoutliers}, \code{\link[FastPCS]{FastPCS}} in \pkg{FastPCS}, \code{\link[robustX]{mvBACON}} in \pkg{robustX}, \code{\link[robustbase]{adjOutlyingness}} in \pkg{robustbase}, \code{\link[cellWise]{DDC}} in \pkg{cellWise}, \code{\link[robustbase]{covMcd}} in \pkg{robustbase}
}
\examples{
a0 <- O3prep(stackloss, method="PCS", tols=0.05, boxplotLimits=3)

b0 <- O3prep(stackloss, method=c("BAC", "adjOut"), k1=2, tols=0.01, boxplotLimits=6)

\dontrun{
a1 <- O3prep(stackloss, method="PCS", tols=c(0.1, 0.05, 0.01), boxplotLimits=c(3, 6, 10))

b1 <- O3prep(stackloss, method=c("HDo", "BAC", "DDC"), tolHDo=0.025, tolBAC=0.01, tolDDC=0.05)
}

}

