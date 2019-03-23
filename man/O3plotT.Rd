\name{O3plotT}
\alias{O3plotT}

\title{
Draws an Overview of Outliers (O3) plot for one method and parallel coordinate plots
}
\description{
Function for drawing Overview of Outliers (O3) plots for one method and up to 3 tolerance levels and for drawing supporting parallel coordinate plots.
}
\usage{
O3plotT(outResults, caseNames=paste0("X", 1:nrow(outResults$data)),
 sortVars=TRUE, coltxtsize=14, O3control=O3plotColours())
}
\arguments{
  \item{outResults}{
a list for each tolerance level, and within that for each variable combination, of the variables used, the indices of cases identified as outliers, and the outlier distances for all cases in the dataset.
}
  \item{caseNames}{
the ID variable used to identify the cases in an O3 plot, the default is the rownames from the dataset (so that they will then just be X7, X11, etc.)
}
  \item{sortVars}{
sort the variable columns by how often the variables occur in combinations, otherwise keep the variable order in the dataset
}
  \item{coltxtsize}{
set the size of text for column names in O3 plots (useful if there are so many columns that names overlap)
}
  \item{O3control}{
A list of colours for O3 plots. If omitted, \code{\link{O3plotColours}} gives the defaults.
}
}
\details{
This function takes the output from \code{\link{O3prep}} and draws an O3 plot with up to 3 different tolerance levels.  The default colours are khaki for the least strict tolerance level, yellow for the next, and red for the strictest.

The two parallel coordinate plots, one using the raw data and one using outlier distances, are examples of what can be done to explore the results in more detail.  If you want these plots with other highlighting then you can use outsTable with either the dataset or the Ds array to draw them using \code{\link[GGally]{ggparcoord}} from \pkg{GGally} or whatever graphics tool you prefer.

The plots produced are ggplot objects so that you can work with them—to some extent—directly.  In particular, plot margins can be set using  \code{+ theme(plot.margin = unit(c(t, r, b, l), ''cm''))}, which is useful when the cases are labelled with the caseNames option.
}

\value{
  \item{nOut}{
numbers of outliers found for the specified tolerance levels
}
  \item{gpcp}{
a parallel coordinate plot of all the data with cases ever found to be outliers coloured red
}
  \item{gO3}{
an O3 plot
}
  \item{gCombs}{
a parallel coordinate plot of the outlier distances calculated for each variable combination for the full dataset with cases found to be outliers at the strictest tolerance level coloured red
}
  \item{outsTable}{
a table of all outliers found by case, variable combination, and tolerance level.  The variable combination labels are a binary coding in the original order of the variables in the dataset.
}
  \item{Ds}{
a three-dimensional array of tolerance levels by variable combinations by cases of the outlier distances calculated.
}

}

\author{
Antony Unwin unwin@math.uni-augsburg.de
}


\seealso{\code{\link{O3plotColours}}, \code{\link[HDoutliers]{HDoutliers}} in \pkg{HDoutliers}, \code{\link[FastPCS]{FastPCS}} in \pkg{FastPCS}, \code{\link[robustX]{mvBACON}} in \pkg{robustX}, \code{\link[robustbase]{adjOutlyingness}} in \pkg{robustbase}, \code{\link[cellWise]{DDC}} in \pkg{cellWise}, \code{\link[robustbase]{covMcd}} in \pkg{robustbase}
}
\examples{
a0 <- O3prep(stackloss, method="PCS", tols=0.05, boxplotLimits=3)
a1 <- O3plotT(a0)
a1$nOut
a1$gO3

b0 <- O3prep(stackloss, method="BAC", k1=2, tols=0.01, boxplotLimits=6)
b1 <- O3plotT(b0)
b1$nOut
b1$gpcp
b1$gO3

\dontrun{
a2 <- O3prep(stackloss, method="PCS", tols=c(0.1, 0.05, 0.01), boxplotLimits=c(3, 6, 10))
a3 <- O3plotT(a2)
a3$nOut
a3$gpcp
a3$gO3
a3$outsTable
}

}

