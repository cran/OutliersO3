\name{O3plotM}
\alias{O3plotM}

\title{
Draws an Overview of Outliers (O3) plot for more than one method and parallel coordinate plots
}
\description{
Function for drawing Overview of Outliers (O3) plots for comparing outlier methods and for drawing supporting parallel coordinate plots.
}
\usage{
O3plotM(outResults, caseNames=paste0("X", 1:nrow(outResults$data)),
 sortVars=TRUE, coltxtsize=14, O3control=O3plotColours())
}
\arguments{
  \item{outResults}{
a list for each method, and within that for each variable combination, of the variables used, the indices of cases identified as outliers, and the outlier distances for all cases in the dataset.
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
This function takes the output from \code{\link{O3prep}} and draws an O3 plot.  If there are only two methods, then the default colours are red if both methods identify the case as an outlier for a variable combination and blue or green if only one method does.  With more than two methods the default colours are red if all methods identify the outlier, orange if all but one method do, and shades of slategray otherwise.

The two parallel coordinate plots, one using the raw data and one using outlier distances, are examples of what can be done to explore the results in more detail.  If you want these plots with other highlighting then you can use outsTable with either the dataset or the Cs array to draw them using \code{\link[GGally]{ggparcoord}} from \pkg{GGally} or whatever graphics tool you prefer.

}

\value{
  \item{nOut}{
numbers of outliers found by each method
}
  \item{gpcp}{
a parallel coordinate plot of the dataset with cases ever found to be outliers coloured red
}
  \item{gO3}{
an O3 plot
}
  \item{gO3x}{
an O3 plot for three or more methods in which outliers identified by only one method for a variable combination are ignored.
}
  \item{gMethods}{
a parallel coordinate plot of the outlier distances calculated by each method for the full dataset with cases ever found to be outliers coloured red
}
  \item{outsTable}{
a table of all outliers found by case, variable combination, and method.  The variable combination labels are a binary coding in the original order of the variables in the dataset.
}
  \item{Cs}{
a three-dimensional array of methods by variable combinations by cases of the outlier distances calculated.
}

}

\author{
Antony Unwin unwin@math.uni-augsburg.de
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{O3plotColours}}, \code{\link[HDoutliers]{HDoutliers}} in \pkg{HDoutliers}, \code{\link[FastPCS]{FastPCS}} in \pkg{FastPCS}, \code{\link[robustX]{mvBACON}} in \pkg{robustX}, \code{\link[robustbase]{adjOutlyingness}} in \pkg{robustbase}, \code{\link[cellWise]{DetectDeviatingCells}} in \pkg{cellWise}, \code{\link[robustbase]{covMcd}} in \pkg{robustbase}
}
\examples{
c1 <- O3prep(stackloss, k1=2, method=c("HDo", "BAC"), tolHDo=0.025, tolBAC=0.01)
c2 <- O3plotM(c1)
c2$nOut
c2$gpcp
c2$gO3

\dontrun{
b1 <- O3prep(stackloss, method=c("HDo", "BAC", "DDC"), tolHDo=0.025, tolBAC=0.01, tolDDC=0.05)
b2 <- O3plotM(b1)
b2$nOut
b2$gpcp
b2$gO3
b2$outsTable
}

# It is advisable with large datasets to check the number of outliers identified (nOut)
# before drawing graphics.  Occasionally methods find very many outliers.
\dontrun{ 
data(diamonds, package="ggplot2")
data <- diamonds[1:5000, c(1, 5, 6, 8:10)]
pPa <- O3prep(data, method=c("PCS", "adjOut"), tolPCS=0.01, toladj=0.01, boxplotLimits=10)
pPx <- O3plotM(pPa)
pPx$nOut
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
