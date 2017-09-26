\name{sortO3}
\alias{sortO3}

\title{
Sorting the outlier indicator matrix ready for plotting
}
\description{
Rows and columns are sorted by counts within structures.
}
\usage{
sortO3(z1, n1=n1, nz=nz)
}
\arguments{
  \item{z1}{
a matrix of variables and outliers in columns and variable combinations in rows
}
  \item{n1}{
number of variables in the dataset
}
  \item{nz}{
number of cases found to be outliers for at least one variable combination
}
}
\details{
z1 only includes columns for cases identified as outliers by some combination and similarly only combinations which find at least one outlier.
Rows are ordered by numbers of variables within combinations and then by numbers of outliers.  
Columns are ordered in two groups, first a column for each variable, then a column for each outlier case (there is a gap column between these two groups of columns).
Variables come first that appear in least combinations; cases are ordered by how often they are identified as outliers.
}
\value{
z1p, a long version of the sorted matrix
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
