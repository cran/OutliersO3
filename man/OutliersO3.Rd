\name{OutliersO3}
\alias{OutliersO3-package}
\docType{package}
\title{
OutliersO3: displays and compares potential outliers
}

\description{
Up to six different methods can be used to identify potential outliers in a dataset.  An O3 (Overview of Outliers) plot and supporting parallel coordinate plots are drawn.
}

\details{
An Overview of Outliers (O3) plot is the result of checking for potential outliers for every possible combination of (numeric) variables.  It shows which cases are identified as outliers for each combination for which any outliers are found.

The available methods are "HDo" \command{HDoutliers} (from  \pkg{HDoutliers}), "PCS" \command{FastPCS} (\pkg{FastPCS}), "BAC" \command{mvBACON} (\pkg{robustX}), "adjOut" \command{adjOutlyingness} (\pkg{robustbase}), "DDC" \command{DectectDeviatingCells} (\pkg{Cellwise}), "MCD" \command{covMcd} (\pkg{robustbase}).  Outlier tolerance levels can be set individually for each method.

Plots can be drawn for a single method with up to 3 tolerance levels, for two methods showing which cases are found by one or both, or for up to 6 methods showing how often a case is identified as an outlier.

If a method cannot be used for single variables, then boxplot limits defined by the argument boxplotLimits are used.

If only one method is specified, then up to three tolerance levels (tols) and boxplot limits (boxplotLimits) can be specified.
}

\author{
Antony Unwin unwin@math.uni-augsburg.de

Thanks are due to Bill Venables for some key coding advice.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
