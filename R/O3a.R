# Function to find the outliers for each combination for one method and one level

O3a <- function(ouF, k1=k1, K=K, mm=mm, Alpha=Alpha, Coef=Coef, n1=n1, n2=n2) {

# This function is required to suppress unnecessary output
DetectDeviantCells <- function(...) {
    utils::capture.output(x <- cellWise::DetectDeviatingCells(...))
    x
  }
  
# Choice of method----------------------------
  if (mm=="HDo") {
    xloop <- function(k) {
      utils::combn(names(ouF), k, FUN = function(vars) {
        list(variables = noquote(vars),
        outlierIndices = as.vector(as.integer(HDoutliers::HDoutliers(ouF[ , vars],
        alpha=Alpha)))) }, simplify = FALSE)
      }
    }
  if (mm=="PCS") {
    xloop <- function(k) {
      utils::combn(names(ouF), k, FUN = function(vars) {
        list(variables = noquote(vars),
        outlierIndices = setdiff(1:n2, FastPCS::FastPCS(ouF[ , vars],
        alpha=1-Alpha)$best)) }, simplify = FALSE)
      }
    }
  if (mm=="BAC") {
    xloop <- function(k) {
      utils::combn(names(ouF), k, FUN = function(vars) {
        list(variables = noquote(vars),
        outlierIndices = as.vector(which(robustX::mvBACON(ouF[ , vars],
        alpha=Alpha, verbose=FALSE)$subset==FALSE)))
        }, simplify = FALSE)
      }
    }
  if (mm=="adjOut") {
  xloop <- function(k) {
    utils::combn(names(ouF), k, FUN = function(vars) {
      list(variables = noquote(vars),
           outlierIndices = as.vector(which(robustbase::adjOutlyingness(
           ouF[ , vars], alpha.cutoff=1-Alpha)$nonOut==FALSE)))
    }, simplify = FALSE)
  }
  }
 if (mm=="DDC") {
    xloop <- function(k) {
      utils::combn(names(ouF), k, FUN = function(vars) {
        list(variables = noquote(vars),
        outlierIndices =
                    DetectDeviantCells(ouF[, vars],
                    DDCpars=list(tolProb=1-Alpha))$indrows)
        }, simplify = FALSE)
      }
    }

# Method for 1-d outliers--------------
  out1d <- function(j) {
    outx <- boxplot.stats(ouF[ ,j], coef=Coef, do.conf=FALSE)$out
    list(variables=noquote(names(ouF)[[j]]),
    outlierIndices=which(ouF[ ,j] %in% outx))
    }

  if (k1 < 2 & !(mm %in% c("HDo", "adjOut"))) {
    #Find 1-d outliers
    suspectsA <- lapply(1:n1, out1d)
    #Find higher-d outliers
    if (K<2) {
      suspects <- suspectsA
      } else {
      suspectsB <- lapply((k1+1):K, xloop)
      #Combine 1-d and higher-d in one list
      suspects <- c(suspectsA, unlist(suspectsB, recursive = FALSE))
      }
    }
  if (mm %in% c("HDo", "adjOut")| k1 > 1) {
    suspects <- lapply(k1:K, xloop)
    suspects <- unlist(suspects, recursive = FALSE)
    }
  list(outM=suspects)
}
