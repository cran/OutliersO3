# Function to find the outliers for each combination for one method and one level

O3a <- function(ouF, k1=k1, K=K, mm=mm, tol=tol, boxplotLimit=boxplotLimit, n1=n1, n2=n2) {

# This function is required to suppress unnecessary output
DetectDeviantCells <- function(...) {
    utils::capture.output(x <- cellWise::DetectDeviatingCells(...))
    x
  }
  
# Choice of method----------------------------
  if (mm == "HDo") {
    xloop <- function(k) {
      utils::combn(names(ouF), k, FUN = function(vars) {
        list(variables = noquote(vars),
        outlierIndices = as.vector(as.integer(HDoutliers::HDoutliers(ouF[ , vars],
        alpha = tol))),
        outDist = rep(NA, n2)) 
        }, simplify = FALSE)
      }
    }
  if (mm == "PCS") {
    xloop <- function(k) {
      utils::combn(names(ouF), k, FUN = function(vars) {
        outV = FastPCS::FastPCS(ouF[ , vars], alpha = 1-tol)
        list(variables = noquote(vars),
        outlierIndices = setdiff(1:n2, outV$best),
        outDist = outV$distance) 
        }, simplify = FALSE)
      }
    }
  if (mm == "BAC") {
    xloop <- function(k) {
      utils::combn(names(ouF), k, FUN = function(vars) {
        outV = robustX::mvBACON(ouF[ , vars],
        alpha = tol, verbose = FALSE, allowSingular = TRUE)
        list(variables = noquote(vars),
        outlierIndices = as.vector(which(outV$subset == FALSE)),
        outDist = outV$dis)
        }, simplify = FALSE)
      }
    }
  if (mm == "adjOut") {
  xloop <- function(k) {
    utils::combn(names(ouF), k, FUN = function(vars) {
      outV = robustbase::adjOutlyingness(ouF[ , vars], alpha.cutoff = 1-tol)
      list(variables = noquote(vars),
           outlierIndices = as.vector(which(outV$nonOut == FALSE)),
           outDist = outV$adjout)
    }, simplify = FALSE)
  }
  }
 if (mm == "DDC") {
    xloop <- function(k) {
      utils::combn(names(ouF), k, FUN = function(vars) {
        outV = DetectDeviantCells(ouF[ , vars],
                    DDCpars = list(tolProb = 1-tol))
        list(variables = noquote(vars),
          outlierIndices = outV$indrows,
          outDist = outV$Ti)
         }, simplify = FALSE)
      }
    }
if (mm == "MCD") {
    xloop <- function(k) {
      utils::combn(names(ouF), k, FUN = function(vars) {
        outV = robustbase::covMcd(ouF[ , vars], alpha = 0.9)
        list(variables = noquote(vars),
           outlierIndices = as.vector(which(outV$mah > qchisq(1-tol, n1))),
           outDist = outV$mah)
        }, simplify = FALSE)
      }
    }

# if (mm == "sHDo") {
#    xloop <- function(k) {
#      utils::combn(names(ouF), k, FUN = function(vars) {
#        list(variables = noquote(vars),
#        outlierIndices = as.vector(as.integer(stray::find_HDoutliers(ouF[ , vars],
#        alpha = tol))),
#        outDist = rep(NA, n2)) 
#        }, simplify = FALSE)
#      }
#    }

# Would using chisq be what the CS people intend?
# if (mm == "hdb") {
#    xloop <- function(k) {
#      utils::combn(names(ouF), k, FUN = function(vars) {
#        outV = dbcan::hdbscan(as.matrix(ouF[ , vars]), minPts = 3),
#        list(variables = noquote(vars),
#        outlierIndices =  as.vector(which(outV$outlier_scores > qchisq(1-tol, n1))),
#        outDist = outV$outlier_scores)
#        }, simplify = FALSE)
#      }
#    }
       
# Method for 1-d outliers--------------
  out1d <- function(j) {
    outx <- boxplot.stats(ouF[ ,j], coef = boxplotLimit, do.conf = FALSE)$out
    list(variables = noquote(names(ouF)[[j]]),
    outlierIndices = which(ouF[ ,j] %in% outx),
    outDist = abs((ouF[ ,j]-median(ouF[ ,j]))/(IQR(ouF[ ,j]))))
    }

#  if (k1 < 2 & !(mm %in% c("HDo", "adjOut", "sHDo"))) {
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
#  if (mm %in% c("HDo", "adjOut", "sHDo")| k1 > 1) {
  if (mm %in% c("HDo", "adjOut")| k1 > 1) {
    suspects <- lapply(k1:K, xloop)
    suspects <- unlist(suspects, recursive = FALSE)
    }
  list(outM = suspects)
}
