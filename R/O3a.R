# Function to find the outliers for each combination for one method and one level

O3a <- function(ouF, k1=k1, K=K, mm=mm, tol=tol, boxplotLimit=boxplotLimit, n1=n1, n2=n2) {

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
        outV = robustX::mvBACON(data.frame(ouF[ , vars]),
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
        outV = cellWise::DetectDeviatingCells(ouF[ , vars],
                    DDCpars = list(tolProb = 1-tol, silent = TRUE))
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
        if (k==1) {
        mu1 <- as.numeric(outV$center)
        var1 <- as.numeric(outV$cov)
        outDist <- (ouF[ , vars]- mu1)^2/var1
        } else {
        outDist <- outV$mah
        }
        list(variables = noquote(vars),
           outlierIndices = as.vector(which(outDist > qchisq(1-tol, k))),
           outDist)
        }, simplify = FALSE)
      }
    }

# SHOULDN'T chisq df always be k?
# WHY didn't I use outV$cutoff or just outV$flagX??  Because I want to choose tol
# if (mm == "adjOutD") {
# xloop <- function(k) {
#   utils::combn(names(ouF), k, FUN = function(vars) {
#     outV = mrfDepth::adjOutl(ouF[ , vars])
#     list(variables = noquote(vars),
#          outlierIndices = as.vector(which(outV$outlyingnessX > (qchisq(1-tol, n1)**0.5)*median(outV$outlyingnessX))),
#          outDist = outV$outlyingnessX)
#     }, simplify = FALSE)
#   }
#  }

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

  if (k1 < 2 & mm %in% c("PCS", "DDC")) {
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

  if (mm %in% c("HDo", "BAC", "adjOut", "MCD")| k1 > 1) {
    suspects <- lapply(k1:K, xloop)
    suspects <- unlist(suspects, recursive = FALSE)
    }
  list(outM = suspects)
}
