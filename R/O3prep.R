# quiets concerns of R CMD check
if(getRversion() >= "2.15.1")  utils::globalVariables(c("ID",
            "outlierIndices", "pID", "psN", "psNx", "sB", "s1B",
            "sN", "sumR", "sumS", "sumV", "xsumR", "oh"))

# Main function--------------

O3prep <- function(data, k1=1, K=ncol(data), method="HDo",
            tols=0.05, boxplotLimits=c(6, 10, 12),
     tolHDo=0.05, tolPCS=0.01, tolBAC=0.001, toladj=0.05, tolDDC=0.01, tolMCD=0.000001) {
  ouF <- data.frame(data) #in case the dataset is a tibble

# Put methods in a standard order and check they are OK
#  mz <- c("HDo", "PCS", "BAC", "adjOut", "DDC", "MCD", "sHDo")
  mz <- c("HDo", "PCS", "BAC", "adjOut", "DDC", "MCD")
  mm <- intersect(mz, method)
  if(!setequal(mm, method))
  stop("unavailable method(s) requested: ", paste(setdiff(method, mm), collapse = ", "))

# Check for missings in the dataset
  varna <- apply(ouF, 2, function(x) any(is.na(x)))
  if (sum(varna) > 0) {
      cat(colnames(ouF)[varna], "\n")
      stop("These variable(s) have missing values and ",
      "(most) outlier methods cannot deal with them.", call.=FALSE)
    }

# Check if all variables are numeric and warn about integer variables
  n1 <- ncol(ouF)
  varnum <- sapply(ouF, is.numeric)
  if (sum(varnum) < n1) {
    cat(colnames(ouF)[!(varnum)], "\n")
    stop("These variable(s) are not numeric.  (Most) outlier",
    "methods cannot deal with non-numeric variables.", call.=FALSE)
    }
#  varint <- sapply(ouF, is.integer)
#  if (sum(varint) > 0) {
#   message("Some variables are of class integer.  Outlier methods may produce poor results, if there is a lot of heaping.")
#    }

  n2 <- nrow(ouF)
  nw <- sum(choose(n1, k1:K)) #no of variable combinations to be analysed
  mx <- length(mm)
  if (nw*mx > 1000) {
    message("There are ", nw, " possible variable combinations and ", mx, " methods, it could take a while.")
    }

# Check parameters are within allowed limits
  stopifnot(K > 0, K <= n1, n1 > 1)
  stopifnot(k1 > 0, k1 <= K)
  stopifnot(min(tols) > 0, max(tols) < 1)
  stopifnot(min(boxplotLimits) > 0)
  stopifnot(length(tols) <= length(boxplotLimits) & length(tols) <= 3)
# stopifnot(length(tols) <= 3)
# stopifnot(length(unique(tols)) == length(unique(boxplotLimits)))
  ua <- length(tols)

  stopifnot(tolHDo >= 0, tolHDo < 1)
  stopifnot(tolPCS >= 0, tolPCS <= 0.5)
  stopifnot(tolBAC >= 0, tolBAC < 1)
  stopifnot(toladj >= 0, toladj < 1)
  stopifnot(tolDDC >= 0, tolDDC < 1)
  stopifnot(tolMCD >= 0, tolMCD < 1)

#Sort tols and boxplotLimits into descending and ascending order respectively
  tols <- sort(tols, decreasing = TRUE)
  boxplotLimits <- sort(boxplotLimits)

# Allow up to 3 significance levels if only 1 method
  stopifnot((mx > 1 & length(tols) == 1)|(mx == 1 & ua <= 3))
  if (mx > 1) {
    boxplotLimit <- boxplotLimits[1]
      Ax <- c(tolHDo*("HDo" %in% mm), tolPCS*("PCS" %in% mm),
      tolBAC*("BAC" %in% mm), toladj*("adjOut" %in% mm),
      tolDDC*("DDC" %in% mm), tolMCD*("MCD" %in% mm))
      tols <- Ax[Ax>0]
      outList <- mapply(O3a, mm, tols, MoreArgs=list(ouF=ouF, K=K,
      k1=k1, boxplotLimit=boxplotLimit, n1=n1, n2=n2))
      return(list(data = ouF, nw = nw, mm = mm, tols = tols, outList = outList))
      } else {
      outList <- mapply(O3a, tols, boxplotLimits[1:length(tols)], MoreArgs=list(ouF=ouF, K=K,
      k1=k1, mm=mm, n1=n1, n2=n2))
    return(list(data = ouF, nw = nw, mm = mm, tols = tols, outList = outList))
    }
   }
