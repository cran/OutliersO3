# quiets concerns of R CMD check
if(getRversion() >= "2.15.1")  utils::globalVariables(c("ID",
            "outlierIndices", "pID", "psN", "psNx", "sB", "s1B",
            "sN", "sumR", "sumS", "sumV", "xsumR", "oh"))

# Main function--------------

O3plot <- function(ouF, k1=1, K=ncol(ouF), mm="HDo", Alphas=0.05, Coefs=6,
     outColours = c("khaki", "yellow", "red", "salmon", "chocolate", "red",
       "slategray1", "slategray2", "slategray3", "slategray4", "red")) {
           
  ouF <- data.frame(ouF) #in case the dataset is a tibble
  
# Check for missings in the dataset
  varna <- apply(ouF, 2, function(x) any(is.na(x)))
  if (sum(varna) > 0) {
      cat(colnames(ouF)[varna])
      stop("These variable(s) have missing values and ",
      "(most) outlier methods cannot deal with them.", call.=FALSE)
    }

# Check if all variables are numeric and warn about integer variables
  n1 <- ncol(ouF)
  varnum <- sapply(ouF, is.numeric)
  if (sum(varnum) < n1) {
    cat(colnames(ouF)[varnum])
    stop("These variable(s) are not numeric.  (Most) outlier",
    "methods cannot deal with non-numeric variables.", call.=FALSE)
    }
  varint <- sapply(ouF, is.integer)
  if (sum(varint) > 0) {
    message("Some variables are of class integer.  Outlier methods may produce poor results, if there is a lot of heaping.")
    }  

  n2 <- nrow(ouF)
  nw <- sum(choose(n1, k1:K)) #no of variable combinations analysed
  if (nw > 1000) {
    message("There are ", nw, " possible variable combinations, it could take a while.")
    }
  stopifnot(K > 0, K <= n1, n1 > 1)
  stopifnot(k1 > 0, k1 <= K)
  mx <- length(mm)
  stopifnot(sum((mm[1:mx] %in% c("HDo", "PCS", "BAC", "adjOut", "DDC")))==mx)
  stopifnot(min(Alphas) > 0, max(Alphas) < 1)
  stopifnot(min(Coefs) > 0)
  ua <- length(unique(Alphas))
  uc <- length(unique(Coefs))
  stopifnot(length(Alphas)==length(Coefs))
  
# Allow up to 3 significance levels if only 1 method
  stopifnot((mx>1 & length(Alphas)==1)|(mx==1 & ua<=3))
  mxm <- mx
  if (mx > 1) {
    Coef <- Coefs[1]
    Alpha <- Alphas[1]
    outList <- lapply(mm, O3a, ouF=ouF, K=K, k1=k1, Coef=Coef,
        Alpha=Alpha, n1=n1, n2=n2)

# There are mx outM's and l2a gives the indices of cases that are ever outliers
    l1a <- vector("list", mx)
    l2a <- NULL
    for  (j in 1:mx) {
      l1a[[j]] <- rlist::list.cases(outList[[j]]$outM, outlierIndices)
      l2a <- union(l2a, l1a[[j]])
      }
    } else {
  mxm <- length(Alphas)
  outList <- mapply(O3a, Alphas, Coefs, MoreArgs=list(ouF=ouF, K=K,
    k1=k1, mm=mm, n1=n1, n2=n2))
  l1a <- vector("list", mxm)
  l2a <- NULL
  for  (j in 1:mxm) {
    l1a[[j]] <- rlist::list.cases(outList[[j]], outlierIndices)
    l2a <- union(l2a, l1a[[j]])
    }
  }

# Stop if no outliers are found--------------
  nz <- length(l2a)
  if (nz < 1) {
    cat("No outliers found for alpha equal to ", max(Alphas))
    }
  if (nz > 0) {
    #Set up an array with the cases ever found to be outliers
    zz <- rep(0, nw*nz*mxm)
    dim(zz) <- c(nw, nz, mxm)
    nam <- paste0("X", l2a)
    colnames(zz) <- nam

# Fill a matrix for each method/alpha level and summarise--------------
    if (mx==1) {
      val <- c(3, 4, 5)
        for (j in 1:mxm) {
        l1 <- rlist::list.ungroup(outList[[j]])
        for (k in 1:nw) {
          if (length(l1[[2*k]])>0) {
            zz[k, paste0("X", l1[[2*k]]), j] <- val[3-mxm+j]
            }
          }
        }
      zs <- apply(zz, c(1,2), max)
      if (mxm==1) {
        nrout <- nz
        names(nrout) <- mm
        }
      if (mxm > 1) {
        zt <- apply(zz, c(2,3), sum)
        zt[zt>0] <- 1
        nrout <- colSums(zt) 
        names(nrout) <- Alphas
        }
      } else {
        for (j in 1:mx) {
          l1 <- rlist::list.ungroup(outList[[j]]$outM)
          #Prepare a comparison plot if there are only two methods
          if (mx==2) {
            for (k in 1:nw) {
            if (length(l1[[2*k]])>0) {
              zz[k, paste0("X", l1[[2*k]]), j] <- 3+j/100
              }
            }
          } else {
          if (mx>2) {
            for (k in 1:nw) {
              if (length(l1[[2*k]])>0) {zz[k, paste0("X", l1[[2*k]]), j] <- 7}
              }
            }
          }
        }

        zs <- apply(zz, c(1,2), sum)
        
# Numbers of outliers found by each method
        zt <- apply(zz, c(2,3), sum)
        zt[zt>0] <- 1
        nrout <- colSums(zt)
        names(nrout) <- mm

# Adjust scores so that outliers found by all methods are red if more than 2
        if (mx>2) {
          zs[zs>0] <- zs[zs>0] + 7*(5-mx)
        }
      }

#Add variable columns and gap column to matrix--------------
    zv <- rep(0, nw*(n1 + 1))
    dim(zv) <- c(nw, n1+1)
    colnames(zv) <- c(names(ouF), "._______.")
    for (k in 1:nw) {
      zv[k, l1[[2*k-1]]] <- 1/1000
      zv[k, n1+1] <- 2/1000
    }
    z1 <- data.frame(cbind(zv, zs))

# Remove combinations with no outliers--------------
    if (nz > 1) {
    z1 <- z1[rowSums(z1[ , (n1+2):(n1+1+nz)]) > 0, ]
    } else {
      if (nz==1)
      z1 <- z1[z1[ , (n1+2)] > 0, ]
    }

# Sort and prepare----------------------------
    z1p <- sortO3(z1, n1=n1, nz=nz)
    
# Plot------------------------------------------

# Prepare blue lines separating numbers of combinations--------------
    nc <- nrow(z1)
    tx <- rep(0,n1)
    t1d <- data.frame(table((nc+1)*rowSums(z1[ ,1:n1])))
    tx[t1d$Var1] <- t1d$Freq
    ty <- cumsum(rev(tx)) + 0.5
    tw <- ty[ty > 0.5 & ty < nc + 0.5]

# Set the colours for the plot------------------
# First add the background colours
    outCols <- c("grey95", "grey60", "white", outColours)
 
    z1p <- z1p %>% mutate(s1B=factor(sB, levels=c("0", "0.001", "0.002", "3",
      "4", "5", "3.01", "3.02", "6.03", "7", "14", "21", "28", "35")))
    names(outCols) <- levels(z1p$s1B)

    ga <- ggplot(z1p, aes(psNx, forcats::fct_rev(pID), group=s1B, fill=s1B)) +
      geom_tile(colour="grey30") + labs(x=NULL, y=NULL) +
      scale_fill_manual(name = "s1B", values = outCols) +
      theme(plot.title=element_text(size = 18, hjust = 0.5),
      legend.position="none", panel.background = element_rect(fill = "white"),
      axis.text.x=element_text(size=14, angle=45, hjust=0, vjust=0),
      axis.ticks.y=element_blank(), axis.text.y=element_blank()) +
      scale_x_discrete(position = "top") +
      geom_hline(yintercept=tw, lty=3, colour="blue")
    if (mx==1 && ua==3) {
      gO3 <- ga  + ggtitle(paste("O3 plot of outliers found by", mm,
        "\n(",outCols[4], " for", Alphas[1], "outliers,", outCols[5], "for", Alphas[2],
        "\n outliers,", outCols[6], "for", Alphas[3], "outliers)"))
    } else {
      if (mx==1 && ua==2) {
      gO3 <- ga  + ggtitle(paste("O3 plot of outliers found by", mm,
      "\n(",outCols[5], " for", Alphas[1], "outliers,\n", outCols[6], "for", Alphas[2], "outliers)"))
      } else {
        if (mx==1 && ua==1) {
        gO3 <- ga  + ggtitle(paste("O3 plot of outliers found by", mm,
          "\n(",outCols[6], "for", Alphas[1], "outliers)"))
        } else {
         if (mx==2) {
         gO3 <- ga  + ggtitle(paste("O3 plot of outliers found by methods \n",
         mm[1], "(",outCols[7], ") and", mm[2], "(",outCols[8],") for alpha=", Alpha,
           "\n(outliers found by both are in ", outCols[9],")"))
         } else {
           gO3 <- ga  + ggtitle(paste("O3 plot of outliers found by",
             mx, "methods for alpha=", Alpha)) }
         }
       }
    }

# Set up PCP with cases ever identified as outliers coloured red--------------
    ouF1 <- ouF
# Create highlighting variable
    ouF1$oh <- rep(0, nrow(ouF))
# Outliers are listed in l2a
    ouF1[l2a, "oh"] <- "A"
    ouF1 <- ouF1 %>% mutate(alev=ifelse(oh==0, 0.5, 1))
    gp <- ggparcoord(ouF1 %>% arrange(oh), scale="uniminmax", columns=1:n1,
      groupColumn="oh", alphaLines="alev") + labs(x=NULL, y=NULL) + 
      scale_colour_manual(values = c("grey70", "red")) +
      theme(plot.title=element_text(size = 18, hjust = 0.5),
      legend.position="none", axis.ticks.y = element_blank(),
      axis.text.y = element_blank())
    if(mx==1) {
      gpcp <- gp + ggtitle(paste("Cases ever found to be outliers \n by method",
        mm, "\n for alpha=", Alphas[1]))
    } else {
      gpcp <- gp  + ggtitle(paste("Cases ever found to be outliers by",
        mx, "methods \n for alpha=", Alpha))
    }
    return(list(nrout=nrout, gpcp=gpcp, gO3=gO3))
  }
}