# quiets concerns of R CMD check
if(getRversion() >= "2.15.1")  utils::globalVariables(c("ID",
            "outlierIndices", "pID", "psN", "psNx", "sB", "s1B",
            "sN", "sumR", "sumS", "sumV", "xsumR", "oh"))

# Main function--------------

O3plotT <- function(outResults, caseNames=as.character(1:nrow(outResults$data)),
                   O3control=O3plotColours()) {

mx <- 1
ouF <- outResults$data
n1 <- ncol(ouF)
n2 <- nrow(ouF)
nw <- outResults$nw
mm <- outResults$mm
if(length(mm) > 1) {
stop("You should use O3PlotM for plotting results for more than one method", call.=FALSE)
    }
tols <- outResults$tols
mxm <- length(tols)
outList <- outResults$outList


# There are mx outM's and l2a gives the indices of cases that are ever outliers
  l1a <- vector("list", mxm)
  l2a <- NULL
  for  (j in 1:mxm) {
    l1a[[j]] <- rlist::list.cases(outList[[j]], outlierIndices)
    l2a <- union(l2a, l1a[[j]])
    }

# Stop if no outliers are found--------------
  nz <- length(l2a)
  if (nz < 1) {
    cat("No outliers found for tol equal to ", max(tols))
    }
  if (nz > 0) {
    #Set up an array with cases ever found to be outliers and one with outlier scores
    zz <- rep(0, nw*nz*mxm)
    dim(zz) <- c(nw, nz, mxm)
    colnames(zz) <- paste0("X", l2a)
    Ds <- rep(0, mxm*nw*n2)
    dim(Ds) <- c(n2, nw, mxm)

# Fill a matrix with outlier indices and one with outlier scores for each tol level
      val <- c(3, 4, 5)
        for (j in 1:mxm) {
        l1 <- rlist::list.ungroup(outList[[j]])
        for (k in 1:nw) {
          if (length(l1[[3*k-1]]) > 0) {
            zz[k, paste0("X", l1[[3*k-1]]), j] <- val[3-mxm+j]
            }
          Ds[, k, j] <- l1[[3*k]]
          }
        }
        
#Summarise the outlier indices array
      zs <- apply(zz, c(1,2), max)
      if (mxm == 1) {
        nOut <- nz
        names(nOut) <- mm
        }
      if (mxm > 1) {
        zt <- apply(zz, c(2,3), sum)
        zt[zt>0] <- 1
        nOut <- colSums(zt)
        names(nOut) <- tols
        }


#Add variable columns and gap column to matrix--------------
    zv <- rep(0, nw*(n1 + 1))
    dim(zv) <- c(nw, n1+1)
    colnames(zv) <- c(names(ouF), "._______.")
    for (k in 1:nw) {
      zv[k, l1[[3*k-2]]] <- 1/1000
      zv[k, n1+1] <- 2/1000
    }
    colnames(zs) <- caseNames[l2a]
    z1 <- data.frame(cbind(zv, zs))

# Remove combinations with no outliers--------------
    if (nz > 1) {
    z1 <- z1[rowSums(z1[ , (n1+2):(n1+1+nz)]) > 0, ]
    } else {
      if (nz == 1)
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
    outCols <- c("grey95", "grey60", "white", O3control$colours)

    z1p <- z1p %>% mutate(s1B = factor(sB, levels = c("0", "0.001", "0.002", "3",
      "4", "5", "3.01", "3.02", "6.03", "7", "14", "21", "28", "35", "42")))
    names(outCols) <- levels(z1p$s1B)

    ga <- ggplot(z1p, aes(psNx, forcats::fct_rev(pID), group = s1B, fill = s1B)) +
      geom_tile(colour = "grey30") + labs(x = NULL, y = NULL) +
      scale_fill_manual(name = "s1B", values = outCols) +
      theme(plot.title = element_text(size = 18, hjust = 0.5),
      legend.position = "none", panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size=14, angle=45, hjust=0, vjust=0),
      axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
      scale_x_discrete(position = "top") +
      geom_hline(yintercept = tw, lty = 3, colour = "blue")
    if (mxm == 3) {
      gO3 <- ga  + ggtitle(paste("O3 plot of outliers found by ", mm,
        "\n(",outCols[4], " for ", tols[1], " outliers, ", outCols[5], " for ", tols[2],
        "\n outliers, ", outCols[6], " for ", tols[3], " outliers)", sep=""))
    } else {
      if (mxm == 2) {
      gO3 <- ga  + ggtitle(paste("O3 plot of outliers found by ", mm, "\n(",outCols[5],
       " for ", tols[1], " outliers,\n", outCols[6], " for ",
        tols[2], " outliers)", sep=""))
      } else {
        if (mxm == 1) {
        gO3 <- ga  + ggtitle(paste("O3 plot of outliers found by ", mm,
          "\n(",outCols[6], " for ", tols[1], " outliers)", sep=""))
        }
       }
    }

# Set up PCP with cases ever identified as outliers coloured red--------------
    ouF1 <- ouF
# Create highlighting variable
    ouF1$oh <- rep(0, nrow(ouF))
# Outliers are listed in l2a
    ouF1[l2a, "oh"] <- "A"
    ouF1 <- ouF1 %>% mutate(alev = ifelse(oh == 0, 0.5, 1))
    gp <- ggparcoord(ouF1 %>% arrange(oh), scale = "uniminmax", columns=1:n1,
      groupColumn="oh", alphaLines="alev") + labs(x = NULL, y = NULL) +
      scale_colour_manual(values = c("grey70", "red")) +
      theme(plot.title = element_text(size = 18, hjust = 0.5),
      legend.position = "none", axis.ticks.y = element_blank(),
      axis.text.y = element_blank())
    gpcp <- gp + ggtitle(paste("Cases ever found to be outliers \n by method",
        mm, "\n for tol=", tols[1]))


# Create dataset of outliers by case, combination and tol--------------
        gg <- data.frame(matrix(zv[,1:n1], nrow = nw))
        if(nw > 1) {
          colnames(gg) <- colnames(zv[,1:n1])
        } else {
          colnames(gg) <- names(zv[,1:n1])
        }
        gg[gg>0] <- 1
        gg <- unite(gg, "Combination", c(1:n1), sep="")
        rownames(zz) <- paste0("c", gg$Combination)
        dimnames(zz)[[3]] <- names(nOut)
        zz[zz>0] <- 1
        outs <- memisc::to.data.frame(zz, as.vars = 0, name = "Outlier")
        colnames(outs)[1:3] <- c("Combination", "Case", "tol")
        outs <- outs[outs$Outlier > 0, 1:3]

# Draw a pcp of method distances (scores) with lowest tol outliers highlighted unless method is HDo
        dimnames(Ds)[[2]] <- rownames(zz)
        dimnames(Ds)[[3]] <- paste0("T", tols)
        if(mm=="HDo") {
           return(list(nOut = nOut, gpcp = gpcp, gO3 = gO3, outsTable = outs, Ds = Ds)) 
        } else {
        Dx <- data.frame(Ds[ , , mxm])
        Dx$oh <- rep(0, n2)
        Dx[l1a[[mxm]], "oh"] <- "A"
        colnames(Dx) <- c(rownames(zz), "oh")
        gd <- ggparcoord(Dx %>% arrange(oh), scale = "uniminmax", columns=1:nw,
          groupColumn="oh") + labs(x = NULL, y = NULL) +
          scale_colour_manual(values = c("grey70", "red")) +
          theme(plot.title = element_text(size = 18, hjust = 0.5),
          legend.position = "none", axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
        gCombs <- gd  + ggtitle(paste("Scores for each combination using\n method", mm, "with outliers highlighted"))

    return(list(nOut = nOut, gpcp = gpcp, gO3 = gO3, gCombs = gCombs, outsTable = outs, Ds = Ds))
    }
  }
}
