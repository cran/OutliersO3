# quiets concerns of R CMD check
if(getRversion() >= "2.15.1")  utils::globalVariables(c("ID",
            "outlierIndices", "pID", "psN", "psNx", "sB", "s1B",
            "sN", "sumR", "sumS", "sumV", "xsumR", "oh", "HDo", "case", "Case", "X", "SsB"))

# Main function--------------

O3plotM <- function(outResults, caseNames=paste0("X", 1:nrow(outResults$data)),
                   sortVars=TRUE, coltxtsize=14, O3control=O3plotColours()) {

ouF <- outResults$data
n1 <- ncol(ouF)
n2 <- nrow(ouF)
nw <- outResults$nw
mm <- outResults$mm
mx <- length(mm)
if(mx < 2) {
stop("You should use O3PlotT for plotting results for a single method", call.=FALSE)
    }
tols <- outResults$tols
outList <- outResults$outList

# There are mx outM's and l2a gives the indices of cases that are ever outliers
  l1a <- vector("list", mx)
  l2a <- NULL
  for  (j in 1:mx) {
    l1a[[j]] <- rlist::list.cases(outList[[j]], outlierIndices)
    l2a <- union(l2a, l1a[[j]])
    }

# Stop if no outliers are found--------------
  nz <- length(l2a)
  if (nz < 1) {
    cat("No outliers found")
    }
  if (nz > 0) {
    #Set up an array with cases ever found to be outliers and one with outlier scores
    zz <- rep(0, nw*nz*mx)
    dim(zz) <- c(nw, nz, mx)
    colnames(zz) <- paste0("X", l2a)
    Cs <- rep(0, mx*nw*n2)
    dim(Cs) <- c(n2, nw, mx)

# Fill a matrix with outlier indices and one with outlier scores for each method
        for (j in 1:mx) {
          l1 <- rlist::list.ungroup(outList[[j]])
          #Prepare a comparison plot if there are only two methods
             for (k in 1:nw) {
               if (mx == 2) {
                 if (length(l1[[3*k-1]]) > 0) {
              zz[k, paste0("X", l1[[3*k-1]]), j] <- 3+j/100
              }
          } else {
          if (mx > 2) {
              if (length(l1[[3*k-1]]) > 0) {zz[k, paste0("X", l1[[3*k-1]]), j] <- 7}
              }
            }
          Cs[, k, j] <- l1[[3*k]] 
          }
        }

        zs <- apply(zz, c(1,2), sum)

# Numbers of outliers found by each method
        zt <- apply(zz, c(2,3), sum)
        zt[zt>0] <- 1
        nOut <- colSums(zt)
        names(nOut) <- mm

# Adjust scores so that outliers found by all methods are red if more than 2
        if (mx > 2) {
          zs[zs>0] <- zs[zs>0] + 7*(6-mx)
        }
      }

#Add variable columns and gap columns to matrix--------------
    zv <- rep(0, nw*(n1 + 2))
    dim(zv) <- c(nw, n1+2)
    colnames(zv) <- c(names(ouF), "gap1", "gap2")
    for (k in 1:nw) {
      zv[k, l1[[3*k-2]]] <- 1/1000
      zv[k, (n1+1):(n1+2)] <- 2/1000
    }
    colnames(zs) <- caseNames[l2a]
    z1 <- data.frame(cbind(zv, zs))

# Remove combinations with no outliers--------------
    if (nz > 1) {
    z1 <- z1[rowSums(z1[ , (n1+3):(n1+2+nz)]) > 0, ]
    } else {
      if (nz == 1)
      z1 <- z1[z1[ , (n1+3)] > 0, ]
    }

# Sort and prepare----------------------------
    z1p <- sortO3(z1, n1=n1, nz=nz, ouF, sortVars)

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
    vc <- c(rep("black", n1), rep("white", 2), rep("black", nz))
    
    ga <- ggplot(z1p, aes(psNx, forcats::fct_rev(pID), group = s1B, fill = s1B)) +
      geom_tile(colour = "grey30") + labs(x = NULL, y = NULL) +
      theme(legend.position = "bottom", panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size=coltxtsize, angle=45, hjust=0, vjust=0, colour=vc),
      axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
      scale_x_discrete(position = "top") +
      geom_hline(yintercept = tw, lty = 3, colour = "blue", size=0.75)
      if (mx == 2) {
         gO3 <- ga  + scale_fill_manual(name = "Outliers identified by", values = outCols,
      breaks=c("3.01", "3.02", "6.03"), labels=c(mm[1], mm[2], "both"))
         } else {
      if (mx == 3) {
      gO3 <- ga + scale_fill_manual(name = "No of methods identifying outliers",
      values = outCols, breaks=c("28", "35", "42"), labels=c("1", "2", "3"))
      } else {
      if (mx == 4) {
      gO3 <- ga + scale_fill_manual(name = "No of methods identifying outliers",
      values = outCols, breaks=c("21", "28", "35", "42"), labels=c("1", "2", "3", "4"))
      } else {
      if (mx == 5) {
      gO3 <- ga + scale_fill_manual(name = "No of methods identifying outliers",
      values = outCols, breaks=c("14", "21", "28", "35", "42"),
      labels=c("1", "2", "3", "4", "5"))
      } else {
           gO3 <- ga  + scale_fill_manual(name = "No of methods identifying outliers",
           values = outCols, breaks=c("7", "14", "21", "28", "35", "42"),
           labels=c("1", "2", "3", "4", "5", "6"))
             }
           }
         }
       }
       
# Draw an O3 plot for 3 or more methods excluding outliers which are only found by 1 method
      if (mx > 2) {
            z1q <- z1p
            z1q[z1q$sB==7*(7-mx), "sB"] <- 0
            z1q <- z1q %>% group_by(psNx) %>% mutate(SsB = sum(sB)) %>% filter(SsB > 0)
            z1q <- z1q %>% mutate(s1B = factor(sB, levels = c("0", "0.001", "0.002", "3",
            "4", "5", "3.01", "3.02", "6.03", "7", "14", "21", "28", "35", "42")))
            names(outCols) <- levels(z1q$s1B)

      gb <- ggplot(z1q, aes(psNx, forcats::fct_rev(pID), group = s1B, fill = s1B)) +
        geom_tile(colour = "grey30") + labs(x = NULL, y = NULL) +
        scale_fill_manual(name = "s1B", values = outCols) +
        theme(plot.title = element_text(size = 18, hjust = 0.5),
        legend.position = "none", panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(size=coltxtsize, angle=45, hjust=0, vjust=0, colour=vc),
        axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
        scale_x_discrete(position = "top") +
        geom_hline(yintercept = tw, lty = 3, colour = "blue")
        gO3x <- gb  + ggtitle(paste("O3 plot of outliers found by at least 2 of",
                mx, "methods"))
                 } else {
                 gO3x <- NULL
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
    gpcp <- gp  + ggtitle(paste("Cases ever found to be outliers by",
        mx, "methods"))

# Create dataset of outliers by case, combination and method--------------
        gg <- data.frame(matrix(zv[,1:n1], nrow = nw))
        if(nw > 1) {
          colnames(gg) <- colnames(zv[,1:n1])
        } else {
          colnames(gg) <- names(zv[,1:n1])
        }
        gg[gg>0] <- 1
        gg <- unite(gg, "Combination", c(1:n1), sep="")
        rownames(zz) <- paste0("c", gg$Combination)
        colnames(zz) <- caseNames[l2a]
        dimnames(zz)[[3]] <- names(nOut)
        zz[zz>0] <- 1
        outs <- memisc::to.data.frame(zz, as.vars = 0, name = "Outlier")
        colnames(outs)[1:3] <- c("Combination", "Name", "Method")
        outs <- outs[outs$Outlier > 0, 1:3]
        Cases <- data.frame(caseNames, Case=1:nrow(outResults$data))
        outy <- merge(outs, Cases, by.x = "Name", by.y = "caseNames")

# Draw a pcp of method distances (scores) for highest combination with outliers highlighted
        if(nw > 1) {
          dimnames(Cs)[[2]] <- rownames(zz)
        } else {
          dimnames(Cs)[[2]] <- list(rownames(zz))
        }
        dimnames(Cs)[[3]] <- mm
        Cx <- data.frame(Cs[ , nw, ])
        Cx$oh <- rep(0, n2)
        Cx[l2a, "oh"] <- "A"
        colnames(Cx) <- c(mm, "oh")
        if("HDo" %in% mm) {
        Cx <- Cx %>% select(-HDo) #HDo provides no distances
        }
        mx <- length(mm)-1*("HDo" %in% mm)
        if(mx > 1) {
        gc <- ggparcoord(Cx %>% arrange(oh), scale = "uniminmax", columns=1:mx,
          groupColumn="oh") + labs(x = NULL, y = NULL) +
          scale_colour_manual(values = c("grey70", "red")) +
          theme(plot.title = element_text(size = 18, hjust = 0.5),
          legend.position = "none", axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
        gMethods <- gc  + ggtitle(paste("Distances for each method with outliers highlighted"))
        } else {
        gMethods <- NULL}
    return(list(nOut = nOut, gpcp = gpcp, gO3 = gO3, gO3x=gO3x, gMethods = gMethods, outsTable = outy, Cs = Cs))
  }
