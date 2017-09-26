# Function to sort the outlier indicator matrix and other plotting preparation
# Order both rows and columns
# Rows are ordered by nos of vars within combinations and then by nos of outliers
# Columns are ordered so that variables come first which appear in most combinations
# Then a gap, then the outliers in order of how often they are found

sortO3 <- function(z1, n1=n1, nz=nz) {
    nc <- nrow(z1)
    n3 <- n1+2
    z1 <- z1 %>% mutate(ID=seq(1:nc), sumV=rowSums(z1[ ,1:n1]),
          sumR=ifelse((nz*nc)>1, rowSums(z1[ ,(n1+2):(nz+n1+1)]), z1[ ,n1+2]))
    z1p <- z1 %>% gather(sN, sB, -ID, -sumV, -sumR)
    z1p <- z1p %>% mutate(pID=factor(ID))
    z1p <- z1p %>% mutate(psN=factor(sN))
    z1p <- z1p %>% group_by(psN) %>% mutate(sumS=sum(sB)) %>% ungroup()
    z1p <- z1p %>% mutate(xsumR=50-sumR)
    z1p <- within(z1p, pID <- reorder(pID, xsumR))
    z1p <- within(z1p, pID <- reorder(pID, sumV))
    z1p <- within(z1p, psNx <- reorder(psN, sumS))
    return(z1p)
}
