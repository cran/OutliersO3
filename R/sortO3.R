# Function to sort the outlier indicator matrix and other plotting preparation
# Order both rows and columns
# Rows are ordered by nos of vars within combinations and then by nos of outliers
# Columns are ordered so that variables come first which appear in fewest combinations
# Then a gap, then the outliers in order of how often they are found

sortO3 <- function(z1, n1=n1, nz=nz, ouF, sortVars) {
    nc <- nrow(z1)
    z1 <- z1 %>% mutate(ID=seq(1:nc), sumV=rowSums(z1[ ,1:n1]))
    if (nz > 1) {
       z1$sumR <- rowSums(z1[ ,(n1+3):(nz+n1+2)])
       } else {
       z1$sumR <- z1[ ,n1+3]
       }
    z1p <- z1 %>% gather(sN, sB, -ID, -sumV, -sumR)
    z1p <- z1p %>% mutate(pID=factor(ID))
    z1p <- z1p %>% mutate(psN=factor(sN))
    z1p <- z1p %>% group_by(psN) %>% mutate(sumS=sum(sB)) %>% ungroup()
    z1p <- z1p %>% mutate(xsumR=50-sumR)
    z1p <- within(z1p, pID <- reorder(pID, xsumR))
    z1p <- within(z1p, pID <- reorder(pID, sumV))
    z1p <- within(z1p, psNx <- reorder(psN, sumS))
    if (sortVars==FALSE) {
       z1p$psNx <- forcats::fct_relevel(z1p$psNx, colnames(ouF)[1:n1])
       }   
    return(z1p)
}
