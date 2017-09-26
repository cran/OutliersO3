## ----fig.width=7, fig.height=8, fig.align='center'-----------------------
library(OutliersO3)
data(Election2005, package="mbgraphic")
ouF <- Election2005[, c(6, 10, 17, 28)]
O3s <- O3plot(ouF, mm="HDo", Alphas=0.05, Coefs=6)
library(gridExtra)
grid.arrange(O3s$gO3, O3s$gpcp, ncol=1)

## ----fig.width=7, fig.height=8, fig.align='center'-----------------------
O3s2 <- O3plot(ouF, mm="HDo", Alphas=c(0.1, 0.05, 0.01), Coefs=c(3, 6, 10))
grid.arrange(O3s2$gO3, O3s2$gpcp, ncol=1)

## ----fig.width=7, fig.height=8, fig.align='center'-----------------------
O3s3 <- O3plot(ouF, mm=c("HDo", "PCS"), Alphas=0.05, Coefs=6)
grid.arrange(O3s3$gO3, O3s3$gpcp, ncol=1)

## ----message=FALSE-------------------------------------------------------
O3s5 <- O3plot(ouF, mm=c("HDo", "PCS", "BAC", "adjOut", "DDC"), Alphas=0.05, Coefs=6)
cx <- data.frame(outlier_method=names(O3s5$nrout), number_of_outliers=O3s5$nrout)
knitr::kable(cx, row.names=FALSE)

## ----fig.width=7, fig.height=8, fig.align='center'-----------------------
grid.arrange(O3s5$gO3, O3s5$gpcp, ncol=1)

