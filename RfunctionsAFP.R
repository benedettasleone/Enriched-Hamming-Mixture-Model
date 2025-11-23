CramerV <- function(dati) {
  p <- ncol(dati)
  n <- nrow(dati)
  Vcramer <- matrix(0, p, p)
  
  for (i in 1:p) {
    for (j in 1:p) {
      tab <- table(dati[, i], dati[, j])
      
      expected <- (rowSums(tab) %*% t(colSums(tab))) / n
      
      chi2.stat <- sum((tab - expected)^2 / expected)
      
      h <- nrow(tab)
      k <- ncol(tab)
      chi2.max <- n * min(h - 1, k - 1)
      
      Vcramer[i, j] <- sqrt(chi2.stat / chi2.max)
    }
  }
  return(Vcramer)
}

myplotpsm2 = function (psm, method = "complete", main = NULL, ...)
{
  library(fields)
  if (any(psm != t(psm)) | any(psm > 1) | any(psm < 0) | sum(diag(psm)) !=
      nrow(psm)) {
    stop("psm must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")
  }
  
  colors <- sort(heat.colors(30), decreasing = TRUE)
  n <- nrow(psm)
  
  par(mar=c(0, 0, 1.5, 0), cex.axis=0.6, cex.main=0.9)
  
  psm_scaled <- psm
  min_val <- min(psm[psm != 1])
  max_val <- max(psm[psm != 1])
  breaks <- seq(0, 1, length.out = 31)
  
  heatmap(psm, Rowv=NA, Colv=NA, revC=TRUE, symm = TRUE,
          cexRow=0.6, cexCol=0.6,
          margins=c(7, 2),
          col = colors,
          breaks = breaks,
          main = main,
          keep.dendro=FALSE)
  
  #box()
}
