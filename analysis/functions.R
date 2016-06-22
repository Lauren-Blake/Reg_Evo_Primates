# This is a list of source functions (functions that are used in multiple analyses)

#PCA function (original code from Julien Roux)

plot_scores <- function(pca, scores, n, m, cols, points=F, pchs =20, legend=F){
  xmin <- min(scores[,n]) - (max(scores[,n]) - min(scores[,n]))*0.05
  if (legend == T){ ## let some room (35%) for a legend                                                                                                                                                 
    xmax <- max(scores[,n]) + (max(scores[,n]) - min(scores[,n]))*0.50
  }
  else {
    xmax <- max(scores[,n]) + (max(scores[,n]) - min(scores[,n]))*0.05
  }
  ymin <- min(scores[,m]) - (max(scores[,m]) - min(scores[,m]))*0.05
  ymax <- max(scores[,m]) + (max(scores[,m]) - min(scores[,m]))*0.05
  plot(scores[,n], scores[,m], xlab=paste("PC", n, ": ", round(summary(pca)$importance[2,n],3)*100, "% variance explained", sep=""), ylab=paste("PC", m, ": ", round(summary(pca)$importance[2,m],3)*100, "% variance explained", sep=""), xlim=c(xmin, xmax), ylim=c(ymin, ymax), type="n")
  if (points == F){
    text(scores[,n],scores[,m], rownames(scores), col=cols, cex=1)
  }
  else {
    points(scores[,n],scores[,m], col=cols, pch=pchs, cex=1.3)
  }
}
