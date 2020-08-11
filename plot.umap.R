## plot the results of a umap run
plot.umap = function(umap) {
    casecol = rgb(1,0,0, maxColorValue=1, alpha = 0.15, names = "red")
    ctrlcol = rgb(0,0,1, maxColorValue=1, alpha = 0.15, names = "blue")
    
    colors = array(dim=(caseNum+ctrlNum))
    colors[cases] = casecol
    colors[ctrls] = ctrlcol
    
    if ((caseNum+ctrlNum)<1000) {
        cex = 2.0
    } else {
        cex = 1.0
    }
    
    plot(umap$layout[,1], umap$layout[,2], col=colors, pch=19, cex=cex, main=prefix, xlab="umap dimension 1", ylab="umap dimension 2")
    
    legend(x="topright", bty="n", pch=19, cex=cex,
           col=c(casecol,ctrlcol),
           legend=c(paste(caseNum,"cases"), paste(ctrlNum,"controls"))
           )
    
    legend(x="topleft", bty="n", cex=0.8,
           legend=c("umap parameters:",
                    paste("method",umap$config$method,sep="="),
                    paste("n_neighbors",umap$config$n_neighbors,sep="="),
                    paste("n_components",umap$config$n_components,sep="="),
                    paste("metric",umap$config$metric,sep="="),
                    paste("n_epochs",umap$config$n_epochs,sep="="),
                    paste("init",umap$config$init,sep="="),
                    paste("min_dist",umap$config$min_dist,sep="="),
                    paste("alpha",umap$config$alpha,sep="="),
                    paste("gamma",umap$config$gamma,sep="="),
                    paste("negative_sample_rate",umap$config$negative_sample_rate,sep="="),
                    paste("a",round(umap$config$a,3),sep="="),
                    paste("b",round(umap$config$b,3),sep="=")
                    )
           )
}
