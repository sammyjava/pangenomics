## plot the results of an Rtsne run
## provide a DF with the colors per path
plot.rtsne = function(rtsne, paths) {
    plot(rtsne$Y, col=paths$color, pch=19, xlab="tSNE dimension 1", ylab="tSNE dimension 2")
    
    ## legend(x="topright", bty="n", pch=19, cex=cex,
    ##        col=c(casecol,ctrlcol),
    ##        c(paste(caseNum,"cases"), paste(ctrlNum,"controls"))
    ##        )
    
    legend(x="topleft", bty="n", cex=0.8, 
           c("Rtsne parameters",
             paste("origD",rtsne$origD,sep="="),
             paste("perplexity",rtsne$perplexity,sep="="),
             paste("theta",round(rtsne$theta,3),sep="="),
             paste("max_iter",rtsne$max_iter,sep="="),
             paste("stop_lying_iter",rtsne$stop_lying_iter,sep="="),
             paste("mom_switch_iter",rtsne$mom_switch_iter,sep="="),
             paste("momentum",rtsne$momentum,sep="="),
             paste("final_momentum",rtsne$final_momentum,sep="="),
             paste("eta",rtsne$eta,sep="="),
             paste("exaggeration_factor",rtsne$exaggeration_factor,sep="=")
             )
           )
}
