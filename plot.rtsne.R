## plot the results of an Rtsne run
plot.rtsne = function(rtsne) {
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
    
    plot(rtsne$Y, col=colors, pch=19, cex=cex, main=prefix, xlab="tSNE dimension 1", ylab="tSNE dimension 2")
    
    legend(x="topright", bty="n", pch=19, cex=cex,
           col=c(casecol,ctrlcol),
           c(paste(caseNum,"cases"), paste(ctrlNum,"controls"))
           )
    
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
