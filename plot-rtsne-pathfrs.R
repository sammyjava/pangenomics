## plot the results of an Rtsne run

casecol = rgb(1,0,0, maxColorValue=1, alpha = 0.5, names = "red")
ctrlcol = rgb(0,0,1, maxColorValue=1, alpha = 0.5, names = "blue")

colors = array(dim=(caseNum+ctrlNum))
colors[cases] = casecol
colors[ctrls] = ctrlcol

if ((caseNum+ctrlNum)<1000) {
    cex = 2.0
} else if ((caseNum+ctrlNum)<2000) {
    cex = 1.0
} else {
    cex = 0.5
}

plot(rtsne.pathfrs$Y, col=colors, pch=19, cex=cex, main=prefix, xlab="Rtsne dimension 1", ylab="Rtsne dimension 2")

legend(x="topright", bty="n", pch=19, cex=cex,
       col=c(casecol,ctrlcol),
       c(paste(caseNum,"cases"), paste(ctrlNum,"controls"))
       )
       
legend(x="topleft", bty="n", cex=0.8, 
       c(paste("origD",rtsne.pathfrs$origD,sep="="),
         paste("perplexity",rtsne.pathfrs$perplexity,sep="="),
         paste("theta",round(rtsne.pathfrs$theta,2),sep="="),
         paste("max_iter",rtsne.pathfrs$max_iter,sep="="),
         paste("stop_lying_iter",rtsne.pathfrs$stop_lying_iter,sep="="),
         paste("mom_switch_iter",rtsne.pathfrs$mom_switch_iter,sep="="),
         paste("momentum",rtsne.pathfrs$momentum,sep="="),
         paste("final_momentum",rtsne.pathfrs$final_momentum,sep="="),
         paste("eta",rtsne.pathfrs$eta,sep="="),
         paste("exaggeration_factor",rtsne.pathfrs$exaggeration_factor,sep="=")
         )
       )
