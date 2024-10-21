##### Nice plots for consensus clustering #####
cclustering_plots = function(results, nk) {
  ##### Create your own nice heatmaps #####
  palsel = c("white", "#6A51A3")
  palsel1 = c(brewer.pal(8, "Dark2")[-3], brewer.pal(3, "Paired"))
  names(palsel1) = as.character(1:10)
  
  for (z in 2:nk) {
    ca = columnAnnotation(Cluster = as.factor(results[[z]]$consensusClass),
                          col=list(Cluster=palsel1[1:z]))
    hm =  Heatmap(results[[z]]$consensusMatrix,
                  name = "Consensus \nvalue",
                  top_annotation = ca,
                  col=palsel, show_row_dend = F)
    
    pdf(paste("heatmap_consensus_k", z, ".pdf", sep=""), height=4, width=5)
    print(hm)
    dev.off()
  }
  
  ##### Create your own nice CDF and Area under CDF #####
  acdf = list()
  dfv = list()
  
  for(j in 2:nk) {
    vals = as.numeric(results[[j]]$consensusMatrix)
    
    vals = sort(vals)
    cvals = list()
    auc = list()
    for(i in 1:length(vals)) {
      cvals[[i]] = length(which(vals<=vals[[i]]))/length(vals)
      if(i>1) {
        auc[[i]] = (vals[[i]]-vals[[i-1]])*cvals[[i]]
      }
    }
    
    acdf[[j]] = sum(unlist(auc))
    dfv[[j]] = as.data.frame(cbind(vals, unlist(cvals)))
    dfv[[j]]$k = j
    
  }
  
  acdf[[3]]-acdf[[2]]
  
  calc = list()
  for(j in 2:10) {
    if(j==2) {
      calc[[j]] = acdf[[j]]
    } else {
      calc[[j]]  = ((acdf[[j]])-acdf[[j-1]])/acdf[[j]]
    }
  }
  
  dfvsum = do.call(rbind, dfv)
  cdfplot = ggplot() +
    geom_line(data=dfvsum, aes(vals, V2, col=as.factor(k)), linewidth=1) +
    theme_classic() +
    scale_y_continuous("CDF", limits=c(0,1), breaks=seq(0,1,0.1)) +
    scale_x_continuous("Consensus index value", breaks=seq(0,1,0.1)) +
    scale_color_manual("k", values=palsel1[-1])
  cdfplot
  
  audf = as.data.frame(cbind(2:10, unlist(calc), unlist(acdf)))
  
  aucplot = ggplot(audf, aes(V1, V2, col = as.factor(V1), fill=as.factor(V1))) +
    geom_line(data=audf, aes(x=V1, y=V2), inherit.aes = F, linewidth=1, col="darkgrey") +
    geom_point(shape=23, size=3) +
    scale_color_manual(guide="none", values=palsel1[-1]) +
    scale_fill_manual("k", values=palsel1[-1]) + 
    scale_x_continuous("k", labels=2:10, breaks=2:10) +
    scale_y_continuous("D (AUCDF)", breaks=seq(0,0.4, 0.1)) +
    theme_classic()
  aucplot
  
  comboplot = ggarrange(cdfplot, aucplot, labels=c("a", "b"), ncol=1, heights = c(0.7,0.3))
  
  ggexport(comboplot, filename="consensus_clustering_CDF.pdf", width=7.25, height=6.5)
  ggexport(comboplot, filename="consensus_clustering_CDF.png", 
           dpi=600, width=7.25, height=6.5)
  return(comboplot)
}