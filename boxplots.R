
boxplot_resection = function(sel.df, tp1, tp2, sdfx, sdfid="Subject_ID", 
                             selected, colpalette=c("#0072B2", "#E69F00", "#1B9E77", "#D95F02"),
                             tps = c("TP1", "TP2"), stratify=c("macroradical_resection", "cluster_nam"), 
                             save=TRUE, filename="default", filedim = c(7, 6),
                             cond=list(c("TP2", 1, "PPG1"), c("TP1", 1, "PPG1"),
                                       c("TP2", 1, "PPG2"), c("TP1", 1, "PPG2"),
                                       c("TP2", 0, "PPG1"), c("TP1", 0, "PPG1"),
                                       c("TP2", 0, "PPG2"), c("TP1", 0, "PPG2"),
                                       c("Δ", 1, "PPG1"), c("Δ", 1, "PPG2"),
                                       c("Δ", 0, "PPG1"), c("Δ", 0, "PPG2")),
                             strata=c("TP2_CRET_PPG1", "TP1_CRET_PPG1",
                                      "TP2_CRET_PPG2", "TP1_CRET_PPG2",
                                      "TP2_Residual_PPG1", "TP1_Residual_PPG1",
                                      "TP2_Residual_PPG2", "TP1_Residual_PPG2",
                                      "Δ_CRET_PPG1", "Δ_CRET_PPG2",
                                      "Δ_Residual_PPG1", "Δ_Residual_PPG2"),
                             ymin=NULL, ymax=NULL, ycut=NULL) {
  tp1 = tp1[which(names(tp1) %in% names(tp2))]
  tp2 = tp2[which(names(tp2) %in% names(tp1))]
  tp2 = tp2[names(tp1)]
  print(paste("All names of TP1 correspond to TP2:", all(names(tp1)==names(tp2))))
  
  sdfx = sdf[which(sdfx[,sdfid] %in% names(tp1)),]
  sdfx1 = sdfx[match(rep(names(tp1), 2), sdfx[,sdfid]),]
  sdfx2 = sdfx[match(names(tp1), sdfx[,sdfid]),]
  
  for (i in selected) {
    xdf1 = cbind("Prot_quant"=sel.df[i,c(tp1,tp2)], 
                 "TP"=c(rep(tps[1], length(tp1)), rep(tps[2], length(tp2))), sdfx1)
    xdf2 = cbind("Prot_quant"=sel.df[i,tp2]-sel.df[i,tp1], 
                 "TP"=rep("Δ", length(tp2)), sdfx2)
    xdf = rbind(xdf1, xdf2)
    xdf = xdf |> arrange("TP") |> arrange(!!sym(stratify[1])) |> 
      arrange(!!sym(stratify[2]))
    
    labels = rep(NA, length(xdf$Prot_quant))
    for(j in 1:length(cond)) {
      labels[which(xdf$TP==cond[[j]][1] & 
                     xdf[,stratify[1]]==cond[[j]][2] & 
                     xdf[,stratify[2]]==cond[[j]][3])] = strata[[j]]
    }
    table(labels)
    
    xdf$Category = as.factor(labels)
    xdf$Groups = as.character(str_extract(xdf$Category, "CR.+|R.+"))
    
    if(is.null(ymin)) {ymin = floor(min(xdf$Prot_quant, na.rm=T))}
    if(is.null(ymax)) {ymax = ceiling(max(xdf$Prot_quant, na.rm=T))}
    if(is.null(ycut)) {ycut = 10}
    
    p = ggboxplot(xdf, x="TP", y="Prot_quant", color="Groups",
                  add="jitter", facet.by="Groups")  + ggtitle(i) +
      xlab("Time point") + ylab("Plasma levels") +
      stat_compare_means(data=xdf[-which(xdf$TP=="Δ"),], aes(group=TP), 
                         method="t.test", paired=T, label.y=ymax/2) +
      scale_y_continuous(breaks=round(seq(ymin, ymax, length.out = ycut), 2)) +
      scale_color_manual("Cancer \nresection", values=colpalette)
    p
    
    if(save==TRUE) {
      if(filename=="default") {
        quartz(type="pdf", file=paste0("boxplots/", paste(stratify, collapse="_"), "_", i, "_levels_", ".pdf"), 
               dpi=600, width=filedim[1], height=filedim[2])
        print(p)
        dev.off()
      } else {
        quartz(type="pdf", file=paste0("boxplots/", filename, i, "_levels_", ".pdf"), 
               dpi=600, width=filedim[1], height=filedim[2])
        print(p)
        dev.off()
      }
    }
    
    print(paste0("Protein: ", i, " done; run: ", match(i, selected)))
    
  }
}
