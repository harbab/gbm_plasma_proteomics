test_confounders = function(sdf, selnum, selcat, variable="cluster_nam", group=c("PPG1", "PPG2"), 
                            save=FALSE, plot=FALSE) {
  t = list()
  vals = list()
  sumry = list()
  testpval = list()
  for (i in selnum) {
    vals[[1]] = sdf[,i]
    vals[[2]] = sdf[which(sdf[,variable]==group[1]),i]
    vals[[3]] = sdf[which(sdf[,variable]==group[2]),i]
    
    for(j in 1:3) {
    sumry[[i]][[j]] = list(paste0(round(mean(vals[[j]], na.rm=T),2), " (", 
                  round(median(vals[[j]], na.rm=T), 2), ")", 
                  " ± ", round(sd(vals[[j]], na.rm=T),2)), 
           paste0("[", 
                  round(min(vals[[j]], na.rm=T),2), "-", round(max(vals[[j]], na.rm=T), 2), "], ",
                  length(which(is.na(vals[[j]]))), " NAs"))
    }
    
    t[[i]] = t.test(vals[[3]], vals[[2]])
    testpval[[i]] = t[[i]]$p.value
    
    print(paste("Variable", i, "done."))
  }
  
  numres = as.data.frame(cbind(selnum, do.call(rbind, lapply(sumry, \(x) unlist(x[[1]]))), 
                               do.call(rbind, lapply(sumry, \(x) unlist(x[[2]]))),
                               do.call(rbind, lapply(sumry, \(x) unlist(x[[3]])))))
  numres$t_pval = unlist(testpval)
  colnames(numres) = c("Variable name", 
                       "All, Mean (median) ± SD", "All, [Range], n of NAs",
                       paste(group[1], c(", Mean (median) ± SD", ", [Range], n of NAs")),
                       paste(group[2], c(", Mean (median) ± SD", ", [Range], n of NAs")), "t-test p-value")
  #numres = apply(numres, 2, as.character)
  
  sumry = list()
  vals = list()
  tb = list()
  f = list()
  
  for (i in selcat) {
    vals[[1]] = sdf[,i]
    vals[[2]] = sdf[which(sdf[,variable]==group[1]),i]
    vals[[3]] = sdf[which(sdf[,variable]==group[2]),i]
    
    for(j in 1:3) {
      t = table(vals[[j]])
      sumry[[i]][[j]] = list(paste(names(t), collapse=" : "),
                              paste0(paste(round(t/sum(t), 2)*100, collapse=" : "), ", ",
                                     length(which(is.na(vals[[j]]))), " NAs"))
    }
    
    tb[[i]] = table(sdf[,variable], sdf[,i])
    f[[i]] = fisher.test(tb[[i]])
    testpval[[i]] = f[[i]]$p.value
    print(paste("Variable", i, "done."))
  }
  
  catres = as.data.frame(cbind(selcat, do.call(rbind, lapply(sumry, \(x) unlist(x[[1]]))), 
                               do.call(rbind, lapply(sumry, \(x) unlist(x[[2]]))), 
                               do.call(rbind, lapply(sumry, \(x) unlist(x[[3]])))))
  catres$f_pval = unlist(lapply(f, function(z) z$p.value))
  colnames(catres) = c("Variable name", 
                       "All, Categories", "All, % n of NAs",
                       paste(group[1], c("PPG1, Categories", "PPG1, % n of NAs")),
                       paste(group[2], c("PPG1, Categories", "PPG1, % n of NAs")),
                       "Fisher's test p-value")
  
  if(save==TRUE) {
    write.csv(numres, paste0("table_numeric_confounders_", paste0(group, collapse="_"), ".csv"))
    write.csv(catres, paste0("table_categorical_confounders_", paste0(group, collapse="_"), ".csv"))
  }
  
  p1 = list()
  p2 = list()
  
  if(plot==TRUE) {
    pldf = sdf[-which(is.na(sdf[,variable])),]
    pldf[,variable] = as.factor(pldf[,variable])
    
    for (i in selnum) {
      p1[[i]] = ggplot(pldf, aes_string(x=sym(variable), y=pldf[,i], col=sym(variable))) +
        geom_jitter(size=2) +
        theme_classic() +
        scale_x_discrete("Treatment-naïve plasma protein levels", 
                         labels=c("PPG1", "PPG2")) +
        scale_y_continuous(i) +
        scale_color_manual(values=c("#009E73", "#D55E00")) +
        annotate("text", x=1.5, y=median(pldf[,i], na.rm=T), 
                 label=paste("t-test:", round(testpval[[i]],3)))
      p1[[i]]
      
      ggsave(paste0(i, "_distribution_proteomic_clusters.pdf"), p1[[i]], dpi=600, width=4, height=4)
    }
    
    for (i in selcat) {
      p2[[i]] = ggplot(pldf, aes_string(x=sym(variable), y=as.factor(pldf[,i]), col=sym(variable))) +
        geom_jitter(size=2) +
        theme_classic() +
        geom_hline(yintercept = 1.5, linetype="dashed") +
        scale_x_discrete("Treatment-naïve plasma protein levels", 
                         labels=c("PPG1", "PPG2")) +
        scale_y_discrete(i) +
        scale_color_manual(values=c("#009E73", "#D55E00")) +
        annotate("text", x=1.5, y=1.7, label=paste("Fisher's test:", round(testpval[[i]],3)))
      p2[[i]]
      
      ggsave(paste0(i, "_distribution_proteomic_clusters.pdf"), p2[[i]], dpi=600, width=4, height=4)
    }
  }
  
  results = list(numres, catres, p1, p2)
  names(results) = c("numeric_summary", "categoric_summary",
                     "numeric_plots", "categoric_plots")
  return(results)
}
