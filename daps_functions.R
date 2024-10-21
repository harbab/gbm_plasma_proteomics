##### Differential analysis of proteins #####
##### Function for analysing DAPs with a t test #####
daps_function = function (group1, group2, matrix, protinfo, comparison, paired=FALSE, adj.pval="fdr",
                          save=FALSE, filename="default") {
  t.result <- apply(matrix, 1, function (i) t.test(i[group1],i[group2], paired=paired))
  results = lapply(t.result, unlist)
  results = as.data.frame(do.call(rbind, results))
  
  results$Gene.Name = rownames(results)
  results$adj.p_values <- p.adjust(results$p.value, method = adj.pval)
  results = results |> 
    left_join(protinfo |> select(Gene.Name, Gene.ID, Protein.ID.s., Description), by="Gene.Name")
  
  if(paired==FALSE) {
    results[,c(1:9)] = apply(results[,c(1:9)], 2, as.numeric)
    colnames(results)[grep("estimate|name", colnames(results))] = c("Mean_Group1", "Mean_Group2", "Comparison")
    results$log2_FC = results$Mean_Group1-results$Mean_Group2
    results = results[,c(13,18,3,14,17,15,16,6,7,1,2,4,5,12,8:11)]
  } else {
    results[,c(1:8)] = apply(results[,c(1:8)], 2, as.numeric)
    colnames(results)[grep("estimate|name", colnames(results))] = c("log2_FC", "Comparison")
    results = results[,c(12,6,3,13,16,14,15,1:2,4,5,11,7:10)]
  }
  
  results$Comparison = comparison
  results <- results[order(results$log2_FC, decreasing = T),]
  
  if(save==TRUE) {
    if(filename=="default") {
      write.csv(results, paste0("t_test_", comparison, "_all_proteins.csv"))
    } else {
      write.csv(results, paste0("t_test_", filename, "_all_proteins.csv"))
    }
  }
  
  return(results)
}

##### Volcano plot for t test ######
ggplot_func = function(ttres, label_up, label_down, colpalette, colpal, overlaps, xnam, ynam, 
                       xmin, xmax, ymax, legend.pos, axis.breaks, other_params, ...) {
  cov = ggplot(ttres, aes(x=log2_FC, y=-log10(adj.p_values), col=alteration)) + 
    geom_point(aes(alpha=as.numeric(shadcat))) + 
    theme_classic() + 
    scale_alpha(range=c(0.1,1), guide="none") +
    geom_hline(yintercept=-log10(0.05), col="black", linetype="dotdash", linewidth=1) + 
    scale_color_manual("Alteration", values=colpalette) +
    theme(legend.position = legend.pos) +
    geom_text_repel(data=label_up, aes(x=log2_FC, y=-log10(adj.p_values), label=Gene.Name), 
                    col=colpal[1], vjust=0.5, hjust=0.1, max.overlaps = overlaps[1]) + 
    geom_text_repel(data=label_down, aes(x=log2_FC, y=-log10(adj.p_values), label=Gene.Name), 
                    col=colpal[2], vjust=0.5, hjust=-0.1, max.overlaps = overlaps[2]) +
    scale_x_continuous(xnam, limits=c(xmin, xmax), 
                       breaks=seq(xmin, xmax, axis.breaks[1])) +
    scale_y_continuous(ynam, limits=c(0, ymax), breaks=seq(0, ymax, axis.breaks[2]))
  
  if(length(other_params) > 0) {
    for (param in other_params) {
      cov <- cov + do.call(param[[1]], param[-1])
    }
  }
  
  return(cov)
  
}

volcanoplot_func = function(ttres, categories = c("Group1", "Group2"), 
                            colpalette=c("#1B9E77", "#D95F02"), save=FALSE, filename="default", filedim = c(7, 5.6), 
                            fcut = c(0,0), shadcut=c(1,1), legend.pos=c(0.2, 0.9), ymax=NULL, xmax=NULL, xmin=NULL,
                            overlaps=c(20,20), axis.breaks=c(0.5,0.5), yaxis="fdr",
                            other_params = list(), ...) {
  ttres$alteration = paste0("Higher in ", categories[1])
  ttres$alteration[which(ttres$log2_FC<0)] = paste0("Higher in ", categories[2])
  
  if(ttres$method[1]=="Paired t-test") {
    xnam = "log2-fold change"
  } else if (ttres$method[1]=="Welch Two Sample t-test") {
    xnam = "log2-mean difference"
  } else {
    xnam = "log2-mean difference"
  }
  
  if(yaxis=="pval") {
    ttres$adj.p_values = ttres$p.value
  }
  
  if(yaxis=="fdr") {
    ynam = expression(paste("-log10 ", italic("q"), " values"))
  } else {
    ynam = expression(paste("-log10 ", italic("p"), " values")) 
    ttres$adj.p_values = ttres$p.value
  }
  
  ttres$shade = abs(ttres$log2_FC)*-log10(ttres$adj.p_values)
  seqv = c(0, unlist(lapply(list(seq(0.1,1,0.1)), 
                            function(x) quantile(ttres$shade, x))))
  ttres$shadcat = cut(ttres$shade, breaks=seqv)
  levels(ttres$shadcat) = c(1:(length(seqv)-1))
  
  label_up <- ttres |> filter(shade>shadcut[1]) |> filter(log2_FC>fcut[1])
  label_down <- ttres |> filter(shade>shadcut[2]) |> filter(log2_FC<fcut[2])
  
  if(is.null(xmin)) {xmin = floor(min(ttres$log2_FC))}
  if(is.null(xmax)) {xmax = ceiling(max(ttres$log2_FC))}
  if(is.null(ymax)) {ymax = ceiling(max(-log10(ttres$adj.p_values)))}
  
  if(all(sort(categories)==categories)) {
    colpal = colpalette
  } else {
    colpal = rev(colpalette)
  }
  
  cov = ggplot_func(ttres, label_up, label_down, colpalette, colpal, overlaps, xnam, ynam, 
                    xmin, xmax, ymax, legend.pos, axis.breaks, other_params, ...)
  
  if(save==TRUE) {
    if(filename=="default") {
      filename = paste0("volcano_plot_", categories[1], "_vs._", categories[2], ".pdf")
      ggsave(filename, cov, width=filedim[1], height=filedim[2], units="in", dpi=600)
    } else {
      ggsave(paste0("volcano_plot_", filename, ".pdf"), cov, width=filedim[1], height=filedim[2], units="in", dpi=600)
    }
  }
  
  return(cov)
  
}

##### Comparing different NMF groups #####
nmf_daps = function(nmf1, pt, cp2, adj.pval="bonferroni") {
  t.result = apply(cp2, 1, function(i) t.test(i[nmf1], i[pt]))
  results = lapply(t.result, unlist)
  results = as.data.frame(do.call(rbind, results))
  results[,c(1:8)] = apply(results[,c(1:8)], 2, as.numeric)
  results$adj.p_values <- p.adjust(results$p.value, method = adj.pval)
  results$gene_id = rownames(results)
  results$log2fc = results$`estimate.mean of x`-results$`estimate.mean of y`
  results <- results[order(results$log2fc, decreasing = T),]
  results = results[,c(14:15,1:13)]
  results
}