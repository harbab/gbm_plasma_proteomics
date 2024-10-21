# Make a plot for protein enrichment
identify_plasma = function(proteins_in_sets, sel.df, x.var, y.var, 
                           point_params = list(), text_params = list(), other_params = list(), ...) {
  if(y.var=="adj.p_values") {
    proteins_in_sets$neglog10.adj.pval = -log10(proteins_in_sets$adj.p_values)
    sel.df$neglog10.adj.pval = -log10(sel.df$adj.p_values)
    y.var = "neglog10.adj.pval"
  }
  xmax = ceiling(max(proteins_in_sets[,x.var]))
  xcut = round(seq(0, xmax, length.out=10))
  ymax = ceiling(max(proteins_in_sets[,y.var]))
  ycut = round(seq(0, ymax, length.out=10))
  
  identplot = ggplot(proteins_in_sets, aes_string(sym(x.var), sym(y.var), 
                                                  alpha=sym("shadcat"), col=sym("database"))) +
    do.call(geom_point, c(list(...), point_params)) +
    theme_classic()  +
    (if(length(other_params) > 0) do.call(other_params[[1]], other_params[-1]) else NULL) +
    do.call(geom_text_repel, c(list(data = sel.df, aes_string(x = x.var, y = y.var, 
                                                              label = "gene_set"), size = 3), text_params)) +
    scale_x_continuous(x.var, 
                       breaks=xcut) +
    scale_y_continuous(y.var, breaks=ycut) +
    scale_alpha_continuous(breaks=seq(0,10,0.5), guide="none")
  return(identplot)
}

plasma_enriched = function(proteins_in_sets, prot.pool, x.var="Coefficient", y.var="adj.p_values",
                           pval.method="max", estimate.method="mean",
                           save=FALSE, plot=TRUE, filename="default", filedim=c(8, 6),
                           prot.cut=3, niter=30, limit=c(10, 0.01), 
                           sel.database=c("GBM", "Glioma", "KEGG", "Brain", "Hallmark", "Reactome"), 
                           point_params = list(), text_params = list(), other_params = list(), ...) {
  proteins_in_sets = as.data.frame(proteins_in_sets)
  proteins_in_sets$N_setprot = as.numeric(lapply(proteins_in_sets$proteins_in_set, length))
  proteins_in_sets$N_plasmaprot = as.numeric(lapply(proteins_in_sets$proteins_in_plasma, length))
  proteins_in_sets = proteins_in_sets |> filter(N_plasmaprot>=prot.cut)
  
  uniqpr = list()
  pvals = list()
  estimate = list()
  for(i in 1:nrow(proteins_in_sets)) {
    uniqpr[[i]] = setdiff(unlist(proteins_in_sets$proteins_in_plasma[i]),
                          unique(unlist(proteins_in_sets$proteins_in_plasma[-i])))
    
    n1 = proteins_in_sets$N_setprot[i]
    n2 = proteins_in_sets$N_plasmaprot[i]
    
    iterations = list()
    fisher.perm = list()
    for(j in 1:niter) {
      backlist = sample(prot.pool, n1)
      simlist = sample(prot.pool, n2)
      
      n3 = length(which(simlist %in% backlist==T))
      t = rbind(c(n2, n1-n2), c(n3, n1-n3))
      if(any(t==0)) {
        t = t+1
      }
      fisher.perm[[j]] = fisher.test(t)
    }
    
    if(pval.method=="max") {
      pvals[[i]] = max(unlist(lapply(fisher.perm, \(x) x$p.value)))
    } else if(pval.method=="min") {
      pvals[[i]] = min(unlist(lapply(fisher.perm, \(x) x$p.value)))
    }
    
    if(estimate.method=="mean") {
      estimate[[i]] = mean(unlist(lapply(fisher.perm, \(x) x$estimate)))
    } else if(estimate.method=="median") {
      estimate[[i]] = median(unlist(lapply(fisher.perm, \(x) x$estimate)))
    } else if(estimate.method=="max") {
      estimate[[i]] = max(unlist(lapply(fisher.perm, \(x) x$estimate)))
    } else if(estimate.method=="min") {
      estimate[[i]] = min(unlist(lapply(fisher.perm, \(x) x$estimate)))
    }
    
    print(paste("Permutation test", i, "finished."))
    
  }
  print("All permutation tests finished.")
  
  proteins_in_sets$unique_per_set = uniqpr
  proteins_in_sets$N_unique = as.numeric(lapply(proteins_in_sets$unique_per_set, length))
  proteins_in_sets$OR = as.numeric(unlist(estimate))
  proteins_in_sets$p.value = as.numeric(unlist(pvals))
  proteins_in_sets$adj.p_values = p.adjust(proteins_in_sets$p.value, method="fdr")
  
  proteins_in_sets$Perc_Plasma_prots = proteins_in_sets$N_plasmaprot/proteins_in_sets$N_setprot*100
  proteins_in_sets$Coefficient = (proteins_in_sets$Perc_Plasma_prots/100)*proteins_in_sets$N_plasmaprot
  proteins_in_sets$ORpval_Coefficient = proteins_in_sets$OR*-log10(proteins_in_sets$adj.p_values)
  
  if(plot==TRUE) {
    seqv = c(0, unlist(lapply(list(seq(0.1,1,0.1)), 
                              function(x) quantile(proteins_in_sets$ORpval_Coefficient, x))))
    proteins_in_sets$shadcat = cut(proteins_in_sets$ORpval_Coefficient, breaks=seqv)
    levels(proteins_in_sets$shadcat) = c(1:(length(seqv)-1))
    proteins_in_sets$shadcat = as.numeric(proteins_in_sets$shadcat)
    
    proteins_in_sets = proteins_in_sets[order(proteins_in_sets$ORpval_Coefficient, decreasing = T),]
    
    sel.df = proteins_in_sets |> filter(!!sym(x.var)>limit[1]) |> filter(!!sym(y.var)<limit[2]) |> 
      filter(database %in% c(sel.database))
    sel.df$gene_set = gsub("REACTOME_|KEGG_|HALLMARK_", "", sel.df$gene_set)
    sel.df$gene_set = gsub("_", " ", sel.df$gene_set)
    sel.df$gene_set = str_to_title(sel.df$gene_set)
    
    identplot = identify_plasma(proteins_in_sets, sel.df, x.var, y.var,
                    point_params, text_params, other_params)
    print(identplot)
  }
  
  proteins_in_sets = as.data.frame(apply(proteins_in_sets, 2, \(x) sapply(x, paste, collapse=", ")))
  proteins_in_sets[,c(5, 6, 8:15)] = apply(proteins_in_sets[,c(5, 6, 8:15)], 2, as.numeric)
  
  if(save==TRUE) {
    
    if(filename=="default") {
      write.csv(proteins_in_sets, paste("pathways_identification", paste(sel.database, collapse="_"), 
                                        estimate.method, x.var, pval.method, y.var, niter, "iterations.csv", sep="_"))
    } else {
      write.csv(proteins_in_sets, file=paste(filename, paste(sel.database, collapse="_"), 
                                             estimate.method, x.var, pval.method, y.var, niter, "iterations.csv", sep="_"))
    }
    
    if(filename=="default") {
      ggsave(paste("pathways_identification", paste(sel.database, collapse="_"), 
                   estimate.method, x.var, pval.method, y.var, niter, "iterations.pdf", sep="_"), 
             identplot, width=filedim[1], height=filedim[2])
    } else {
      ggsave(paste(filename, paste(sel.database, collapse="_"), estimate.method, x.var, 
                   pval.method, y.var, niter, "iterations.pdf", sep="_"), 
             identplot, width=filedim[1], height=filedim[2])
    }
  }
  
  return(proteins_in_sets)
}
