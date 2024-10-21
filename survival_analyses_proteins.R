##### Survival analysis on proteins and pathways #####
##### Cox model for proteins, with and without adjustment #####
survival_proteins = function(osdf, sdfx, outcome_time="Survival_days", outcome_event="death", fdr.cut = 0.05,
                             categorical=FALSE, quantile.cut=0.5,
                             protein_list=NULL, adj.vars=NULL, save=FALSE, filename="default",
                             interact_protein=NULL) {
  if(!(all(rownames(osdf)==rownames(sdfx)))) {
    stop("The subject IDs in the protein matrix differ from the ones in the clinical dataset.\n  Recheck the datasets again.")
  }
  
  cox_m = list()
  if(is.null(protein_list)) {
    proteins = colnames(osdf)
  } else {
    proteins = protein_list[protein_list %in% colnames(osdf)]
  }
  
  protnams = proteins
  proteins = gsub("[-]", ".", proteins)
  colnames(osdf) = gsub("[-]", ".", colnames(osdf))
  
  combdf = cbind(osdf, sdfx)
  
  for(i in proteins) {
    if(categorical==TRUE) {
      combdf[,i] <- ifelse((combdf[,i]>quantile(combdf[,i], quantile.cut, na.rm=T))==T, 1, 0)
    }
    
    if(is.null(adj.vars)) {
      my.formula = as.formula(paste("Surv(", outcome_time, ", ", outcome_event, ") ~ ", i))
    } else {
      if(is.null(interact_protein)) {
        my.formula = as.formula(paste("Surv(Survival_days, death) ~ ", i, "+", 
                                      paste(adj.vars, collapse = "+")))
      } else {
        my.formula = as.formula(paste("Surv(Survival_days, death) ~ ", i, "+", 
                                      paste(adj.vars, collapse = "+"), "+",
                                      interact_protein, ":", i))
      }
    }
    
    cox_m[[i]] = coxph(my.formula,
                       data=combdf)
  }
  sumcox = lapply(cox_m, function(x) summary(x))
  coefs = lapply(sumcox, function(x) c(as.numeric(x$coefficients), as.numeric(x$conf.int), as.numeric(x$logtest)))
  
  coxres = as.data.frame(cbind(protnams, do.call(rbind, coefs)))
  suffix = c(colnames(sumcox[[1]]$coefficients), colnames(sumcox[[1]]$conf.int))
  if(!is.null(interact_protein)) {
    pref = c("Protein", adj.vars, paste0(interact_protein, ":", "Protein"))
  } else {
    pref = c("Protein", adj.vars)
  }
  suffix[6] = "HR"
  nams = list()
  for(i in suffix) {nams[[i]] = paste0(pref, "_", i)}
  colnames(coxres)[-1] = c(unlist(nams), paste0("Loglik_", names(sumcox[[1]]$logtest)))
  colnames(coxres)[1] = "Protein"
  
  coxres[,-1] = apply(coxres[,-1], 2, as.numeric)
  coxres$adj.p_values = p.adjust(coxres$Loglik_pvalue, method="fdr")
  coxres = coxres[order(coxres$`Protein_exp(coef)`, decreasing = T),]
  
  if(!is.null(adj.vars)) {
    pvals = t(apply(coxres[,grep("Pr[(]", colnames(coxres))], 1, function(x) x<0.05))
    coxres$Pval_AllSign = apply(pvals, 1, function(x) all(x==T))
    coxres_sign = coxres |> filter(adj.p_values<fdr.cut) |> filter(Pval_AllSign==T)
  } else {
    coxres_sign = coxres |> filter(adj.p_values<fdr.cut) |> filter(`Protein_Pr(>|z|)`<0.05)
  }
  
  if(save==TRUE) {
    if(filename=="default") {
      write.csv(coxres, paste0("cox_models_", outcome_time, "_", outcome_event, "_", paste(pref, collapse="_"), "_all.csv"))
      write.csv(coxres_sign, paste0("cox_models_", outcome_time, "_", outcome_event, "_", paste(pref, collapse="_"), 
                                    "_sign_at_", fdr.cut*100,"_%FDR.csv"))
    } else {
      write.csv(coxres, paste0("cox_models_", filename, "_all.csv"))
      write.csv(coxres_sign, paste0("cox_models_", filename, "_sign_at_", fdr.cut*100,"_%FDR.csv"))
    }
  }
  
  return(coxres)
}

##### Cox model for pathways, with and without adjustment #####
survival_pathways = function(osdf, sdfx, db.sets, sets=NULL, min.prot=3, score.type="Set_Higher", quant.cut=0.5, 
                             outcome_time="Survival_days", outcome_event="death", adj.vars=NULL, interact_protein=NULL,
                             fdr.cut = 0.05, save=FALSE, filename="default") {
  if(!(all(rownames(osdf)==rownames(sdfx)))) {
    stop("The subject IDs in the protein matrix differ from the ones in the clinical dataset.\n  Recheck the datasets again.")
  }
  
  if(is.null(sets)) {
    sets = unique(db.sets$gs_name)
  }
  
  # Define the lists for saving results
  cox_m = list()
  protset_in = list()
  patscore = list()
  patscore_cat = list()
  
  for(i in sets) {
    sel = db.sets$gene_symbol[which(db.sets$gs_name==i)]
    xdf = as.data.frame(osdf[,which(colnames(osdf) %in% sel)])
    if(ncol(xdf)==1) {
      colnames(xdf) = colnames(osdf)[which(colnames(osdf) %in% sel)]
    }
    
    protset_in[[i]] = list(sel, colnames(xdf))
    
    if(is_empty(xdf)) {
      cox_m[[i]] = "No protein found"
    } else if(ncol(xdf)>=min.prot) {
      patscore[[i]] <- apply(xdf, 1, sum)
      patscore_cat[[i]] <- ifelse((patscore[[i]]>quantile(patscore[[i]], quant.cut))==TRUE, 1, 0)
      
      kmdata = as.data.frame(cbind(sdfx, as.numeric(unlist(patscore[[i]])), as.numeric(unlist(patscore_cat[[i]]))))
      ncol(kmdata)
      colnames(kmdata)[(ncol(kmdata)-1):ncol(kmdata)] = c("Score", "Set_Higher")
      
      if(is.null(adj.vars)) {
        my.formula = as.formula(paste("Surv(", outcome_time, ", ", outcome_event, ") ~ ", score.type))
        cox_m[[i]] = coxph(my.formula, data=kmdata)
      } else {
        if(is.null(interact_protein)) {
          my.formula2 = as.formula(paste("Surv(Survival_days, death) ~ ", score.type, "+", 
                                         paste(adj.vars, collapse = "+")))
        } else {
          my.formula2 = as.formula(paste("Surv(Survival_days, death) ~ ", score.type, "+", 
                                         paste(adj.vars, collapse = "+"), "+",
                                         interact_protein, ":", score.type))
        }
        
        cox_m[[i]] = coxph(my.formula2, data=kmdata)
      }
    }
  }
  cox_m = cox_m[-grep("No protein", cox_m)]
  
  sumcox = lapply(cox_m, function(x) summary(x))
  coefs = lapply(sumcox, function(x) c(as.numeric(x$coefficients), as.numeric(x$conf.int), as.numeric(x$logtest)))
  
  coxres = as.data.frame(cbind(names(cox_m), do.call(rbind, coefs)))
  suffix = c(colnames(sumcox[[1]]$coefficients), colnames(sumcox[[1]]$conf.int))
  if(is.null(interact_protein)) {
    pref = c("Pathway", adj.vars)
  } else {
    pref = c("Pathway", adj.vars, paste0(interact_protein, ":", "Pathway"))
  }
  suffix[6] = "HR"
  nams = list()
  for(i in suffix) {nams[[i]] = paste0(pref, "_", i)}
  colnames(coxres)[-1] = c(unlist(nams), paste0("Loglik_", names(sumcox[[1]]$logtest)))
  colnames(coxres)[1] = "Pathway"
  
  coxres[,-1] = apply(coxres[,-1], 2, as.numeric)
  coxres$adj.p_values = p.adjust(coxres$Loglik_pvalue, method="fdr")
  coxres = coxres[order(coxres$`Pathway_exp(coef)`, decreasing = T),]
  
  if(!is.null(adj.vars)) {
    pvals = t(apply(coxres[,grep("Pr[(]", colnames(coxres))], 1, function(x) x<0.05))
    coxres$Pval_AllSign = apply(pvals, 1, function(x) all(x==T))
    coxres_sign = coxres |> filter(adj.p_values<fdr.cut) |> filter(Pval_AllSign==T)
  } else {
    coxres_sign = coxres |> filter(adj.p_values<fdr.cut) |> filter(`Pathway_Pr(>|z|)`<0.05)
  }
  
  if(save==TRUE) {
    if(filename=="default") {
      write.csv(coxres, paste0("cox_models_", outcome_time, "_", outcome_event, "_", paste(pref, collapse="_"), "_all.csv"))
      write.csv(coxres_sign, paste0("cox_models_", outcome_time, "_", outcome_event, "_", paste(pref, collapse="_"), 
                                    "_sign_at_", fdr.cut*100,"_%FDR.csv"))
    } else {
      write.csv(coxres, paste0("cox_models_", filename, "_all.csv"))
      write.csv(coxres_sign, paste0("cox_models_", filename, "_sign_at_", fdr.cut*100,"_%FDR.csv"))
    }
  }
  
  patscore.df = as.data.frame(do.call(cbind, patscore))
  patscore_cat.df = as.data.frame(do.call(cbind, patscore_cat))
  
  all_results = list(protset_in, patscore.df, patscore_cat.df, coxres)
  names(all_results) = c("Proteins_per_Set", "Pathways_Score", "Pathways_Category", "Cox_models_result")
  
  return(all_results)
}

##### Kaplan-Meier plot for pathways, with and without adjustment #####
km_pathways = function(patscore_cat.df, sdfx, coxlabels, pathways=NULL, stratify=NULL, adj.vars=NULL, interact_protein=NULL, 
                       timeunit="days", timeconvert="none", outcome_time="Survival_days", outcome_event="death",
                       xnam="Time to event", xcut=NULL, 
                       filedim=c(7,6), filedim.strat=c(7,6), filedim.multi=c(7,6), ...) {
  if(!(all(rownames(osdf)==rownames(sdfx)))) {
    stop("The subject IDs in the protein matrix differ from the ones in the clinical dataset.\n  Recheck the datasets again.")
  }
  
  if(is.null(pathways)) {
    pathways = colnames(patscore_cat.df)
  }
  
  if(timeunit=="days" & timeconvert=="months") {
    sdfx[,outcome_time] = sdfx[,outcome_time]/(((365/12*3)+(366/12))/4) # Adjusting for the leap year
  } else if (timeunit=="days" & timeconvert=="weeks") {
    sdfx[,outcome_time] = sdfx[,outcome_time]/7
  } else if(timeunit=="days" & timeconvert=="years") {
    sdfx[,outcome_time] = sdfx[,outcome_time]/((365*3+366)/4) # Adjusting for the leap year
  } else if(timeunit=="months" & timeconvert=="years") {
    sdfx[,outcome_time] = sdfx[,outcome_time]/12
  } else if(timeunit=="months" & timeconvert=="weeks") {
    sdfx[,outcome_time] = sdfx[,outcome_time]*(((365/12*3)+(366/12))/4)/7 # Adjusting for the leap year
  } else if(timeunit=="months" & timeconvert=="days") {
    sdfx[,outcome_time] = sdfx[,outcome_time]*(((365/12*3)+(366/12))/4) # Adjusting for the leap year
  } else if(timeunit=="years" & timeconvert=="months") {
    sdfx[,outcome_time] = sdfx[,outcome_time]*12
  } else if(timeunit=="years" & timeconvert=="weeks") {
    sdfx[,outcome_time] = sdfx[,outcome_time]*365.25/7
  } else if(timeunit=="years" & timeconvert=="days") {
    sdfx[,outcome_time] = sdfx[,outcome_time]*365.25
  }
  
  for(i in pathways) {
    kmdata = as.data.frame(cbind(sdfx, as.numeric(patscore_cat.df[,i])))
    colnames(kmdata)[ncol(kmdata)] = "Set_Higher"
    
    my.formula = as.formula(paste("Surv(", outcome_time, ", ", outcome_event, ") ~ Set_Higher"))
    fit <- survfit(my.formula, data = kmdata)
    
    call2 = paste0("survfit(formula = Surv(", outcome_time, ", ", outcome_event, ") ~  Set_Higher", ", data = kmdata)")
    fit$call = parse(text = call2)[[1]]
    
    xmax = ceiling(max(as.numeric(kmdata[,outcome_time])))
    pval.coord = c(xmax/2, 0.25)
    
    if(is.null(xcut)) {
      xcut = round(seq(0, xmax, length.out=10)[2])
    }
    
    km_os <- ggsurvplot(
      fit,                     # survfit object with calculated statistics.
      data = kmdata,             # data used to fit survival curves.
      risk.table = TRUE,       # show risk table.
      pval = TRUE,             # show p-value of log-rank test.
      xlab = xnam,   # customize X axis label.
      xlim = c(0, xmax),         # present narrower X axis, but not affect
      # survival estimates.
      break.time.by = xcut,     # break X axis in time intervals by 500.
      ggtheme = theme_classic(), # customize plot and risk table with a theme.
      risk.table.y.text.col = T,# colour risk table text annotations.
      risk.table.height = 0.25, # the height of the risk table
      risk.table.y.text = FALSE,# show bars instead of names in text annotations
      surv.median.line = "hv",  # add the median survival pointer.
      pval.coord = pval.coord,
      ...
    )
    
    pdf(paste0(i, "_KM_univariable.pdf"), width=filedim[1], height=filedim[2])
    print(km_os, newpage=FALSE)
    dev.off()
    
    print(paste0("Univariate survival plot for set: --- ", i, " --- done!"))
    
    if(!is.null(stratify)) {
      my.formula3 = as.formula(paste("Surv(", outcome_time, ", ", outcome_event, ") ~ Set_Higher +", stratify))
      fit <- survfit(my.formula3, data = kmdata)
      
      call2 = paste0("survfit(formula = Surv(", outcome_time, ", ", outcome_event, ") ~  Set_Higher + ", stratify, ", data = kmdata)")
      fit$call = parse(text = call2)[[1]]
      
      xmax = ceiling(max(kmdata[,outcome_time]))
      pval.coord = c(xmax/2, 0.5)
      if(is.null(xcut)) {
        xcut = round(seq(0, xmax, length.out=10)[2])
      }
      
      km_os.strat <- ggsurvplot(
        fit,                     # survfit object with calculated statistics.
        data = kmdata,             # data used to fit survival curves.
        risk.table = TRUE,       # show risk table.
        pval = TRUE,             # show p-value of log-rank test.
        xlab = xnam,   # customize X axis label.
        xlim = c(0, xmax),         # present narrower X axis, but not affect
        # survival estimates.
        break.time.by = xcut,     # break X axis in time intervals by 500.
        ggtheme = theme_classic(), # customize plot and risk table with a theme.
        risk.table.y.text.col = T,# colour risk table text annotations.
        risk.table.height = 0.25, # the height of the risk table
        risk.table.y.text = FALSE,# show bars instead of names in text annotations
        surv.median.line = "hv",  # add the median survival pointer.
        pval.coord = pval.coord
      )
      
      pdf(paste0(i, "_KM_stratified_for_", stratify,"_.pdf"), 
          width=filedim.strat[1], height=filedim.strat[2])
      print(km_os.strat, newpage=FALSE)
      dev.off()
      
      print(paste0("Stratified survival plot for set: --- ", i, " --- done!"))
    }
    
    if(!is.null(adj.vars)) {
      if(is.null(interact_protein)) {
        my.formula2 = as.formula(paste("Surv(Survival_days, death) ~ Set_Higher +", 
                                       paste(adj.vars, collapse = "+")))
      } else {
        my.formula2 = as.formula(paste("Surv(Survival_days, death) ~ Set_Higher +", 
                                       paste(adj.vars, collapse = "+"), "+",
                                       interact_protein, ":Set_Higher"))
      }
      
      cox_multi = summary(coxph(my.formula2, data=kmdata))
      
      cdf = as.data.frame(cox_multi$conf.int)
      var.nams = c("Set_Higher", unlist(coxlabels[match(rownames(cdf), names(coxlabels))]))
      cdf$Variable = var.nams
      cdf$Variable = as.factor(cdf$Variable)
      cdf$Variable = ordered(cdf$Variable, levels=rev(var.nams))
      cdf$Label = paste0("HR = ", round(cdf$`exp(coef)`, 3), " [", round(cdf$`lower .95`, 3), 
                         "-", round(cdf$`upper .95`, 3), "]")
      vmax = ceiling(max(cdf$`upper .95`))
      slen = round(seq(0, vmax, length.out=10), 1)
      
      cm_p = ggplot(cdf, aes(`exp(coef)`, Variable)) +
        geom_vline(xintercept=1, col="black", linetype="dashed") +
        geom_point(size=1.5) +
        geom_linerange(xmin=cdf$`lower .95`, xmax=cdf$`upper .95`) +
        geom_text(aes(x=cdf$`exp(coef)`, y=cdf$Variable),
                  label=cdf$Label, hjust=0, vjust=-0.5, size=2.5) +
        theme_classic() +
        ylab("") +
        scale_x_continuous("HR with 95% CI", 
                           limits=c(0,vmax), breaks=slen)  +
        annotate("text", x=(vmax/2), y=0.75, size=2.5, label=paste("LRT, p =",
                                                                   formatC(cox_multi$logtest[3], 
                                                                           format = "e", digits = 3))) +
        theme(axis.text.y=element_text(hjust=1), 
              text = element_text(size=10)) +
        theme(plot.background = element_rect(color = "black"))
      
      tosav = km_os$plot + annotation_custom(ggplotGrob(cm_p), xmin=30, xmax=130, ymin=0.5)
      tosav2 = ggarrange(tosav, km_os$table, nrow=2, ncol=1, heights = c(2, 0.7), align = "v")
      pdf(paste0(i, "_KM_univariable_and_CoxPH.pdf"), width=filedim.multi[1], height=filedim.multi[2])
      print(tosav2, newpage=FALSE)
      dev.off()
      
      print(paste0("Multivariable survival plot for set: --- ", i, " --- done!"))
    }
    print(paste0("Run: ", match(i, pathways), " out of ", length(pathways), " done!"))
  }
}
