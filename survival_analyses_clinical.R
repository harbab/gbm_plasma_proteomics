##### Survival analysis on clinical data #####
##### Kaplan-Meier curves plot #####
km_plot = function(kmdata, var_test, timeunit="days", timeconvert="none", 
                   outcome_time="Survival_days", outcome_event="death",
                   xnam="Time to event", xcut=NULL, 
                   plot=FALSE, filename="default", filedim=c(7,6), ...) {
  kmdata[,outcome_time] = as.numeric(kmdata[,outcome_time])
  
  if(timeunit=="days" & timeconvert=="months") {
    kmdata[,outcome_time] = kmdata[,outcome_time]/(((365/12*3)+(366/12))/4) # Adjusting for the leap year
  } else if (timeunit=="days" & timeconvert=="weeks") {
    kmdata[,outcome_time] = kmdata[,outcome_time]/7
  } else if(timeunit=="days" & timeconvert=="years") {
    kmdata[,outcome_time] = kmdata[,outcome_time]/((365*3+366)/4) # Adjusting for the leap year
  } else if(timeunit=="months" & timeconvert=="years") {
    kmdata[,outcome_time] = kmdata[,outcome_time]/12
  } else if(timeunit=="months" & timeconvert=="weeks") {
    kmdata[,outcome_time] = kmdata[,outcome_time]*(((365/12*3)+(366/12))/4)/7 # Adjusting for the leap year
  } else if(timeunit=="months" & timeconvert=="days") {
    kmdata[,outcome_time] = kmdata[,outcome_time]*(((365/12*3)+(366/12))/4) # Adjusting for the leap year
  } else if(timeunit=="years" & timeconvert=="months") {
    kmdata[,outcome_time] = kmdata[,outcome_time]*12
  } else if(timeunit=="years" & timeconvert=="weeks") {
    kmdata[,outcome_time] = kmdata[,outcome_time]*365.25/7
  } else if(timeunit=="years" & timeconvert=="days") {
    kmdata[,outcome_time] = kmdata[,outcome_time]*365.25
  } 
  
  if(length(var_test)==1) {
    my.formula = as.formula(paste("Surv(", outcome_time, ",", outcome_event, ") ~", var_test))
  } else {
    my.formula = as.formula(paste("Surv(", outcome_time, ", ", outcome_event, ") ~ ", 
                                  paste(var_test, collapse=" + ")))
  }
  
  fit <- survfit(my.formula, data = kmdata)
  if(length(var_test)==1) {
    call2 = paste0("survfit(formula = Surv(", outcome_time, ", ", outcome_event, ") ~ ", var_test, ", data = kmdata)")
  } else {
    call2 = paste0("survfit(formula = Surv(", outcome_time, ", ", outcome_event, ") ~ ", 
                   paste0(var_test, collapse = " + "), ", data = kmdata)")
  }
  
  fit$call = parse(text = call2)[[1]]
  
  xmax = ceiling(max(kmdata[,outcome_time]))
  pval.coord = c(xmax/2, 0.5)
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
  if(plot==T) {
    if(filename=="default") {
      filename = paste("KM", outcome_time, outcome_event, paste0(var_test, collapse="_"), sep="_")
    }
    
    pdf(paste0(filename, ".pdf"), width=filedim[1], height=filedim[2])
    print(km_os, newpage=FALSE)
    dev.off()
  }
  return(km_os)
}

##### Univariate Cox models for clinical variables #####
univariate_cox = function(sdfx, vars, outcome_time, outcome_event, save=FALSE, filename="default") {
  cox_r = list()
  cox_sum = list()
  for(i in vars) {
    my.formula = as.formula(paste("Surv(", outcome_time, ",", outcome_event, ") ~ ", i))
    cox_r[[i]] = coxph(my.formula, sdfx)
    c2 = summary(cox_r[[i]])
    cox_sum[[i]] = c(as.numeric(c2$coefficients), as.numeric(c2$logtest), as.numeric(c2$conf.int))
  }
  
  cox_res = as.data.frame(do.call(rbind, cox_sum))
  colnames(cox_res) = c(colnames(c2$coefficients), names(c2$logtest), colnames(c2$conf.int))
  cox_res$adj.pval = p.adjust(cox_res$`Pr(>|z|)`, method="fdr")
  cox_res$pval.sign = cox_res$`Pr(>|z|)`<0.05
  cox_res$adj.pval.sign = cox_res$adj.pval<0.05
  
  if(filename=="default" & save==TRUE) {
    filename = paste("univariate_cox_models", outcome_time, outcome_event ,"clinical_variables.csv", sep="_")
    write.csv(cox_res, filename)
  }
  return(cox_res)
}

univariate_cox_inter = function(sdfx, vars, outcome_time, outcome_event, save=FALSE, filename="default") {
  terms = list()
  cox_sumi = list()
  for(i in 1:length(vars)) {
    # Change * to : if you want only the interaction terms
    terms[[i]] = paste0(vars[i], ":", vars[-i])
    
    resi = list()
    for(z in 1:length(terms[[i]])) {
      my.formula = as.formula(paste("Surv(", outcome_time, ",", outcome_event, ") ~ ", terms[[i]][[z]]))
      c1 =  coxph(my.formula, data=sdfx)
      c2 = summary(c1)
      resi[[z]] = c(as.numeric(c2$coefficients), as.numeric(c2$logtest), as.numeric(c2$conf.int))
    }
    
    cox_sumi[[i]] = as.data.frame(do.call(rbind, resi))
    
  }
  
  # Extract results for the interaction terms without the variables
  for(i in 1:length(vars)) {
    rownames(cox_sumi[[i]]) = terms[[i]]
    pref = c("Var1:Var2")
    suf = colnames(c2$coefficients)
    suf2 = colnames(c2$conf.int)
    colnames(cox_sumi[[i]]) = c(paste0(pref, "_", suf[1]), paste0(pref, "_", suf[2]), 
                                paste0(pref, "_", suf[3]), paste0(pref, "_", suf[4]),
                                paste0(pref, "_", suf[5]), names(c2$logtest),
                                paste0(pref, "_", "HR"), paste0(pref, "_", suf2[2]), 
                                paste0(pref, "_", suf2[3]), paste0(pref, "_", suf2[4]))
  }
  
  cox_intres = as.data.frame(do.call(rbind, cox_sumi))
  cox_intres$adj.pval = p.adjust(cox_intres$pvalue, method="fdr")
  cox_intres$pval.sign = cox_intres$`Var1:Var2_Pr(>|z|)`<0.05 & cox_intres$pvalue<0.05
  cox_intres$adj.pval.sign = cox_intres$adj.pval<0.05
  
  if(filename=="default" & save==TRUE) {
    filename = paste("interactions_cox_models", outcome_time, outcome_event ,"clinical_variables.csv", sep="_")
    write.csv(cox_intres, filename)
  }
  return(cox_intres)
}


##### Multivariate Cox model function #####
multivariate_cox = function(sdfx, var_test, pcutoff=0.05, outcome_time="Survival_days", outcome_event="death") {
  
  my.formula = as.formula(paste("Surv(", outcome_time, ", ", outcome_event, ") ~ ", 
                                paste(var_test, collapse=" + ")))
  
  cox_m = coxph(my.formula, data=sdfx)
  s = summary(cox_m)
  
  repeat{
    s2 = s$coefficients[,5];
    if(any(s2>pcutoff)) {
      torem = names(which(s2==max(s2)))
      var_test = var_test[-which(var_test==torem)]
      
      my.formula = as.formula(paste0("Surv(", outcome_time, ", ", outcome_event, ") ~ ", 
                                     paste(var_test, collapse=" + ")))
      cox_m = coxph(my.formula, data=sdfx)
      s = summary(cox_m)
    } else {
      cat("The significant variables are: \n", paste(
      rownames(s$coefficients[which(s$coefficients[,5]<0.05),]), collapse=", "), "\n    
          ", "\n The Cox regression results are as follows: \n")
      print(s$coefficients)
      return(s)
      break;}
  }
}

##### Cox model forest plot #####
forestplot_cox = function(coxres, filename, coxlabels,
                          filedim=c(6,2.5), xbreak=1, pval.pos=c(3,1.5), xmax=NULL) {
  coxres.sel = coxres$conf.int[,c(1, 3:4)]
  
  if(is.null(nrow(coxres.sel))) {
    cmres = as.data.frame(t(coxres.sel))
  } else {
    cmres = coxres.sel
    cmres = as.data.frame(apply(cmres, 2, as.numeric))
  }

  cmres$variable = rownames(coxres$coefficients)
  cmres  = cmres[order(cmres$`exp(coef)`, decreasing = T),]
  cmres$order = c(nrow(cmres):1)
  cmres$labs = paste(round(cmres$`exp(coef)`, 3), " [", round(cmres$`lower .95`, 3), 
                     "-", round(cmres$`upper .95`, 3), "]", sep="")
  cmres$label = unlist(coxlabels[match(cmres$variable, names(coxlabels))])
  lrt = paste("LRT, p =", formatC(coxres$logtest[3], format = "e", digits = 3))
  if(is.null(xmax)) {xmax = ceiling(max(cmres$`upper .95`))}
  
  cm_plot = ggplot(cmres, aes(y=as.factor(order), x=`exp(coef)`)) +
    geom_vline(xintercept=1, col="black", linetype="dashed") +
    geom_point(size=2.5) +
    geom_text(data=cmres, aes(y=as.factor(order), x=`exp(coef)`),
              label=cmres$labs, hjust=0, vjust=-0.5) +
    scale_y_discrete("variable", labels=rev(cmres$label)) +
    geom_linerange(xmin=cmres$`lower .95`, xmax=cmres$`upper .95`) +
    theme_classic() +
    scale_x_continuous("HR with 95% CI", limits=c(0,xmax), breaks=seq(0, xmax, xbreak))  +
    annotate("text", x=pval.pos[1], y=pval.pos[2], label=lrt) +
    theme(axis.text.y=element_text(hjust=0))
  
  ggsave(paste0(filename, ".pdf"), cm_plot, dpi=600, width=filedim[1], height=filedim[2])
  cm_plot
}
