##### Change over time #####
##### Change over time in protein levels #####
##### Boxplots of proteins or pathways #####
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

##### Change pathways sum scores #####
calculate_path.score = function(xdf, db.sets, selected) {
  dpath = list()
  nprot = list()
  for(i in selected) {
    sel = db.sets$gene_symbol[which(db.sets$gs_name==i)]
    
    xdf.sel = as.data.frame(xdf[,which(colnames(xdf) %in% sel)])
    xdf.sel = na.omit(xdf.sel)
    if(ncol(xdf.sel)==1) {
      colnames(xdf.sel) = colnames(xdf)[which(colnames(xdf) %in% sel)]
    }
    
    if(is_empty(xdf.sel)) {
      dpath[[i]] = "No protein found"
    } else if(ncol(xdf.sel)==1) {
      dpath[[i]] = unlist(xdf.sel)
    } else {
      dpath[[i]] = apply(xdf.sel, 1, sum)
    }
    nprot[[i]] = ncol(xdf.sel)
    
    print(paste0("Run: ", match(i, selected), " done."))
  }
  
  xdf.path = as.data.frame(do.call(cbind, dpath[-grep("No protein", dpath)]))
  xdf.res = list(xdf.path, dpath[grep("No protein", dpath)], nprot)
  names(xdf.res) = c("Sum_scores", "Pathways_not_detected", "N_prot_per_pathway")
  return(xdf.res)
}

##### Change over time at all timepoints #####
trend_plasma = function(xdf, sdfx, selected, filename="default", filedim=c(6, 4.25), 
                        save=TRUE, plot=TRUE, stratify=c("cluster_nam", "macroradical_resection"), 
                        groups=list(c("PPG1", "PPG2"), c(0, 1)), 
                        error_bars = "std.error", conf.int = 0.95,
                        colpalette=c("#1B9E77", "#0072B2", "#D95F02", "#E69F00"), 
                        sdfid="Swim_ID", shapevals = c(24,21,24,21,25,22,25,22), xbreaks=seq(0, 460, 30),
                        tpdiffvars=c("diff_tp2_tp1", "diff_tp3_tp1", "diff_tp4_tp1"),
                        labids=c("Lab_ID_TP1", "Lab_ID_TP2", "Lab_ID_TP3", "Lab_ID_TP4"),
                        reference=FALSE, ref_ids=NULL, ref.chars=c(16, "#87CEEB")) {
  for(m in selected) {
    d = xdf[,m]
    names(d) = rownames(xdf)
    
    ntps = length(labids)
    tpnams = paste0("TP", 1:ntps)
    
    tplabels = list()
    for(j in tpnams) {
      tplabels[[j]] = rep(j, nrow(sdfx))
    }
    
    tptime = list()
    for(j in 1:length(tpdiffvars)) {
      tptime[[1]] = rep(0, nrow(sdfx))
      tptime[[1+j]] = sdfx[,tpdiffvars[[j]]]
    }
    
    tpvals = list()
    for(j in labids) {
      tpvals[[j]] = d[sdfx[,j]]
    }
    
    tpstrata = list()
    for(j in stratify) {
      tpstrata[[j]] = rep(sdfx[,j], ntps)
    }
    
    
    tpd = as.data.frame(cbind(rep(sdfx[,sdfid], ntps), unlist(tplabels),
                              unlist(tptime), unlist(tpvals),
                              do.call(cbind, tpstrata))
    )
    
    colnames(tpd) = c(sdfid, "TP", "Days_from_TP1", "Protein_value", stratify)
    tpd$Protein_value = as.numeric(tpd$Protein_value)
    tpd$Days_from_TP1 = as.numeric(tpd$Days_from_TP1)
    tpd[,sdfid] = factor(tpd[,sdfid], levels=unique(tpd[,sdfid]))
    
    times = tpnams
    vals = list()
    tpvals = list()
    med = list()
    st.dev = list()
    st.err = list()
    tpmed = list()
    tpst.dev = list()
    tpst.err = list()
    combos = list()
    cat = list()
    for(t in times) {
      if(length(stratify)>1) {
        combos[[t]] = apply(expand.grid(groups), 1, \(x) paste0(x, collapse="_&_"))
      } else {
        combos[[t]] = groups[[1]]
      }
      
      k = paste0(t, "_", 1:length(combos[[t]]))
      names(combos[[t]]) = k
      
      for(l in k) {
        if(length(stratify)>1) {
          # When stratified for two conditions
          cond_a = gsub("_&.+", "", combos[[t]][[l]])
          cond_b = gsub(".+_+", "", combos[[t]][[l]])
          
          vals[[l]] = tpd |> filter(TP==t) |> 
            filter(!!sym(stratify[1])==cond_a) |>
            filter(!!sym(stratify[2])==cond_b) |>
            select(Protein_value)
          
          tpvals[[l]] = tpd |> filter(TP==t) |> 
            filter(!!sym(stratify[1])==cond_a) |>
            filter(!!sym(stratify[2])==cond_b) |>
            select(Days_from_TP1)
          
          cat[[l]] = rep(combos[[t]][[l]], length(unlist(vals[[l]])))
          
        } else {
          # When stratified for one condition
          vals[[l]] = tpd |> filter(TP==t) |> 
            filter(!!sym(stratify[1])==combos[[t]][l]) |> 
            select(Protein_value)
          
          tpvals[[l]] = tpd |> filter(TP==t) |> 
            filter(!!sym(stratify[1])==combos[[t]][l]) |> 
            select(Days_from_TP1)
          
          cat[[l]] = rep(combos[[t]][[l]], length(unlist(vals[[l]])))
          
        }
      }
      
      for(j in k) {
        med[[j]] = mean(as.numeric(unlist(vals[[j]])), na.rm=T)
        st.dev[[j]] = sd(as.numeric(unlist(vals[[j]])), na.rm=T)
        st.err[[j]] = st.dev[[j]]/sqrt(length(na.omit(unlist(vals[[j]]))))
        tpmed[[j]] = mean(as.numeric(unlist(tpvals[[j]])), na.rm=T)
        tpst.dev[[j]] = sd(as.numeric(unlist(tpvals[[j]])), na.rm=T)
        tpst.err[[j]] = tpst.dev[[j]]/sqrt(length(na.omit(unlist(tpvals[[j]]))))
      }
    }
    
    tpd2 = as.data.frame(cbind(names(med), unlist(med), unlist(st.dev), unlist(st.err),
                               unlist(tpmed), unlist(tpst.dev), unlist(tpst.err)))
    tpd2[,-1] = apply(tpd2[,-1], 2, as.numeric)
    colnames(tpd2) = c("Category", "Prot_Mean", "Prot_SD", "Prot_SE",
                       "Time_Mean", "Time_SD", "Time_SE")
    tpd2$Group = unlist(combos)
    tpd$Group = unlist(cat)
    
    if(error_bars == "std.error") {
      ### Standard error of the mean
      tpd2$y_li = tpd2$Prot_Mean-tpd2$Prot_SE
      tpd2$y_ui = tpd2$Prot_Mean+tpd2$Prot_SE
      tpd2$x_li = tpd2$Time_Mean-tpd2$Time_SE
      tpd2$x_ui = tpd2$Time_Mean+tpd2$Time_SE
    } else if(error_bars == "conf.int") {
      zval = qnorm((1-conf.int)/2+conf.int)
      tpd2$y_li = tpd2$Prot_Mean-(zval*tpd2$Prot_SD)
      tpd2$y_ui = tpd2$Prot_Mean+(zval*tpd2$Prot_SD)
      tpd2$x_li = tpd2$Time_Mean-(zval*tpd2$Time_SD)
      tpd2$x_ui = tpd2$Time_Mean+(zval*tpd2$Time_SD)
    }
    
    s = summary(tpd$Protein_value)
    forder = as.numeric(as.factor(tpd$Group[1:nrow(sdfx)]))
    
    if(length(stratify)>1) {
      forder = c(forder, 1:4)
    } else {
      forder = c(forder, 1:2)
    }
    
    if(reference==TRUE) {
      nvals = xdf[ref_ids,m]
      ref.tpd = as.data.frame(cbind(c(0, max(xbreaks)), rep(mean(nvals, na.rm=T), 2), 
                          rep(sd(nvals, na.rm=T), 2), 
                      rep(sd(nvals, na.rm=T)/sqrt(length(na.omit(nvals))), 2)
                      ))
      colnames(ref.tpd) = c("Time", "Prot_Mean", "Prot_SD", "Prot_SE")
      
      if(error_bars == "std.error") {
        ### Standard error of the mean
        ref.tpd$y_li = ref.tpd$Prot_Mean-ref.tpd$Prot_SE
        ref.tpd$y_ui = ref.tpd$Prot_Mean+ref.tpd$Prot_SE
      } else if(error_bars == "conf.int") {
        zval = qnorm((1-conf.int)/2+conf.int)
        ref.tpd$y_li = ref.tpd$Prot_Mean-(zval*ref.tpd$Prot_SD)
        ref.tpd$y_ui = ref.tpd$Prot_Mean+(zval*ref.tpd$Prot_SD)
      }
      
    }
    
    
    # Fix the colours and the variable names
    trendplot = ggplot(tpd, aes(Days_from_TP1, Protein_value, 
                                shape=Group, fill=Swim_ID, col=Group)) +
      geom_point(size=1.5, alpha=0.2) +
      geom_line(linewidth=0.5, alpha=0.2) +
      theme_classic() +
      scale_color_manual(values=rep(colpalette, 2)) +
      scale_shape_manual(values=shapevals) +
      geom_ribbon(data=tpd2, aes(x=Time_Mean, ymin=y_li, ymax=y_ui, col=Group, fill=Group), 
                  inherit.aes = FALSE, alpha=0.4) +
      geom_point(data=tpd2, aes(x=Time_Mean, y=Prot_Mean, shape=Group, fill=Group, col=Group), 
                 size=3)  + 
      geom_line(data=tpd2, aes(x=Time_Mean, y=Prot_Mean, shape=Group, fill=Group, col=Group)) +
      scale_fill_manual(values=colpalette[forder]) +
      geom_linerange(data=tpd2, aes(x=Time_Mean, ymin=y_li, ymax=y_ui, col=Group), inherit.aes = FALSE) +
      geom_linerange(data=tpd2, aes(y=Prot_Mean, xmin=x_li, xmax=x_ui), inherit.aes = FALSE) +
      ggtitle(m) +
      scale_x_continuous("Days from TP1", breaks=xbreaks) +
      scale_y_continuous("Relative Protein levels (log2)", 
                         breaks=round(seq(s[1], s[6], length.out=10), 1)) +
      guides(fill = "none") +
      theme(legend.position = c(0.85, 0.85))
    
    if(reference==TRUE) {
      trendplot = trendplot +
        geom_point(data=ref.tpd, aes(x=Time, y=Prot_Mean),
                   shape=as.numeric(ref.chars[1]), fill=ref.chars[2], col=ref.chars[2],
                   size=3)  + 
        geom_line(data=ref.tpd, aes(x=Time, y=Prot_Mean),
                  shape=as.numeric(ref.chars[1]), fill=ref.chars[2], col=ref.chars[2],
                  linetype="dashed", inherit.aes = FALSE) +
        geom_ribbon(data=ref.tpd, aes(x=Time, ymin=y_li, ymax=y_ui), linetype="dashed", 
                    inherit.aes = FALSE, fill=ref.chars[2], col=ref.chars[2], alpha=0.4)
    }
    
    
    if(plot==TRUE) {
      print(trendplot)
    }
    
    if(save==TRUE) {
      stratnam = paste(stratify, collapse="_")
      if(filename=="default") {
        ggsave(paste0(m, "_protein_trend_per_", stratnam, ".pdf"), trendplot, 
               dpi=600, width=filedim[1], height=filedim[2])
      } else {
        ggsave(paste0(filename, m, "_protein_trend_per_", stratnam, ".pdf"), trendplot, 
               dpi=600, width=filedim[1], height=filedim[2])
      }
    }
    
    print(paste0("Run: ", match(m, selected), " is done. Protein: ", m))
  }
}
