
##### Radiology analyses functions #####
##### Explore proteins' relations with radiological parameters #####
radio_proteome = function(sel.df, sel.var, adj.var) {
  corl = list()
  glm.res = list()
  glm.res_sign = list()
  
  for(j in sel.var) {
    
    corl = list()
    for(i in colnames(xdfsel)) {
      sum = summary(glm(as.formula(paste0(j, " ~ ", paste(adj.var, collapse="+"), " + ", i)), data=sel.df))
      corl[[i]] = as.numeric(sum$coefficients)
    }
    corl_res = do.call(rbind, corl)
    corl_res = as.data.frame(apply(corl_res, 2, as.numeric))
    
    pref = c(rownames(sum$coefficients)[-nrow(sum$coefficients)], "Protein")
    suffix = colnames(sum$coefficients)
    nams = list()
    for(i in suffix) {nams[[i]] = paste0(pref, "_", i)}
    
    colnames(corl_res) = unlist(nams)
    corl_res$adj.p_values = p.adjust(as.numeric(corl_res$`Protein_Pr(>|t|)`), method="fdr")
    corl_res$Gene.Name = colnames(xdf)
    s = which(corl_res$adj.p_values<0.05)
    glm.res[[j]] = corl_res
    glm.res_sign[[j]] = corl_res[s,]
    print(paste("Run", j, "done."))
  }
  
  all_res = list(glm.res, glm.res_sign)
  names(all_res) = c("All_results", "Significant_5%FDR")
  
  return(all_res)
}

##### Proteins adjusted for each other #####
glm.radio_proteome = function(rdf.sel, xdfsel, sel.var, sdfx=NULL, adj.var=NULL, nfold=NULL) {
  cv.fit = list()
  coefs = list()
  predictions = list()
  observations = list()
  cors = list()
  
  for(j in sel.var) {
    outcomes = rdf.sel |> filter(!is.na(!!sym(j)))
    
    traind = xdfsel
    
    if(is.null(adj.var)) {
      traind = traind[which(rownames(traind) %in% rownames(outcomes)),]
    } else {
      traind = cbind(xdfsel, sdfx[,adj.var])
      traind = traind[which(rownames(traind) %in% rownames(outcomes)),]
    }
    
    if(is.null(nfold)) {
      nfold = nrow(traind)
    }
    
    cv.fit[[j]] = cv.glmnet(as.matrix(traind), 
                            outcomes[,j],
                            nfolds=nfold
    )
    plot(cv.fit[[j]])
    
    coefs[[j]] = coef(cv.fit[[j]], s="lambda.min")
    observations[[j]] = outcomes[,j]
    predictions[[j]] = predict(cv.fit[[j]], newx=as.matrix(traind), newy = observations[[j]], 
                               s="lambda.min")
    plot(predictions[[j]], observations[[j]])
    cors[[j]] = cor.test(predictions[[j]], observations[[j]], method="spearman")
  }
  
  all_res = list(cv.fit, coefs, predictions, observations, cors)
  names(all_res) = c("GLM", "coefficients", "predictions", "observations", "correlations")
  
  return(all_res)
}
