##### Missing values #####
##### Plot missing values in a matrix #####
ggplot_missing <- function(matrix){
  matrix %>% 
    is.na %>%
    melt %>%
    ggplot(data = .,
           aes(x = as.factor(Var2),
               y = Var1)) +
    geom_raster(aes(fill = value)) + scale_fill_manual(name = "",
                                                       labels = c("Present","Missing"), values=c("black", "white")) +
    scale_y_discrete(labels=NULL) +
    scale_x_discrete(labels=NULL) +
    theme_minimal() + 
    theme(axis.text.x  = element_text(angle=90, hjust=1)) + 
    labs(x = "Columns / Samples",
         y = "Rows / Proteins")
}

##### Calculate % of missing values in a matrix ##### 
miss_func = function(matrix, what=1) {
  apply(matrix, what, \(x) length(which(is.na(x)))/length(x)*100)
}