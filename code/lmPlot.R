lmPlot <- function(d, o, c, v = FALSE, p = NULL, f = NULL){
  # d = data df
  # c = the col names of the variables I want to include in the linear model
  # o = name of the outcome (y) variable column
  # runs a linear model
  # prints summary of the linear model if v = TRUE
  # returns a plot of the coefficients of the linear models
  # if p supplied, plot only contains variables (column names)
  # that have p in column name (via grep command)
  # if f supplied, sub f chars with blank string in covariate names
  # on plot
  
  # subset variables for model
  y <- d[,o]
  d <- d[,c]
  d <- cbind(d,y)
  
  # run model
  mod <- lm(y ~ ., data = d)
  
  # print mod summary if v == TRUE
  if (isTRUE(v)){
    print(summary(mod))
  }
  
  # make a plot of normalized coefficients with CI
  
  # first get coefficients and CIs
  coefs <- summary(mod)$coefficients
  xvals <- coefs[grep(p, rownames(coefs)),1]
  CIlow <- xvals - (1.96*coefs[grep(p, rownames(coefs)),2])
  CIhigh <- xvals + (1.96*coefs[grep(p, rownames(coefs)),2])
  pval <- coefs[grep(p, rownames(coefs)),4]
  yvals <- gsub(f, "", names(xvals), fixed = T)
  plt.df <- data.frame(yvals = factor(yvals, levels = yvals),
                       xvals = xvals,
                       CIlow = CIlow,
                       CIhigh = CIhigh,
                       pval = pval)
  
  # make plot using ggplot
  plotcols <- rep("black", nrow(plt.df))
  plotcols[plt.df$pval <= 0.05] <- "red"
  p <- ggplot(plt.df, aes(x = xvals, y = yvals, size = -log10(pval))) + 
    geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
    geom_errorbarh(aes(xmax = CIhigh, xmin = CIlow), size = .5, height = 
                     .2, color = "gray50") +
    geom_point(color = plotcols) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    ylab("") +
    xlab("Standardized Coefficient")
    #ggtitle("Feeding method and risk of obesity in cats")
  
  return(p)
  
}