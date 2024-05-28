indices_normality <- function(rich, nrow, ncol) {
  
  ### p-value < 0.05 means data failed normality test
  
  par(mfrow = c(nrow, ncol))
  
  for (i in names(rich)) {
    shap <- shapiro.test(rich[, i])
    qqnorm(rich[, i], main = i, sub = shap$p.value)
    qqline(rich[, i])
  }
  
  par(mfrow = c(1, 1))
}

rankabuncomp<-function (x, y = NULL, factor = NULL, return.data = T, specnames = c(1:3), 
                        scale = "abundance", scaledx = F, type = "o", rainbow = T, 
                        legendpos = "topright", xlim = c(1, max1), ylim = c(0, max2), ...) 
{
  groups <- table(y[, factor])
  levels <- names(groups)
  m <- length(groups)
  max1 <- max(diversitycomp(x, y, factor1 = factor, index = "richness", 
                            method = "pooled")[, 2])
  if (scaledx == T) {
    xlim <- c(0, 100)
  }
  max2 <- max.2 <- 0
  for (i in 1:m) {
    if (scale == "abundance") {
      max.2 <- rankabundance(x, y, factor, levels[i])[1, 
                                                      "abundance"]
    }
    if (scale == "logabun") {
      max.2 <- rankabundance(x, y, factor, levels[i])[1, 
                                                      "abundance"]
    }
    if (scale == "proportion") {
      max.2 <- rankabundance(x, y, factor, levels[i])[1, 
                                                      "proportion"]
    }
    if (max.2 > max2) {
      max2 <- max.2
    }
  }
  if (scale == "accumfreq") {
    max2 <- 100
  }
  max2 <- as.numeric(max2)
  if (rainbow == F) {
    if (scale == "logabun" && all.equal(ylim, c(0, max2)) == 
        T) {
      ylim <- c(1, max2)
    }
    rankabunplot(rankabundance(x, y, factor, levels[1]), 
                 scale = scale, scaledx = scaledx, type = type, labels = levels[1], 
                 xlim = xlim, ylim = ylim, pch = 1, specnames = NULL, 
                 ...)
    for (i in 2:m) {
      rankabunplot(rankabundance(x, y, factor, levels[i]), 
                   addit = T, scale = scale, scaledx = scaledx, 
                   type = type, labels = levels[i], pch = i, specnames = NULL, 
                   ...)
    }
    legend(legendpos, legend = levels, pch = c(1:m))
  }
  else {
    grDevices::palette(colorspace::rainbow_hcl(m, c = 90, 
                                               l = 50))
    if (scale == "logabun" && all.equal(ylim, c(0, max2)) == 
        T) {
      ylim <- c(1, max2)
    }
    rankabunplot(rankabundance(x, y, factor, levels[1]), 
                 scale = scale, scaledx = scaledx, type = type, labels = levels[1], 
                 xlim = xlim, ylim = ylim, col = 1, pch = 1, specnames = NULL, 
                 ...)
    for (i in 2:m) {
      rankabunplot(rankabundance(x, y, factor, levels[i]), 
                   addit = T, scale = scale, scaledx = scaledx, 
                   type = type, labels = levels[i], col = i, pch = i, 
                   specnames = NULL, ...)
    }
    legend(legendpos, legend = levels, pch = c(1:m), 
           col = c(1:m))
    grDevices::palette("default")
  }
  if (return.data == T) {
    for (i in 1:m) {
      resulti <- data.frame(rankabundance(x, y, factor, 
                                          levels[i]))
      resulti <- data.frame(Grouping = rep(levels[i], nrow(resulti)), 
                            species = rownames(resulti), labelit = rep(FALSE, 
                                                                       nrow(resulti)), resulti)
      spec.max <- min(max(specnames), nrow(resulti))
      resulti[c(1:spec.max), "labelit"] <- as.logical(1)
      rownames(resulti) <- NULL
      if (i == 1) {
        result <- resulti
      }
      else {
        result <- rbind(result, resulti)
      }
    }
    return(result)
  }
}
