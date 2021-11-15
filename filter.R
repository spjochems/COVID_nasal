rm(list = ls())

library('flowCore')
library(ggplot2)
library("gridExtra")

setwd('C:/Users/spjochems/Day5/1.Cytofclean/')


files <- list.files(pattern  = 'fcs')

arcsinh5 <- function(x){
  return(asinh(x/5))  
}

dat <- read.FCS(files[1], transform = F)
labels <- dat@parameters@data

for(j in 1:length(files)){
  dat <- read.FCS(files[j], transform = F)

  
  #####################  
  #DNA gating  
  #####################
  data <- dat@exprs
  
  #DNA selection
  DNA_min = 4.7
  DNA_max = 5.8
  
  #dataframe
  marker1 <- data.frame(data[,c(colnames(dat@exprs) == 'Ir191Di' | colnames(dat@exprs) == 'Time')])
  marker1[,2] <- arcsinh5(marker1[,2]) #asin transformation
  keep1 <- marker1[,2] > DNA_min & marker1[,2] < DNA_max #cells to keep
  keep_perc <- table(keep1)[2] /nrow(marker1)*100 #% of cells to keep
  
  #downsample
  if(nrow(marker1) > 100000){
    marker1 <- marker1[sample(nrow(marker1), 100000), ]
  }
  
  #plot
  colnames(marker1) <- c('Time', 'DNA1')
  a <- ggplot(data = marker1, aes(x = Time, y = DNA1)) + 
    geom_point(alpha = 0.1) + 
    theme_bw() + 
    ggtitle(paste0(files[j], '      ', round(keep_perc, digits = 1), '% Kept')) + 
    geom_hline(yintercept = DNA_max, linetype = 'dashed', colour = 'blue', size = 1.2) + 
    geom_hline(yintercept = DNA_min, linetype = 'dashed', colour = 'blue', size = 1.2)
  a
  
  #subset flowframe on gate
  dat2 <- Subset(dat, keep1)
  
  
  
  
  
  setwd("C:/Users/spjochems/Day5/2.Filter_DNA")
  dir.create('Pics')
  
  #make pictures
  png(file = paste0('Pics/Filter_DNA_', files[j], '.png'), width = 600, height = 600)
  print(a)
  dev.off()
  
  #make files
  write.FCS(dat2, filename = paste0('Filter_', files[j]))

  setwd('../1.Cytofclean/')
  
}
