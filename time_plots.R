library('flowCore')
library(ggplot2)
library("gridExtra")

setwd('C:/Users/spjochems/Day5/2.Filter_DNA/')

files <- list.files(pattern  = 'fcs')

arcsinh5 <- function(x){
  return(asinh(x/5))  
}

dat <- read.FCS(files[1], transform = F)
labels <- dat@parameters@data

for(i in 3:60){
  desc <- labels$desc[i]
  name <- labels$name[i]

  if(length(strsplit(desc, '_')[[1]])>1){
    
    for(j in 1:length(files)){
      
      dat <- read.FCS(files[j], transform = F)
      data <- dat@exprs
      marker <- data.frame(data[,c(colnames(dat@exprs) == name | colnames(dat@exprs) == 'Time')])

      if(nrow(marker) > 20000){

        marker <- marker[sample(nrow(marker), 20000), ]
      }
      
      marker[,2] <- arcsinh5(marker[,2]) 
      marker[,2][marker[,2] > 6] <- 6
      colnames(marker)[2] <- 'Intensity'
      a <- ggplot(data = marker, aes(x = Time, y = Intensity)) + 
           geom_point(alpha = 0.1) + 
           theme_bw() + 
           ylab(desc) + 
          ggtitle(files[j]) + 
           ylim(c(0,6))


      assign(paste0('plot', j), a)
    }

    #create a png file with name, plot all plots and close file
    png(file = paste0('../3.Timeplots/Time_', desc, '.png'), width = 2400, height = 1800)
    grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8,
                 plot9, plot10, plot11, plot12, ncol=4)
    dev.off() 

    rm(list = ls(pattern = 'plot')) #remove the created plots
  }
}
