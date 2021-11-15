rm(list=ls())

library('CATALYST')
library(flowCore)
library('pheatmap')
library(reshape2)
library(ggplot2)
library(SingleCellExperiment)
library(ggcyto)      
library(openCyto)     
library(flowWorkspace) 
library(cowplot)

setwd('C:/Users/spjochems/Dropbox/LUMC/Results/CyTOF/BEAT-COVID/Day5/4.Compensate/')

arcsinh5 <- function(x){
  return(asinh(x/5))  
}



#################################################
#Read in beads and make compensation matrix
####################################################3
sm4 <- read.delim('../../Comp/compmatrix_BC_adapted_D1.txt', sep = '\t', header = 1, row.names = 1)


#compensate all files together and split up into pops for debarcoding
files <- list.files('../2.Filter_DNA/', pattern = 'fcs')

data <- prepData(x = paste0('../2.Filter_DNA/' ,files)) #import data



data <- compCytof(data, sm4, method = "nnls", overwrite = F) #do compensation



#set the names of the B2Ms to be unique
rownames(data)[c(3,5,6,7,9,11)] <- gsub('Di', '',paste0('B2M_',channels(data)[c(3,5,6,7,9,11)]))



#make plots to show the effect of compensation
#############################################################
channels_exp <- data.frame('Channel'= channels(data))
channels_exp$metal <- do.call(c, lapply(strsplit(channels_exp[,1], split = "[0-9]"), function(x) x[1]))
channels_exp$number <- substring(channels_exp[,1], regexpr("[0-9]", channels_exp[,1]))
channels_exp$number <- gsub('Di', '', channels_exp$number)
channels_exp$number  <- as.numeric(channels_exp$number )

unused_ch2 <- c("Cd108Di",  "Cd113Di", "I127Di", "Xe131Di", "Cs133Di", "Ba138Di", "Ce140Di", 
                "BCKG190Di", "Ir191Di", "Ir193Di", "Pb208Di")


dir.create('pics')

for(i in 1:length(channels(data))){
  chan <- channels(data)[i]
  
  #only show used channels
  if(!chan %in% unused_ch2){
    
    #based on oxidation or +/- 1 in mass
    num <- c(channels_exp[i,3] - 1, channels_exp[i,3] + 1, channels_exp[i,3] + 16)
    chan1 <- channels_exp[channels_exp[,3] %in% num,1]
    
    #based on impurities
    met <- channels_exp[i,2]
    chan2 <- channels_exp[channels_exp[,2] == met,1] 
    
    #all comparisons of interest
    all_chan <- c(chan, chan1, chan2)
    all_chan <- all_chan[!(duplicated(all_chan) | all_chan %in% unused_ch2)]
    
    #make plots
    
    if(length(all_chan)>1){
      
      pdf(paste0('pics/', chan, '.pdf'), width = 12, height = 6, onefile = T)
      
      
      for(j in 2:length(all_chan)){
        chs <- c(all_chan[1], all_chan[j])
        as <- c("exprs", "compexprs")
        ps <- lapply(as, function(a) 
          plotScatter(data, chs, assay = a))
        p <- plot_grid(plotlist = ps, nrow = 1, labels =c('Uncompensated', 'Compensated'))
        print(p)
      }
      
      dev.off()
    }
    
  }
  
}




##################################################################33
#write to fcs files
fs <- sce2fcs(data, split_by = "sample_id", assay = 'compcounts')

#rm(data) #remove the SCE in case it gives a memory error

ids <- fsApply(fs, identifier)
for (id in ids) {
  ff <- fs[[id]]                     # subset 'flowFrame'
  fn <- sprintf("compensated_%s", id) # specify output name that includes ID
  write.FCS(ff, fn)                  # write frame to FCS
}

