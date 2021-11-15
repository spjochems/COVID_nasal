rm(list = ls())
setwd("C:/Users/spjochems/BatchCorrection/Days1-5/BatchAdjust/")

source("BatchAdjust.R")
library('xlsx')

dir.create('output')
setwd('../DataPreNormalization/')
files <- list.files(pattern = 'fcs')

#make the new names that can be used for batchadjust
#load in data
md <- read.xlsx(paste0("../sample_info_2021_01_27_filtered.xlsx"), 1)

#add the extra SB versus MQ split from first sample
md <- rbind(md, md[1:20,])
md[1:20,1] <- paste0('MQ_', md[1:20,1])
md[81:100,1] <- paste0('SB_', md[81:100,1])
md <- md[c(21:80, 1:20, 81:100),]



md$file_name <- files
md$filename2 <- md$file_name
md$filename2 <- gsub('Imm_Dyn_Sample_', '', md$filename2)
md$key <- '_'
md$key[md$Group == 'REF'] <- '_REF'
md$batch <- 'Run1'
md$batch[grepl('MQ', md$CyTOF_ID)] <- 'Run2'
md$batch[grepl('SB', md$CyTOF_ID)] <- 'Run3'
md$batch[grepl('C', md$CyTOF_ID)] <- 'Run4'
md$batch[grepl('D', md$CyTOF_ID)] <- 'Run5'
md$batch[grepl('E', md$CyTOF_ID)] <- 'Run6'
md$filename3 <- paste0(md$batch, md$key, md$filename2)

#rename the files
file.rename(from = files, to = md$filename3 )

setwd('C:/Users/spjochems/BatchCorrection/Days1-5/BatchAdjust/')

BatchAdjust(
  basedir="C:/Users/spjochems/BatchCorrection/Days1-5/DataPreNormalization/",
  outdir="C:/Users/spjochems/BatchCorrection/Days1-5/BatchAdjust/output",
  channelsFile = "ChannelsToAdjust.txt",
  batchKeyword="Run",
  anchorKeyword = "REF",
  method="99p", #the mormalization method
  transformation=T, #whetehr data should be transfomred before nomralization or not
  addExt=NULL,
  plotDiagnostics=TRUE)

setwd("output/")
file.rename(to = files, from = md$filename3 )
setwd("C:/Users/spjochems/BatchCorrection/Days1-5/DataPreNormalization/")
file.rename(to = files, from = md$filename3 )
