rm(list = ls())

library('flowCore')
library(ggplot2)
library(reshape2)
library("gridExtra")
library('reshape2')

setwd('C:\\Users\\spjochems\\Day5/4.Compensate/')

files <- list.files(pattern  = 'fcs')


arcsinh5 <- function(x){
  return(asinh(x/5))  
}

dat <- read.FCS(files[1], transform = F)
labels <- dat@parameters@data


yields <- matrix(ncol = 2, nrow = length(files))
rownames(yields) <- files
colnames(yields) <- c('Epithelial', 'Immune')


for(j in 1:length(files)){
  dat <- read.FCS(files[j], transform = F, truncate_max_range = F)




  #####################  
  #cPARP  
  #####################
  data <- dat@exprs
  
  #L/D selection
  cPARP_co = 2
  
  
  #dataframe
  names <- c('Time', 'Nd143Di')
  marker1 <- data.frame(data[,c(colnames(dat@exprs) %in% names)])
  marker1[,2] <- arcsinh5(marker1[,2]) #asin transformation

  
  #keep immune cells, based on cd45 or cd66b and then split
  keep <-  marker1[,2] <  cPARP_co 
  keep_perc <- table(keep)[2] /nrow(marker1)*100 #% of cells to keep
  
  #downsample
  if(nrow(marker1) > 100000){
    marker1 <- marker1[sample(nrow(marker1), 100000), ]
  }
  
  #plot
  colnames(marker1) <- c('Time', 'cPARP')
  a <- ggplot(data = marker1, aes(x = Time, y = cPARP)) + 
    geom_point(alpha = 0.1) + 
    theme_bw() + 
    ggtitle(paste0(files[j], '     ', round(keep_perc, digits = 1), '% Alive')) + 
    geom_hline(yintercept = cPARP_co, linetype = 'dashed', colour = 'blue', size = 1.2) 
  a
  
  
dat1 <- Subset(dat, keep)  

#####################  
#CD45/Epcam gating  
#####################
data <- dat1@exprs

#L/D selection
Epcam_co = 2.1
CD45_co = 2.5


#dataframe
names <- c('Y89Di', 'Gd160Di')
marker1 <- data.frame(data[,c(colnames(dat@exprs) %in% names)])
marker1[,1:2] <- arcsinh5(marker1[,1:2]) #asin transformation

#keep immune cells, based on cd45 or cd66b and then split
keepimm <-  marker1[,1] > CD45_co 
keepCD45_perc <- table(keepimm)[2] /nrow(marker1)*100 #% of cells to keep

keepep <-  marker1[,2] > Epcam_co & marker1[,1] < CD45_co  #cells to keep
keepep_perc <- table(keepep)[2] /nrow(marker1)*100 #% of cells to keep

#downsample
if(nrow(marker1) > 100000){
  marker1 <- marker1[sample(nrow(marker1), 100000), ]
}

#plot
colnames(marker1) <- c('CD45', 'Epcam')
b <- ggplot(data = marker1, aes(x = CD45, y = Epcam)) + 
  geom_point(alpha = 0.1) + 
  theme_bw() + 
  ggtitle(paste0(files[j], '     ', round(keepCD45_perc, digits = 1), '% Immune; ',
                                    round(keepep_perc, digits = 1), '% Epcam+')) + 
  geom_hline(yintercept = Epcam_co, linetype = 'dashed', colour = 'blue', size = 1.2) +
  geom_vline(xintercept = CD45_co, linetype = 'dashed', colour = 'blue', size = 1.2)
b



datCD45 <- Subset(dat1, keepimm)
datEp <- Subset(dat1, keepep)






#####################  
#Remove doublets  
#####################
data <- datCD45@exprs

#L/D selection
CD66b_co = 4
CD3_co = 4
CD14_co = 2

#dataframe
names <- c('Nd144Di', 'Nd148Di', 'Er170Di')
marker1 <- data.frame(data[,c(colnames(dat@exprs) %in% names)])
marker1[,1:3] <- arcsinh5(marker1[,1:3]) #asin transformation

#keep immune cells, based on cd45 or cd66b and then split
singlets1 <-  !(marker1[,1] > CD66b_co &  marker1[,3] > CD3_co)
singlets2 <-  !(marker1[,1] > CD66b_co &  marker1[,2] > CD14_co)
singlets3 <-  !(marker1[,2] > CD14_co &  marker1[,3] > CD3_co)
singlets <- singlets1 & singlets2 & singlets3
keep_perc <- table(singlets)[2] /nrow(marker1)*100 #% of cells to keep

#downsample
if(nrow(marker1) > 100000){
  marker1 <- marker1[sample(nrow(marker1), 100000), ]
}

#plot
colnames(marker1) <- c('CD66b', 'CD14', 'CD3')
c <- ggplot(data = marker1, aes(x = CD66b, y = CD3)) + 
  geom_point(alpha = 0.1) + 
  theme_bw() + 
  ggtitle(paste0(files[j], '     ', round(keep_perc, digits = 1), '% singlets')) + 
  geom_hline(yintercept = CD3_co, linetype = 'dashed', colour = 'blue', size = 1.2) +
  geom_vline(xintercept = CD66b_co, linetype = 'dashed', colour = 'blue', size = 1.2)
c


d <- ggplot(data = marker1, aes(x = CD66b, y = CD14)) + 
  geom_point(alpha = 0.1) + 
  theme_bw() + 
  ggtitle(paste0(files[j], '     ', round(keep_perc, digits = 1), '% singlets')) + 
  geom_hline(yintercept = CD14_co, linetype = 'dashed', colour = 'blue', size = 1.2) +
  geom_vline(xintercept = CD66b_co, linetype = 'dashed', colour = 'blue', size = 1.2)
d


e <- ggplot(data = marker1, aes(x = CD3, y = CD14)) + 
  geom_point(alpha = 0.1) + 
  theme_bw() + 
  ggtitle(paste0(files[j], '     ', round(keep_perc, digits = 1), '% singlets')) + 
  geom_hline(yintercept = CD14_co, linetype = 'dashed', colour = 'blue', size = 1.2) +
  geom_vline(xintercept = CD3_co, linetype = 'dashed', colour = 'blue', size = 1.2)
e




datImm <- Subset(datCD45, singlets)




setwd('C:\\Users\\spjochems\\Dropbox\\LUMC\\Results\\CyTOF\\BEAT-COVID/Day4/5.Split_cPARP/')

#make pictures
png(file = paste0('Filter_', files[j], '.png'), width = 2000, height = 400)
grid.arrange(a,b,c,d,e, ncol=5)
dev.off()

#make files
write.FCS(datImm, filename = paste0('Imm_', files[j]))
write.FCS(datEp, filename = paste0('Ep_', files[j]))

setwd('C:\\Users\\spjochems\\Day5/4.Compensate/')

yields[j,] <- c(keepep_perc, keep_perc)

}


setwd('C:\\Users\\spjochems\\Day5/5.Split_cPARP/')
write.table(yields, 'Distribution_per_file.txt', sep = '\t')
yields <- data.frame(yields)
yields$ID <- rownames(yields)
yields2 <- melt(yields, id = 'ID')
ggplot(yields2, aes(x=variable, y = value)) + 
  geom_point(aes(colour = ID)) + 
  theme_bw() + 
  ylab('Percentage of cells')
ggsave('Yields_per_sample.pdf', width = 8, height = 5)


