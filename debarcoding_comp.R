rm(list=ls())

library('CATALYST')
library(flowCore)
library('pheatmap')
library(reshape2)
library(ggplot2)
library(cowplot)

setwd('C:\\Users\\spjochems\\Dropbox\\LUMC\\Results\\CyTOF\\BEAT-COVID/Day5/6.Debarcode_filtered/')

arcsinh5 <- function(x){
  return(asinh(x/5))  
}

#load in barcoding file
sample_key <- read.delim("barcode_file.txt", row.names=1, header=T)
colnames(sample_key) <- gsub('X', '', colnames(sample_key))



for(pattern in c('Imm', 'Ep')){

  dir.create(pattern)
  
#merge together all files per cell type
files <- list.files('C:\\Users\\spjochems\\Dropbox\\LUMC\\Results\\CyTOF\\BEAT-COVID/Day5/5.Split_cPARP/standard_better_USETHIS/', pattern  = pattern)
data <- prepData(x = paste0('C:\\Users\\spjochems\\Dropbox\\LUMC\\Results\\CyTOF\\BEAT-COVID/Day5/5.Split_cPARP/standard_better_USETHIS/' ,files))


setwd(pattern)


#do debarcoding
data <- assignPrelim(data, sample_key, verbose=F)
# estimate separation cutoffs
data <- estCutoffs(data)
# view separation cutoff estimates
metadata(data)$sep_cutoffs
plotYields(data)
sce2 <- applyCutoffs(data)
sce3 <- applyCutoffs(data, sep_cutoffs = mean(metadata(data)$sep_cutoffs))
if(pattern == 'Ep'){
  sce3 <- applyCutoffs(data, sep_cutoffs = mean(metadata(data)$sep_cutoffs[1:19])) #remove the PBMC
}

table(sce2$bc_id)
table(sce3$bc_id)

#print assigned events per sample
all <- table(sce2$bc_id)
all.m <- melt(all)
write.table(all.m, paste0('Dynamic_Assigned_events_counts.txt'), row.names = F, quote = F, sep = '\t')
a <- ggplot(all.m, aes(x=Var1, y=value, label=round(value, 1))) + 
  geom_bar(stat='identity', position="stack", aes(fill = Var1)) + 
  geom_text(size = 3, position = position_stack(vjust = 1.1)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_log10() + 
  theme_bw() +
  theme(legend.position = 'none')
#percentages
all <- table(sce2$bc_id)/sum(table(sce2$bc_id))*100
all.m <- melt(all)
write.table(all.m, paste0('Dynamic_Assigned_events_perc.txt'), row.names = F, quote = F, sep = '\t')
b <- ggplot(all.m, aes(x=Var1, y=value, label=round(value, 1))) + 
  geom_bar(stat='identity', position="stack", aes(fill = Var1)) + 
  geom_text(size = 3, position = position_stack(vjust = 1.1)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme_bw() +
  theme(legend.position = 'none')
plot_grid(a,b, ncol=1)
ggsave(paste0('Dynamic_assigned_', pattern, '.pdf'), width=35, height=35, units='cm')


#with one cutoff
all <- table(sce3$bc_id)
all.m <- melt(all)
write.table(all.m, paste0('Static_Assigned_events_counts.txt'), row.names = F, quote = F, sep = '\t')
a <- ggplot(all.m, aes(x=Var1, y=value, label=round(value, 1))) + 
  geom_bar(stat='identity', position="stack", aes(fill = Var1)) + 
  geom_text(size = 3, position = position_stack(vjust = 1.1)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_log10() + 
  theme_bw() +
  theme(legend.position = 'none')
#percentages
all <- table(sce3$bc_id)/sum(table(sce3$bc_id))*100
all.m <- melt(all)
write.table(all.m, paste0('Static_Assigned_events_perc.txt'), row.names = F, quote = F, sep = '\t')
b <- ggplot(all.m, aes(x=Var1, y=value, label=round(value, 1))) + 
  geom_bar(stat='identity', position="stack", aes(fill = Var1)) + 
  geom_text(size = 3, position = position_stack(vjust = 1.1)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme_bw() +
  theme(legend.position = 'none')
plot_grid(a,b, ncol=1, 
          labels = paste0('Cutoff = ', round(mean(metadata(data)$sep_cutoffs), 2)))
ggsave(paste0('Static_assigned_', pattern, '.pdf'), width=35, height=35, units='cm')


try(dev.off())




#plot the expression of the barcodes
normalized <- data@assays@data$exprs[grepl('B2M', rownames(data)),1:5000] #columns match the barcodes


top <- sce3$bc_id[1:5000]
top2 <- top[top>0]
normalized2 <- t(normalized)[top>0,]

ord <- order(top2, decreasing=F)
pdf(paste0('static_heatmap_barcodes_', pattern, '.pdf'), width=6, height=6)
pheatmap(normalized2[ord,], cluster_cols = F, cluster_rows = F)
dev.off()

#write output files
sce2 <- sce2[, sce2$bc_id != 0] #remove unassigned and wrong
fs <- sce2fcs(sce2, split_by = "bc_id")
all(c(fsApply(fs, nrow)) == table(sce2$bc_id))

dir.create('Dynamic')
# get sample identifiers
ids <- fsApply(fs, identifier)
for (id in ids) {
  ff <- fs[[id]]                     # subset 'flowFrame'
  fn <- paste0(pattern, "_Dyn_Sample_", id, '.fcs') # specify output name that includes ID
  fn <- file.path("Dynamic", fn)         # construct output path
  write.FCS(ff, fn)                  # write frame to FCS
}


#general cutoff
sce3 <- sce3[, sce3$bc_id != 0] #remove unassigned and wrong
fs <- sce2fcs(sce3, split_by = "bc_id")

dir.create('Static')
# get sample identifiers
ids <- fsApply(fs, identifier)
for (id in ids) {
  ff <- fs[[id]]                     # subset 'flowFrame'
  fn <- paste0(pattern, "_Sample_", id, '.fcs') # specify output name that includes ID
  fn <- file.path("Static", fn)         # construct output path
  write.FCS(ff, fn)                  # write frame to FCS
}



setwd('../')

}
