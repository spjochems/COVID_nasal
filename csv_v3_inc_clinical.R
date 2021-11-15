#this is with 3 granulocyte populations and with the myeloid as tsne


rm(list=ls())

library(reshape2)
library(ggplot2)
library(pheatmap)
library(xlsx)
library(robustrank)
library(stringr)
library(mixOmics)
library(Hmisc)
library(plyr)
library(dplyr)
library(lmerTest)
library(lme4)
library(emmeans)
library(RColorBrewer)
library("factoextra")
library(viridis)
library(cowplot)

setwd('C:/Users/spjochems/Dropbox/LUMC/Results/CyTOF/BEAT-COVID/Analysis/2021_02_05_D1toD5/csv v3/')


info <- read.delim('../sample_info_2021_02_05.txt')
info2 <- info[info$Group == 'BEAT-COVID' & info$EpithelialD>99 & info$ImmuneD>99,]

files <- list.files('../csv v3/', pattern  = 'csv')

###############################################################
#get names of samples from one file
############################################################
data <- read.csv(files[1])
names <- info2$CyTOF_ID

#initialize an empty matrix and name rows and columns
percentages <- matrix(ncol = nrow(info2), nrow = length(files))
rownames(percentages) <- files
rownames(percentages) <- gsub('.csv', '',rownames(percentages))
rownames(percentages) <- gsub('Nasal', '', rownames(percentages))
rownames(percentages) <- gsub(' ', '_', rownames(percentages))
colnames(percentages) <- names
names2 <- data.frame(names)

#load in counts per population
for(i in 1:length(files)){
  data <- read.csv(files[i])
  data$Filename <- gsub('.fcs', '', data$Filename)
  data$Filename <- str_sub(data$Filename, -3, -1)
  data <- merge(names2, data, by.x = 'names', by.y='Filename', all.x=T)
  percentages[i,] <- data[,4]
}

percentages[is.na(percentages)] <- 0


#make into percentage of immune data
counts <- percentages







for(i in 1:ncol(percentages)){
  percentages[,i] <- percentages[,i]/sum(percentages[,i])*100
}
#add in metadata
perc <- data.frame(t(percentages))
perc2 <- merge(info, perc, by.x='CyTOF_ID', by.y=0)

#remove volunteers with less than 100 cells
perc2 <- perc2[perc2$EpithelialD>99 & perc2$ImmuneD>99,]

perc2$Granulocytes <- perc2$CD16Dim_Neutrophils +
  perc2$CD16Hi_Neutrophils + perc2$CD16Neg_Granulocytes + 
  perc2$CD69_Neutrophils + perc2$Mast_cells
perc2$CD4T <- perc2$CD4T_CD161_CD4_TRM + perc2$CD4T_TEM +
  perc2$CD4T_Naive + perc2$CD4T_Treg + perc2$CD4T_CD4_TRM 
perc2$CD8T <- perc2$CD8T_CD8T_EM + perc2$CD8T_CD8T_EMRA + 
  perc2$CD8T_CD161Hi_CD8_T +
  perc2$CD8T_CD8T_Naive + perc2$CD8T_CD8T_TRM 
perc2$Monocytes <- perc2$Myeloid_CD163_Mono + perc2$Myeloid_CD206_CD163_Mono +
  perc2$Myeloid_Classical_Mono 
perc2$mDC <- perc2$Myeloid_CD206_DC + perc2$Myeloid_CD163_DC
perc2$NK <- perc2$ILC_ILC + perc2$ILC_NK_CD11c +
  perc2$ILC_NK_CD69


colnames(perc2) <- gsub('CD4T_DP', 'DP_T', colnames(perc2))
colnames(perc2) <- gsub('CD4T_CD161_CD4_TRM', 'CD4T_CD161_TRM', colnames(perc2))
colnames(perc2) <- gsub('CD4T_TEM', 'CD4T_EM', colnames(perc2))
colnames(perc2) <- gsub('CD4T_Naive', 'CD4T_Naive', colnames(perc2))
colnames(perc2) <- gsub('CD4T_Treg', 'CD4T_Treg', colnames(perc2))
colnames(perc2) <- gsub('CD4T_CD4_TRM', 'CD4T_TRM', colnames(perc2))
colnames(perc2) <- gsub('CD8T_CD161Hi_CD8_T', 'CD8T_CD161Hi', colnames(perc2))
colnames(perc2) <- gsub('CD8T_CD8T_EM', 'CD8T_EM', colnames(perc2))
colnames(perc2) <- gsub('CD8T_CD8T_EMRA', 'CD8T_EMRA', colnames(perc2))
colnames(perc2) <- gsub('CD8T_CD8T_Naive', 'CD8T_Naive', colnames(perc2))
colnames(perc2) <- gsub('CD8T_CD8T_TRM', 'CD8T_TRM', colnames(perc2))
colnames(perc2) <- gsub('CD8T_DN_T', 'DN_T', colnames(perc2))
colnames(perc2) <- gsub('ILC_ILC', 'NK_ILC', colnames(perc2))
colnames(perc2) <- gsub('ILC_NK_CD69', 'NK_CD69', colnames(perc2))
colnames(perc2) <- gsub('ILC_NK_CD11c', 'NK_CD11c', colnames(perc2))
colnames(perc2) <- gsub('Myeloid_CD163_DC', 'DC', colnames(perc2))
colnames(perc2) <- gsub('Myeloid_CD206_DC', 'DC_CD206', colnames(perc2))
colnames(perc2) <- gsub('Myeloid_CD163_Mono', 'Mono_CD163', colnames(perc2))
colnames(perc2) <- gsub('Myeloid_CD206_CD163_Mono', 'Mono_CD163_CD206', colnames(perc2))
colnames(perc2) <- gsub('Myeloid_Classical_Mono', 'Mono', colnames(perc2))
colnames(perc2) <- gsub('CD16Dim_Neutrophils', 'Gran_CD16Dim', colnames(perc2))
colnames(perc2) <- gsub('CD16Hi_Neutrophils', 'Gran_CD16Hi', colnames(perc2))
colnames(perc2) <- gsub('CD16Neg_Granulocytes', 'Gran_CD16Neg', colnames(perc2))
colnames(perc2) <- gsub('CD69_Neutrophils', 'Gran_CD69', colnames(perc2))



write.table(perc2, 'Percentages.txt', sep = '\t')


perc2$label <- paste(perc2$Donor, perc2$Info, sep = '; ')
perc3 <- melt(perc2, id = colnames(perc2)[c(1:ncol(info),ncol(perc2))])

#plot
perc5 <- perc3[perc3$variable %in% c('Granulocytes', 'Monocytes', 'CD4T', 
                                     'CD8T', 'mDC', 'pDC', 
                                     'B_cells',  'NK'),]

perc7 <- perc5[perc5$Sample_Number != 'Acute2' & 
                 perc5$Sample_Number != 'Acute3' & 
                 perc5$Sample_Number != 'Acute4' & 
                 perc5$Sample_Number != 'HD2',]
#perc7 <- perc5
perc7$Sample <- paste0(perc7$Donor, '_', perc7$Sample_Number)



ggplot(perc7, aes(x = Sample, y = value)) + 
  geom_bar(stat = 'identity', aes(fill = variable)) + 
  facet_grid(.~Timepoint, scales = 'free_x', space= 'free') + 
  coord_cartesian(ylim = c(0,100), expand = F) + 
  theme_bw() +
  scale_fill_brewer(palette = 'Set1')
ggsave('Percentages_bars_1timepoint.pdf', width = 25, height = 10, units ='cm', useDingbats = F)


info3 <- read.delim('../info_samples_2021_03_02.txt')
perc77 <- merge(info3, perc7, by = 'CyTOF_ID')
perc77$Group2 <- factor(perc77$Group2, 
                     levels = c('Acute', 'ERS', 'Convalescent', 'HD'))
ggplot(perc77, aes(x = Sample, y = value)) + 
  geom_bar(stat = 'identity', aes(fill = variable)) + 
  facet_grid(.~Group2, scales = 'free_x', space= 'free') + 
  coord_cartesian(ylim = c(0,100), expand = F) + 
  theme_bw() +
  scale_fill_brewer(palette = 'Set1')
ggsave('Percentages_bars_1timepoint_4g.pdf', width = 25, height = 10, units ='cm', useDingbats = F)



perc6 <- perc2[perc2$Sample_Number != 'Acute2' & 
                 perc2$Sample_Number != 'Acute3' & 
                 perc2$Sample_Number != 'Acute4' & 
                 perc2$Sample_Number != 'HD2',]

perc66 <- merge(info3, perc2, by = 'CyTOF_ID') 

ddply(perc66[,
             c(8,which(colnames(perc66) == 'Monocytes'))], .(Group2), summary)




#make ratio to epithelial
Ratio <- data.frame(t(counts))
Ratio2 <- merge(info, Ratio, by.x='CyTOF_ID', by.y=0)

#remove volunteers with less than 100 cells
Ratio2 <- Ratio2[Ratio2$EpithelialD>99 & Ratio2$ImmuneD>99,]
#remove beat 040



Ratio2$Granulocytes <- Ratio2$CD16Dim_Neutrophils +
  Ratio2$CD16Hi_Neutrophils + Ratio2$CD16Neg_Granulocytes + 
  Ratio2$CD69_Neutrophils + Ratio2$Mast_cells
Ratio2$CD4T <- Ratio2$CD4T_CD161_CD4_TRM + Ratio2$CD4T_TEM +
  Ratio2$CD4T_Naive + Ratio2$CD4T_Treg + Ratio2$CD4T_CD4_TRM 
Ratio2$CD8T <- Ratio2$CD8T_CD8T_EM + Ratio2$CD8T_CD8T_EMRA + Ratio2$CD8T_CD161Hi_CD8_T +
  Ratio2$CD8T_CD8T_Naive + Ratio2$CD8T_CD8T_TRM 
Ratio2$Monocytes <- Ratio2$Myeloid_CD163_Mono + Ratio2$Myeloid_CD206_CD163_Mono +
  Ratio2$Myeloid_Classical_Mono 
Ratio2$mDC <- Ratio2$Myeloid_CD206_DC + Ratio2$Myeloid_CD163_DC
Ratio2$NK <- Ratio2$ILC_ILC + Ratio2$ILC_NK_CD11c +
  Ratio2$ILC_NK_CD69


colnames(Ratio2) <- gsub('CD4T_DP', 'DP_T', colnames(Ratio2))
colnames(Ratio2) <- gsub('CD4T_CD161_CD4_TRM', 'CD4T_CD161_TRM', colnames(Ratio2))
colnames(Ratio2) <- gsub('CD4T_TEM', 'CD4T_EM', colnames(Ratio2))
colnames(Ratio2) <- gsub('CD4T_Naive', 'CD4T_Naive', colnames(Ratio2))
colnames(Ratio2) <- gsub('CD4T_Treg', 'CD4T_Treg', colnames(Ratio2))
colnames(Ratio2) <- gsub('CD4T_CD4_TRM', 'CD4T_TRM', colnames(Ratio2))
colnames(Ratio2) <- gsub('CD8T_CD161Hi_CD8_T', 'CD8T_CD161Hi', colnames(Ratio2))
colnames(Ratio2) <- gsub('CD8T_CD8T_EM', 'CD8T_EM', colnames(Ratio2))
colnames(Ratio2) <- gsub('CD8T_CD8T_EMRA', 'CD8T_EMRA', colnames(Ratio2))
colnames(Ratio2) <- gsub('CD8T_CD8T_Naive', 'CD8T_Naive', colnames(Ratio2))
colnames(Ratio2) <- gsub('CD8T_CD8T_TRM', 'CD8T_TRM', colnames(Ratio2))
colnames(Ratio2) <- gsub('CD8T_DN_T', 'DN_T', colnames(Ratio2))
colnames(Ratio2) <- gsub('ILC_ILC', 'NK_ILC', colnames(Ratio2))
colnames(Ratio2) <- gsub('ILC_NK_CD69', 'NK_CD69', colnames(Ratio2))
colnames(Ratio2) <- gsub('ILC_NK_CD11c', 'NK_CD11c', colnames(Ratio2))
colnames(Ratio2) <- gsub('Myeloid_CD163_DC', 'DC', colnames(Ratio2))
colnames(Ratio2) <- gsub('Myeloid_CD206_DC', 'DC_CD206', colnames(Ratio2))
colnames(Ratio2) <- gsub('Myeloid_CD163_Mono', 'Mono_CD163', colnames(Ratio2))
colnames(Ratio2) <- gsub('Myeloid_CD206_CD163_Mono', 'Mono_CD163_CD206', colnames(Ratio2))
colnames(Ratio2) <- gsub('Myeloid_Classical_Mono', 'Mono', colnames(Ratio2))
colnames(Ratio2) <- gsub('CD16Dim_Neutrophils', 'Gran_CD16Dim', colnames(Ratio2))
colnames(Ratio2) <- gsub('CD16Hi_Neutrophils', 'Gran_CD16Hi', colnames(Ratio2))
colnames(Ratio2) <- gsub('CD16Neg_Granulocytes', 'Gran_CD16Neg', colnames(Ratio2))
colnames(Ratio2) <- gsub('CD69_Neutrophils', 'Gran_CD69', colnames(Ratio2))




#normalize to epithelial cells
for(i in (ncol(info)+1):ncol(Ratio2)){
  Ratio2[,i] <- Ratio2[,i]/Ratio2[,ncol(info)]
}

Ratio22 <- Ratio2

#do log transformation
for(i in (ncol(info)+1):ncol(Ratio2)){
  m <- min(Ratio2[,i])
  if(m == 0){
    m1 <- min(Ratio2[,i][Ratio2[,i]>0])
    Ratio2[,i] <- Ratio2[,i] + m1/2
  }
  Ratio2[,i] <- log10(Ratio2[,i])
}


#add labels and metl
write.table(Ratio2, 'Ratio.txt', sep = '\t')

Ratio2 <- Ratio2[Ratio2$Sample_Number != 'HD2',]

load("R:/Para-CIH/CIH Group member folders/Simon/Results/BEAT-COVID database/HermelijnSmits.Rdata")


virus <- rbind(BC_Labresults_nasalcytof[,c(1,2,29,30)], BC_Labresults_wbs[,c(1,2,24,25)])
virus <- virus[!is.na(virus$viral_posneg),]  
virus <- virus[!duplicated(virus$bc_id),]  

virus2 <- merge(BC_Info, virus, by = 'bc_id', all.y=T)  

#check if duplicates
sev <- BC_Severity_Score_Adhoc
Ratio2$ID <- paste(Ratio2$Donor, Ratio2$Date, sep = '_')

dups <- sev$beat_id[grepl('_' , sev$beat_id)]
dups2 <- sev[sev$beat_id %in% dups,]
dups3 <- gsub('_1', '', dups)
Ratio2$Donor[Ratio2$Donor %in% dups3]
#beat 113 and beat 119 are admitted twice, get the right info


BC_Severity_Score_Adhoc$beat_id <- gsub('_2', '', BC_Severity_Score_Adhoc$beat_id)
BC_Severity_Score_Adhoc$beat_id <- gsub('_1', '', BC_Severity_Score_Adhoc$beat_id)
BC_Severity_Score_Adhoc$ID <- paste(BC_Severity_Score_Adhoc$beat_id, 
                                    BC_Severity_Score_Adhoc$day_date, sep = '_')

Ratio4 <- merge(BC_Severity_Score_Adhoc, Ratio2, by = 'ID', all.y=T)

#add patient ID for the controls and convalescents
BC_Info$beat_id <- gsub('H', 'H0', BC_Info$beat_id)
for(i in 1:nrow(Ratio4)){
  if(is.na(Ratio4$patient_id[i])){
    bi <- Ratio4$Donor[i]
    Ratio4$patient_id[i] <- BC_Info[grepl(bi, BC_Info$beat_id),]$patient_id[1]
    Ratio4$beat_id[i] <- BC_Info[grepl(bi, BC_Info$beat_id),]$beat_id[1]
  }
}
BC_Info$beat_id <- gsub('H0', 'H', BC_Info$beat_id)
Ratio4$beat_id <- gsub('H0', 'H', Ratio4$beat_id)


Ratio5 <- merge(C19_BASELINE_EA[!duplicated(C19_BASELINE_EA$patient_id),], Ratio4, by = 'patient_id', all.y=T)

#fix mistakes
Ratio5[Ratio5$beat_id == 'BEAT-139',]$onset_date <- '2020-10-29'

Ratio5[Ratio5$beat_id == 'BEAT-116',]$admission_date <- '2020-11-01'
Ratio5[Ratio5$beat_id == 'BEAT-116',]$Outcome_dt <- '2020-11-20'
Ratio5$hospitalisation_total[Ratio5$Donor == 'BEAT-116'] <- '20'

Ratio5$Outcome_dt[Ratio5$Donor == 'BEAT-113'] <- '2020-11-17'
Ratio5$hospitalisation_total[Ratio5$Donor == 'BEAT-113'] <- '23'

p <- BC_Severity_Score_Adhoc[BC_Severity_Score_Adhoc$beat_id == 'BEAT-119',]$patient_id[1]
pp <- C19_BASELINE_EA[C19_BASELINE_EA$patient_id == p,]
Ratio5[Ratio5$beat_id == 'BEAT-119',]$Outcome_dt <- '2020-11-11'
Ratio5[Ratio5$beat_id == 'BEAT-119',]$hospitalisation_total <- '14'
Ratio5[Ratio5$beat_id == 'BEAT-113',]$Outcome_dt <- '2020-11-18'
Ratio5[Ratio5$beat_id == 'BEAT-030',]$Outcome_dt <- '2020-05-12'
Ratio5[Ratio5$beat_id == 'BEAT-030',]$hospitalisation_total <- '12'
Ratio5[Ratio5$beat_id == 'BEAT-113',]$hospitalisation_total <- '23'

#find ICUs and add the yes/no + dates
Ratio5$ICU <- NA
Ratio5$ICU_start <- NA
Ratio5$ICU_stop <- NA

#add the ICU data
persons <- Ratio5$beat_id[!duplicated(Ratio5$beat_id)]
for(person in persons){
  
  a <- BC_Severity_Score_Adhoc[BC_Severity_Score_Adhoc$beat_id == person,]
  if(nrow(a) > 0){
    if(sum(a$icu_yn) == 0){Ratio5$ICU[Ratio5$beat_id == person] <- 'No'
    } else {Ratio5$ICU[Ratio5$beat_id == person] <- 'Yes'
    Ratio5$ICU_start[Ratio5$beat_id == person] <- as.character(a$day_date[a$icu_yn == 1][1])
    Ratio5$ICU_stop[Ratio5$beat_id == person] <- as.character(tail(a$day_date[a$icu_yn == 1], n=1))
    }
  }
}

#add the medication data
Ratio5$Medication <- NA

persons <- Ratio5$patient_id[!duplicated(Ratio5$patient_id)]
for(person in persons){
  
  a <- DIG_HixMedicatietoedieningDeellijst[DIG_HixMedicatietoedieningDeellijst$patient_id == person,]
  if(nrow(a) > 0){
    Ratio5$Medication[Ratio5$patient_id == person] <- paste(names(table(a$Medication_code5_ATC_display)), collapse = '|')
  }
}


#add data
Ratio5$days_since_admission <- as.Date(Ratio5$Date) - as.Date(Ratio5$admission_date)
Ratio5$days_since_onset <- as.Date(Ratio5$Date) - as.Date(Ratio5$onset_date)
Ratio5$days_since_discharge <- as.Date(Ratio5$Date) - as.Date(Ratio5$Outcome_dt)





#reorder to have last
Ratio6 <- Ratio32 <-Ratio5[,c(1:73, 108:114, 74:107)]


conv <- Ratio6[Ratio6$Timepoint == 'Convalescent',]
write.xlsx(conv[,c(1,62)], 'convalescent patients.xlsx', row.names = F)


Ratio6 <- Ratio6[,c(59:65, 55:57,68, 72:74,78:114)]
write.table(Ratio6[,1:ncol(Ratio6)], 'Ratio_annotated2.txt', sep = '\t', row.names = F)


perc7 <- perc5[perc5$Sample_Number != 'Acute2' & 
                 perc5$Sample_Number != 'Acute3' & 
                 perc5$Sample_Number != 'Acute4' & 
                 perc5$Sample_Number != 'HD2',]
#perc7 <- perc5
perc7 <- merge(perc7, Ratio6[,c(1,15)], by = 'CyTOF_ID')
perc7$Sample <- paste0(perc7$Donor, '_', 
                       perc7$Sample_Number)

orde <- perc7$Sample[order(perc7$days_since_admission)]
perc7$Sample <- factor(perc7$Sample, 
                       levels = orde[!duplicated(orde)])


ggplot(perc7, aes(x = Sample, y = value)) + 
  geom_bar(stat = 'identity', aes(fill = variable)) + 
  facet_grid(.~Timepoint, scales = 'free_x', space= 'free') + 
  coord_cartesian(ylim = c(0,100), expand = F) + 
  theme_bw() +
  scale_fill_brewer(palette = 'Set1')
ggsave('Percentages_bars_1timepoint.pdf', width = 25, height = 10, units ='cm', useDingbats = F)



########################################################################
#Ratio plots
#########################################################################
Ratio6 <- read.delim('Ratio_annotated2.txt')


cor <- Ratio6
cor$Acute <- 'No'
cor$Acute[cor$Timepoint == 'Acute'] <- 'Yes'
cor$Group2 <- cor$Timepoint
cor$Group2[as.numeric(cor$days_since_admission) > 11] <- 'ERS'
cor$Group2[cor$Timepoint == 'Convalescent'] <- 'Convalescent'

ac <- cor[cor$Acute == 'Yes',]
ac$ImmuneRank <- rank(ac$ImmuneD)
ac$EpRank <- rank(ac$EpithelialD)
non <- cor[cor$Acute == 'No',]
non$ImmuneRank <- rank(as.numeric(non$ImmuneD))
non$EpRank <- rank(non$EpithelialD)


a<-ggplot(ac, aes(x=EpRank, y=ImmuneRank)) + 
  geom_point(aes(colour = Group2)) + 
  scale_colour_manual(values = c('#e34a33', '#fdbb84' )) +
  stat_smooth(method = 'lm') + 
  theme_bw() + ggtitle('Hospital patients') + 
  theme(aspect.ratio = 1, legend.position = 'none') 

b<-ggplot(non, aes(x=EpRank, y=ImmuneRank)) + 
  geom_point(aes(colour = Group2)) + 
  scale_colour_manual(values = c('#fa9fb5',  '#2b8cbe')) +
  stat_smooth(method = 'lm') + 
  theme_bw() + ggtitle('Recovered and HD')+ 
  theme(aspect.ratio = 1, legend.position = 'none')


plot_grid(a, b)
ggsave('Correlations Ep_Imm.pdf', width = 4, height = 2)

cor.test(non$EpithelialD, non$ImmuneD, method = 'spearman')
cor.test(ac$EpithelialD, ac$ImmuneD, method = 'spearman')



acute <- cor[cor$Group2 == 'Acute',]
acute$ImmuneRank <- rank(acute$ImmuneD)
acute$EpRank <- rank(acute$EpithelialD)
ers <- cor[cor$Group2 == 'ERS',]
ers$ImmuneRank <- rank(ers$ImmuneD)
ers$EpRank <- rank(ers$EpithelialD)
conv <- cor[cor$Group2 == 'Convalescent',]
conv$ImmuneRank <- rank(conv$ImmuneD)
conv$EpRank <- rank(conv$EpithelialD)
health <- cor[cor$Group2 == 'HD',]
health$ImmuneRank <- rank(health$ImmuneD)
health$EpRank <- rank(health$EpithelialD)

a<-ggplot(acute, aes(x=EpRank, y=ImmuneRank)) + 
  geom_point() + 
  stat_smooth(method = 'lm') + 
  theme_bw() + ggtitle('Acute')

b<-ggplot(ers, aes(x=EpRank, y=ImmuneRank)) + 
  geom_point() + 
  stat_smooth(method = 'lm') + 
  theme_bw() + ggtitle('ERS')

c<-ggplot(conv, aes(x=EpRank, y=ImmuneRank)) + 
  geom_point() + 
  stat_smooth(method = 'lm') + 
  theme_bw() + ggtitle('REcovered')

d<-ggplot(health, aes(x=EpRank, y=ImmuneRank)) + 
  geom_point() + 
  stat_smooth(method = 'lm') + 
  theme_bw() + ggtitle('HD')


plot_grid(a, b,c, d)



cor.test(acute$EpithelialD, acute$ImmuneD, method = 'spearman')
cor.test(ers$EpithelialD, ers$ImmuneD, method = 'spearman')
cor.test(conv$EpithelialD, conv$ImmuneD, method = 'spearman')
cor.test(health$EpithelialD, health$ImmuneD, method = 'spearman')



non <- cor[cor$Acute == 'No',]
non$ImmuneRank <- rank(as.numeric(non$ImmuneD))
non$EpRank <- rank(non$EpithelialD)






ggplot(cor, aes(x=log10(EpithelialD), y=log10(ImmuneD))) + 
  geom_point() + 
  facet_grid(~Timepoint) + 
  stat_smooth(method = 'lm') + 
  theme_bw()


Ratio6.melt <- melt(Ratio6, id = colnames(Ratio6)[c(1:17)])



#plots one timepoint
Ratio44 <- Ratio6.melt[Ratio6.melt$Sample_Number != 'Acute2' & 
                   Ratio6.melt$Sample_Number != 'Acute3' & 
                   Ratio6.melt$Sample_Number != 'Acute4' & 
                   Ratio6.melt$Sample_Number != 'HD2',]



Ratio55 <- Ratio44[Ratio44$variable %in% c('Granulocytes', 'Monocytes', 'CD4T', 
                                        'CD8T', 'mDC', 'pDC', 
                                        'B_cells',  'NK'),]
ggplot(Ratio55, aes(x=Timepoint, y=value,  group = Donor, colour = Timepoint)) + 
  geom_boxplot(outlier.colour = NA, aes(group = Timepoint)) + 
  geom_point() + 
  geom_line(colour = 'grey') +
  facet_wrap(~variable, scales = 'free_y', ncol = 4) + 
  theme_bw() +
  ylab('Log10 Ratio to epithelial cells') +
  theme(legend.position = 'none') +
  scale_colour_manual(values = c('#e34a33', '#fa9fb5', '#2b8cbe')) +
  theme(strip.background = element_rect(color=NA, fill=NA))
ggsave('Ratio to Ep_lineages_1acutepoint_colours.pdf', width = 13, height = 12, units = 'cm', useDingbats = F)



Ratio6.melt$days_since_admission[Ratio6.melt$Timepoint == 'HD'] <- 125
Ratio6.melt$patient <- 'yes'
Ratio6.melt$patient[Ratio6$Timepoint == 'HD'] <- 'no'

ggplot(Ratio6.melt, aes(x=days_since_admission, y=value)) + 
  geom_point(aes(colour = Timepoint)) + 
  facet_wrap(~variable, scales = 'free_y') + 
  stat_smooth(method = 'lm',aes(group = Timepoint, colour= Timepoint)) +
  theme_bw() +
  ylab('Log10 Ratio to epithelial cells') +
scale_colour_manual(values = c('#e34a33', '#fa9fb5', '#2b8cbe')) 
#ggsave('Ratio to Ep_admission_incHD_pgroup.pdf', width = 20, height = 20)


Ratio6.melt$lines <- paste(Ratio6$Donor, Ratio6$Timepoint, sep = '_')
ggplot(Ratio6.melt, aes(x=days_since_admission, y=value)) + 
  geom_point(aes(colour = Timepoint)) + 
  geom_line(aes(group = Donor), colour = 'grey') + 
  facet_wrap(~variable, scales = 'free_y') + 
  stat_smooth(colour = 'black') +
  theme_bw() +
  ylab('Log10 Ratio to epithelial cells') +
  scale_colour_manual(values = c('#e34a33', '#fa9fb5', '#2b8cbe')) 
ggsave('Ratio to Ep_admission_incHD_plines.pdf', width = 10, height = 10)

ggplot(Ratio6.melt, aes(x=days_since_admission, y=value)) + 
  geom_point(aes(colour = Timepoint)) + 
  facet_wrap(~variable, scales = 'free_y') + 
  stat_smooth(colour = 'black', aes(group = patient)) +
  theme_bw() +
  ylab('Log10 Ratio to epithelial cells') +
  theme(strip.background = element_rect(color=NA, fill=NA)) +
  scale_colour_manual(values = c('#e34a33', '#fa9fb5', '#2b8cbe')) + 
  geom_vline(xintercept = 120, linetype = 'dashed')
ggsave('Ratio to Ep_admission_incHD_nolines.pdf', width = 10, height = 10)


#plot per donor
ggplot(Ratio6.melt, aes(x=Sample_Number, y=value, group = Donor, colour = Timepoint)) + 
  geom_boxplot(outlier.colour = NA, aes(group = Timepoint)) + 
  geom_point() + 
  geom_line(colour = 'grey') +
  facet_wrap(~variable, scales = 'free_y') + 
  theme_bw() +
  ylab('Log10 Ratio to epithelial cells') +
  theme(legend.position = 'none') +
  scale_colour_manual(values = c('#e34a33', '#fa9fb5', '#2b8cbe')) +
  theme(strip.background = element_rect(color=NA, fill=NA))
#ggsave('ratio to ep_acutepoint_colours_incGran.pdf', width = 15, height = 15, useDingbats = F)


a<- ddply(Ratio6[,
           c(52,which(colnames(Ratio6) == 'Gran_CD16Hi'))], .(Group2), summary)
a$number <- do.call(c, lapply(strsplit(a[,3], split = ':'), function(x) x[2]))
a$number <- 10^as.numeric(gsub(' ', '', a$number))
a[3,4] / a[21,4] #acute fc
a[15,4] / a[21,4] #ERS fc
a[9,4] / a[21,4] #conv fc

a<- ddply(Ratio6[,
                 c(52,which(colnames(Ratio6) == 'Gran_CD16Dim'))], .(Group2), summary)
a$number <- do.call(c, lapply(strsplit(a[,3], split = ':'), function(x) x[2]))
a$number <- 10^as.numeric(gsub(' ', '', a$number))
a[3,4] / a[21,4] #acute fc
a[15,4] / a[21,4] #ERS fc
a[9,4] / a[21,4] #conv fc


a<- ddply(Ratio6[,
                 c(52,which(colnames(Ratio6) == 'Gran_CD16Neg'))], .(Group2), summary)
a$number <- do.call(c, lapply(strsplit(a[,3], split = ':'), function(x) x[2]))
a$number <- 10^as.numeric(gsub(' ', '', a$number))
a[3,4] / a[21,4] #acute fc
a[15,4] / a[21,4] #ERS fc
a[9,4] / a[21,4] #conv fc


a<- ddply(Ratio6[,
                 c(52,which(colnames(Ratio6) == 'CD4T_EM'))], .(Group2), summary)
a$number <- do.call(c, lapply(strsplit(a[,3], split = ':'), function(x) x[2]))
a$number <- 10^as.numeric(gsub(' ', '', a$number))
a[3,4] / a[21,4] #acute fc
a[15,4] / a[21,4] #ERS fc
a[9,4] / a[21,4] #conv fc


a<- ddply(Ratio6[,
                 c(52,which(colnames(Ratio6) == 'CD8T_EMRA'))], .(Group2), summary)
a$number <- do.call(c, lapply(strsplit(a[,3], split = ':'), function(x) x[2]))
a$number <- 10^as.numeric(gsub(' ', '', a$number))
a[3,4] / a[21,4] #acute fc
a[15,4] / a[21,4] #ERS fc
a[9,4] / a[21,4] #conv fc




#do stats and median FC, make the matrix with all timepoints #do with lmer and emmeans
#ratio 6 is with the repeat samples
stats2 <- matrix(nrow = 34, ncol = 3)
rownames(stats2)  <- colnames(Ratio6)[18:51]
colnames(stats2) <- c('Ac_Conv', 'Ac_HD', 'Conv_HD')
for(i in 18:ncol(Ratio6)){
  mixed.lmer <- lmerTest::lmer(Ratio6[,i] ~ Timepoint + (1|Donor), 
                            data = Ratio6)

  sig <- emmeans(mixed.lmer, pairwise ~ Timepoint, adjust = "tukey")
  stats2[i-17,] <- summary(sig)[[2]][,6]
}


stats2 <- data.frame(stats2)


stats2.subsets <- stats2[!rownames(stats2) %in% c('Granulocytes', 'Monocytes', 'CD4T', 
                                                'CD8T', 'mDC', 'pDC', 
                                                'B_cells',  'NK'),]
stats2.subsets$AC_Conv_adj <- p.adjust(stats2.subsets$Ac_Conv, method = 'BH')
stats2.subsets$AC_HD_adj <- p.adjust(stats2.subsets$Ac_HD, method = 'BH')
stats2.subsets$Conv_HD_adj <- p.adjust(stats2.subsets$Conv_HD, method = 'BH')

stats2.lineages <- stats2[rownames(stats2) %in% c('Granulocytes', 'Monocytes', 'CD4T', 
                                                  'CD8T', 'mDC', 'pDC', 
                                                  'B_cells',  'NK'),]
stats2.lineages$AC_Conv_adj <- p.adjust(stats2.lineages$Ac_Conv, method = 'BH')
stats2.lineages$AC_HD_adj <- p.adjust(stats2.lineages$Ac_HD, method = 'BH')
stats2.lineages$Conv_HD_adj <- p.adjust(stats2.lineages$Conv_HD, method = 'BH')

stats3 <- rbind(stats2.lineages, stats2.subsets)
write.xlsx(stats3, 'Stats_acute_together_LME.xlsx')





#do stats and median FC separating acute and ERS, make the matrix with all timepoints #do with lmer and emmeans
#ratio 6 is with the repeat samples
#use acute as base
Ratio6$Group2 <- Ratio6$Timepoint
Ratio6$Group2[as.numeric(Ratio6$days_since_admission) > 11] <- 'ERS'
Ratio6$Group2[Ratio6$Timepoint == 'Convalescent'] <- 'Convalescent'
stats2 <- matrix(nrow = 34, ncol = 6)
rownames(stats2)  <- colnames(Ratio6)[18:51]
for(i in 18:(ncol(Ratio6)-1)){
  mixed.lmer <- lmerTest::lmer(Ratio6[,i] ~ Group2 + (1|Donor), 
                               data = Ratio6)
  
  sig <- emmeans(mixed.lmer, pairwise ~ Group2, adjust = "tukey")
  stats2[i-17,] <- summary(sig)[[2]][,6]
  colnames(stats2) <- summary(sig)[[2]][,1]

}


stats2 <- data.frame(stats2)
write.xlsx(stats2, 'Stats_acute_ERS_LME_non_corrected.xlsx')

stats2.subsets <- stats2[!rownames(stats2) %in% c('Granulocytes', 'Monocytes', 'CD4T', 
                                                  'CD8T', 'mDC', 'pDC', 
                                                  'B_cells',  'NK'),]
for(i in 1:ncol(stats2.subsets)){
  stats2.subsets[,i] <-  p.adjust(stats2.subsets[,i], method = 'BH')
}


stats2.lineages <- stats2[rownames(stats2) %in% c('Granulocytes', 'Monocytes', 'CD4T', 
                                                  'CD8T', 'mDC', 'pDC', 
                                                  'B_cells',  'NK'),]
for(i in 1:ncol(stats2.lineages)){
  stats2.lineages[,i] <-  p.adjust(stats2.lineages[,i], method = 'BH')
}


stats3 <- rbind(stats2.lineages, stats2.subsets)
write.xlsx(stats3, 'Stats_acute_ERS_LME_.xlsx')




################################################################################
#laternative stats
#stats with just the model and no post test, comparing to healthy donors
##############################################################################


#do stats with days since admission, no HD
#ratio 6 is with the repeat samples
#use acute as base
stats2 <- matrix(nrow = 34, ncol = 5)
rownames(stats2)  <- colnames(Ratio6)[18:51]
for(i in 18:(ncol(Ratio6)-1)){
  mixed.lmer <- lmerTest::lmer(Ratio6[,i] ~ days_since_admission + (1|Donor), 
                               data = Ratio6)
  stats2[i-17,] <- summary(mixed.lmer)[[10]][2,]
  colnames(stats2) <- colnames(summary(mixed.lmer)[[10]])
  
}


stats2 <- data.frame(stats2)
stats2.subsets <- stats2[!rownames(stats2) %in% c('Granulocytes', 'Monocytes', 'CD4T', 
                                                  'CD8T', 'mDC', 'pDC', 
                                                  'B_cells',  'NK'),]
for(i in 1:ncol(stats2.subsets)){
  stats2.subsets[,i] <-  p.adjust(stats2.subsets[,i], method = 'BH')
}


stats2.lineages <- stats2[rownames(stats2) %in% c('Granulocytes', 'Monocytes', 'CD4T', 
                                                  'CD8T', 'mDC', 'pDC', 
                                                  'B_cells',  'NK'),]
for(i in 1:ncol(stats2.lineages)){
  stats2.lineages[,i] <-  p.adjust(stats2.lineages[,i], method = 'BH')
}


stats3 <- rbind(stats2.lineages, stats2.subsets)
write.xlsx(stats3, 'Stats_acute_early_LME_linear_Acbase.xlsx')





#do stats and median FC separating acute and ERS, make the matrix with all timepoints #do with lmer and emmeans
#ratio 6 is with the repeat samples
#use HC as base and get out the estimate and std errors of models
stats2 <- matrix(nrow = 34, ncol = 12)
rownames(stats2)  <- colnames(Ratio6)[18:51]
colnames(stats2) <- c('HD Est', 'Acute Est', 'ERS Est', 'Convalescent Est',
                      'HD StdError', 'Acute StdError', 'ERS StdError', 'Convalescent StdError',
                      'HD Pval', 'Acute Pval', 'ERS Pval', 'Convalescent Pval')
Ratio60 <- Ratio6
Ratio60$Group2 <- factor(Ratio60$Group2, levels = c('HD', 'Acute', 'ERS', 'Convalescent'))
for(i in 18:(ncol(Ratio6)-1)){
  mixed.lmer <- lmerTest::lmer(Ratio60[,i] ~ Group2 + (1|Donor), 
                               data = Ratio60)
  
  stats2[i-17,c(1:4)] <- summary(mixed.lmer)[[10]][1:4,1]
  stats2[i-17,c(5:8)] <- summary(mixed.lmer)[[10]][1:4,2]
  stats2[i-17,c(9:12)] <- summary(mixed.lmer)[[10]][1:4,5]
  
}

stats2 <- data.frame(stats2)
write.xlsx(stats2, 'model_4groups_HDbase.xlsx')

stats2.subsets <- stats2[!rownames(stats2) %in% c('Granulocytes', 'Monocytes', 'CD4T', 
                                                  'CD8T', 'mDC', 'pDC', 
                                                  'B_cells',  'NK'),]
for(i in c(2,4,6)){
  stats2.subsets[,i] <-  p.adjust(stats2.subsets[,i], method = 'BH')
}
stats2.lineages <- stats2[rownames(stats2) %in% c('Granulocytes', 'Monocytes', 'CD4T', 
                                                'CD8T', 'mDC', 'pDC', 
                                                  'B_cells',  'NK'),]
for(i in c(2,4,6)){
  stats2.lineages[,i] <-  p.adjust(stats2.lineages[,i], method = 'BH')
}
stats3 <- rbind(stats2.lineages, stats2.subsets)
stats4 <- merge(stats2, stats3[,c(2,4,6)], by = 0)
colnames(stats4)[8:10] <- c('Acute_padj', 'ERS_padj', 'Conv_padj')
stats4 <- stats4[,c(1:3, 8, 4:5, 9, 6:7, 10)]
colnames(stats4) <- gsub('.x', '', colnames(stats4))
write.xlsx(stats4, 'model_4groups_HDbase.xlsx')





#do stats with lmer and emmeans, only 3 groups
#ratio 6 is with the repeat samples
#use HC as base
stats2 <- matrix(nrow = 34, ncol = 9)
rownames(stats2)  <- colnames(Ratio6)[18:51]
colnames(stats2) <- c('HD Est', 'Acute Est',  'Convalescent Est', 
                      'HD StdDev', 'Acute StdDev',  'Convalescent StdDev',
                      'HD pval', 'Acute pval',  'Convalescent pval')
Ratio60 <- Ratio6
Ratio60$Timepoint <- factor(Ratio60$Timepoint, levels = c('HD', 'Acute', 'Convalescent'))
for(i in 18:(ncol(Ratio6)-1)){
  mixed.lmer <- lmerTest::lmer(Ratio60[,i] ~ Timepoint + (1|Donor), 
                               data = Ratio60)
  
  stats2[i-17,c(1:3)] <- summary(mixed.lmer)[[10]][1:3,1]
  stats2[i-17,c(4:6)] <- summary(mixed.lmer)[[10]][1:3,2]
  stats2[i-17,c(7:9)] <- summary(mixed.lmer)[[10]][1:3,5]
  
  
}


stats2 <- data.frame(stats2)
write.xlsx(stats2, 'model_3groups_HDbase_nocorrection.xlsx')
stats2.subsets <- stats2[!rownames(stats2) %in% c('Granulocytes', 'Monocytes', 'CD4T', 
                                                  'CD8T', 'mDC', 'pDC', 
                                                  'B_cells',  'NK'),]
for(i in c(2,4)){
  stats2.subsets[,i] <-  p.adjust(stats2.subsets[,i], method = 'BH')
}
stats2.lineages <- stats2[rownames(stats2) %in% c('Granulocytes', 'Monocytes', 'CD4T', 
                                                  'CD8T', 'mDC', 'pDC', 
                                                  'B_cells',  'NK'),]
for(i in c(2,4)){
  stats2.lineages[,i] <-  p.adjust(stats2.lineages[,i], method = 'BH')
}
stats3 <- rbind(stats2.lineages, stats2.subsets)
stats4 <- merge(stats2, stats3[,c(2,4)], by = 0)
colnames(stats4)[6:7] <- c('Acute_padj',  'Conv_padj')
stats4 <- stats4[,c(1:3, 6, 4:5, 7)]
colnames(stats4) <- gsub('.x', '', colnames(stats4))
write.xlsx(stats4, 'model_3groups_HDbase.xlsx')



##############################################################################
#plot per donor
Ratio6$Group2 <- factor(Ratio6$Group2, levels = c( 'Acute', 'ERS', 'Convalescent', 'HD'))
Ratio6.arr <- Ratio6[order(Ratio6$Group2),]
heatmap <- t(scale(Ratio6.arr[,18:45]))
colnames(heatmap) <- Ratio6.arr$CyTOF_ID
heatmap[heatmap > 2] <- 2
heatmap[heatmap < -2] <- -2

ann <- data.frame(Group  = Ratio6.arr$Group2)
rownames(ann) <- Ratio6.arr$CyTOF_ID
dev.off()
pdf('heatmap_subsets.pdf', width = 10, height = 4)
pheatmap(heatmap,
         annotation_col = ann, 
         color             = inferno(10),
         cluster_cols = F,
         treeheight_row = 0
                  )
dev.off()




Ratio66.melt <- melt(Ratio6, id = colnames(Ratio6)[c(1:17,52)])
Ratio66.melt$Group2 <- factor(Ratio66.melt$Group2, 
                              levels = c('Acute', 'ERS', 'Convalescent', 'HD'))
ggplot(Ratio66.melt, aes(x=Group2, y=value, group = Donor, colour = Group2)) + 
  geom_boxplot(outlier.colour = NA, aes(group = Group2)) + 
  geom_point() + 
  geom_line(colour = 'grey') +
  facet_wrap(~variable, scales = 'free_y') + 
  theme_bw() +
  ylab('Log10 Ratio to epithelial cells') +
  theme(legend.position = 'none') +
  theme(strip.background = element_rect(color=NA, fill=NA)) +
  scale_colour_manual(values = c('#e34a33', '#fdbb84', '#fa9fb5',  '#2b8cbe')) 
ggsave('ratio to ep_2acute_groups_colours_incGran.pdf', width = 9, height = 9, useDingbats = F)



Ratio66.melt$days_since_admission[Ratio66.melt$Timepoint == 'HD'] <- 125
Ratio66.melt$patient <- 'yes'
Ratio66.melt$patient[Ratio6$Timepoint == 'HD'] <- 'no'

ggplot(Ratio66.melt, aes(x=days_since_admission, y=value)) + 
  geom_point(aes(colour = Group2)) + 
  geom_line(aes(group = Donor), colour = 'grey') + 
  facet_wrap(~variable, scales = 'free_y') + 
  stat_smooth(colour = 'black') +
  theme_bw() +
  ylab('Log10 Ratio to epithelial cells') +
  theme(strip.background = element_rect(color=NA, fill=NA)) +
  scale_colour_manual(values = c('#e34a33', '#fdbb84', '#fa9fb5',  '#2b8cbe'))
ggsave('Ratio to Ep_admission_incHD_4groups_plines.pdf', width = 10, height = 10)

ggplot(Ratio66.melt, aes(x=days_since_admission, y=value)) + 
  geom_point(aes(colour = Group2)) + 
  facet_wrap(~variable, scales = 'free_y') + 
  stat_smooth(colour = 'black', aes(group = patient)) +
  theme_bw() +
  ylab('Log10 Ratio to epithelial cells') +
  theme(strip.background = element_rect(color=NA, fill=NA)) +
  scale_colour_manual(values = c('#e34a33', '#fdbb84', '#fa9fb5',  '#2b8cbe')) + 
  geom_vline(xintercept = 120, linetype = 'dashed')
ggsave('Ratio to Ep_admission_incHD_4groups_nolines.pdf', width = 10, height = 10)


ggplot(Ratio66.melt, aes(x=days_since_admission, y=value)) + 
  geom_point(aes(colour = Group2)) + 
  facet_wrap(~variable, scales = 'free_y') + 
  stat_smooth(colour = 'black') +
  theme_bw() +
  ylab('Log10 Ratio to epithelial cells') +
  theme(strip.background = element_rect(color=NA, fill=NA)) +
  scale_colour_manual(values = c('#e34a33', '#fdbb84', '#fa9fb5',  '#2b8cbe')) + 
  geom_vline(xintercept = 120, linetype = 'dashed')
ggsave('Ratio to Ep_admission_incHD_4groups_nolines_linetoend.pdf', width = 10, height = 10)

#only 1 timepoint per donor
Ratio76.melt <- Ratio66.melt[Ratio66.melt$Sample_Number != 'Acute2' & 
                               Ratio66.melt$Sample_Number != 'Acute3' &
                               Ratio66.melt$Sample_Number != 'Acute4',]
ggplot(Ratio76.melt, aes(x=Group2, y=value, group = Donor, colour = Group2)) + 
  geom_boxplot(outlier.colour = NA, aes(group = Group2)) + 
  geom_point() + 
  geom_line(colour = 'grey') +
  facet_wrap(~variable, scales = 'free_y') + 
  theme_bw() +
  ylab('Log10 Ratio to epithelial cells') +
  theme(legend.position = 'none') +
  theme(strip.background = element_rect(color=NA, fill=NA))+
  scale_colour_manual(values = c('#e34a33', '#fdbb84', '#fa9fb5',  '#2b8cbe')) 
ggsave('ratio to ep_2acute_groups_colours_1sample_incGran.pdf', width = 9, 
       height = 9, useDingbats = F)


Ratio86 <- Ratio76.melt[Ratio76.melt$variable %in% c('Granulocytes', 'Monocytes', 'CD4T', 
                                        'CD8T', 'mDC', 'pDC', 
                                        'B_cells',  'NK'),]
ggplot(Ratio86, aes(x=Group2, y=value, group = Donor, colour = Group2)) + 
  geom_boxplot(outlier.colour = NA, aes(group = Group2)) + 
  geom_point() + 
  geom_line(colour = 'grey') +
  facet_wrap(~variable, scales = 'free_y', ncol = 4) + 
  theme_bw() +
  ylab('Log10 Ratio to epithelial cells') +
  theme(legend.position = 'none') +
  theme(strip.background = element_rect(color=NA, fill=NA))+
  scale_colour_manual(values = c('#e34a33', '#fdbb84', '#fa9fb5',  '#2b8cbe')) 
ggsave('ratio to ep_2acute_groups_lineages.pdf', width = 5, 
       height = 4, useDingbats = F)










ggplot(Ratio66.melt, aes(x=severity_score, y=value)) + 
  geom_point(aes(colour = Group2)) + 
  facet_wrap(~variable, scales = 'free_y') + 
  stat_smooth(colour = 'black', aes(group = patient)) +
  theme_bw() +
  ylab('Log10 Ratio to epithelial cells') +
  theme(strip.background = element_rect(color=NA, fill=NA)) +
  scale_colour_manual(values = c('#e34a33', '#fdbb84', '#fa9fb5',  '#2b8cbe')) 

Ratio666.melt <- Ratio66.melt[Ratio66.melt$Timepoint == 'Acute',]
ggplot(Ratio666.melt, aes(x=days_since_discharge, y=value)) + 
  geom_point(aes(colour = Group2)) + 
  facet_wrap(~variable, scales = 'free_y') + 
  stat_smooth(colour = 'black', aes(group = patient)) +
  theme_bw() +
  ylab('Log10 Ratio to epithelial cells') +
  theme(strip.background = element_rect(color=NA, fill=NA)) +
  scale_colour_manual(values = c('#e34a33', '#fdbb84', '#fa9fb5',  '#2b8cbe')) 



#MDS
Ratio32 <- Ratio5[,c(1:73, 108:114, 74:107)]

ctrls <- read.xlsx('C:/Users/spjochems/Dropbox/LUMC/Results/CyTOF/BEAT-COVID/Demographics healthy controls BEAT_COVID.xlsx',1)
ctrls$gender <- 'male'
ctrls$gender[ctrls$Geslacht == 'Vrouw'] <- 'female'

Ratio32$Age <- 2020 - Ratio32$pat_yob
Ratio32$gender[Ratio32$gender == 1] <- 'male'
Ratio32$gender[Ratio32$gender == 2] <- 'female'


for(i in 1:nrow(Ratio32)){
  if(Ratio32$Timepoint[i] == 'HD'){
    d <- Ratio32$Donor[i]
    Ratio32$Age[i] <- ctrls$Age[ctrls$BEAT.COVID.nr == d]
    Ratio32$gender[i] <- ctrls$gender[ctrls$BEAT.COVID.nr == d]
  }
}




Ratio32$Group2 <- Ratio32$Timepoint
Ratio32$Group2[as.numeric(Ratio32$days_since_admission) > 11] <- 'ERS'
Ratio32$Group2[Ratio32$Timepoint == 'Convalescent'] <- 'Convalescent'
Ratio32$Group2 <- factor(Ratio32$Group2, levels = c( 'Acute', 'ERS', 'Convalescent', 'HD'))

d <- dist(Ratio32[,81:108]) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
fit # view results
x <- fit$points[,1]
y <- fit$points[,2]
fitted <- data.frame(fit$points)
fit2 <- cbind(fitted, Ratio32)
fit2$Link <- paste(fit2$Donor, fit2$Timepoint, sep = '_')

centroids <- cbind(aggregate(fit2$X1~fit2$Group2, FUN=mean), aggregate(fit2$X2~fit2$Group2, FUN=mean))
centroids <- centroids[,2:4]
colnames(centroids) <- c('Av1', 'Group2', 'Av2')

ggplot(fit2, aes(x=X1, y=X2, Group = Group2, colour = Group2)) +
  geom_point(aes(shape = Sample_Number)) +
  geom_line(aes(Group = Donor), colour = 'grey')+
  geom_point(data = centroids, mapping = aes(x = Av1, y = Av2), shape = 'diamond', size = 3) +
  theme_bw() + 
  xlab('MDS_1') +
  ylab('MDS_2') +
  geom_text(aes(label = Donor)) +
  scale_colour_manual(values = c('#e34a33', '#fdbb84', '#fa9fb5',  '#2b8cbe')) 
#ggsave('MDS_Ratio.pdf', width = 8, height = 5, useDingbats=F)



ggplot(fit2, aes(x=X1, y=X2, Group = Group2, colour = Group2)) +
  geom_point(aes(shape = Sample_Number), size = 2) +
  geom_line(aes(Group = Donor), colour = 'grey', alpha = 0.8)+
  geom_point(data = centroids, mapping = aes(x = Av1, y = Av2, fill = Group2), 
             pch = 23, colour ='black', size = 4) +
  theme_bw() + 
  xlab('MDS_1') +
  ylab('MDS_2')  +
  theme(aspect.ratio = 1) + 
  scale_colour_manual(values = c('#e34a33', '#fdbb84', '#fa9fb5',  '#2b8cbe')) +
  scale_fill_manual(values = c('#e34a33', '#fdbb84', '#fa9fb5',  '#2b8cbe')) 
ggsave('MDS_Ratio_nonames.pdf', width = 11, height = 11, units = 'cm', useDingbats=F)


fit2$gender[fit2$Donor == 'BEAT-116'] <- 'female'
ggplot(fit2, aes(x=X1, y=X2, Group = Group2, colour = gender)) +
  geom_point(size = 2) +
  geom_line(aes(Group = Donor), colour = 'grey', alpha = 0.8)+
  theme_bw() + 
  xlab('MDS_1') +
  ylab('MDS_2')   +
  facet_grid(.~Group2) +
  theme(aspect.ratio = 1) 
ggsave('MDS_Ratio_sex.pdf', width = 20, height = 8, units = 'cm', useDingbats=F)


fit2$Age[fit2$Donor == 'BEAT-116'] <- 70
ggplot(fit2, aes(x=X1, y=X2, Group = Group2, colour = Age)) +
  geom_point(size = 2) +
  geom_line(aes(Group = Donor), colour = 'grey', alpha = 0.8)+
  theme_bw() + 
  xlab('MDS_1') +
  ylab('MDS_2')   +
  facet_grid(.~Group2) +
  theme(aspect.ratio = 1)+
  scale_colour_viridis_c()
ggsave('MDS_Ratio_age.pdf', width = 20, height = 8, units = 'cm', useDingbats=F)


fit2$admission_bmi[fit2$Donor == 'BEAT-116'] <- 30
ggplot(fit2, aes(x=X1, y=X2, Group = Group2, colour = as.numeric(admission_bmi))) +
  geom_point(size = 2) +
  geom_line(aes(Group = Donor), colour = 'grey', alpha = 0.8)+
  theme_bw() + 
  xlab('MDS_1') +
  ylab('MDS_2')   +
  facet_grid(.~Group2) +
  theme(aspect.ratio = 1)+
  scale_colour_viridis_c()
ggsave('MDS_Ratio_BMI.pdf', width = 20, height = 8, units = 'cm', useDingbats=F)


ggplot(fit2, aes(x=X1, y=X2, Group = Group2, colour = as.numeric(days_since_discharge))) +
  geom_point(size = 2) +
  geom_line(aes(Group = Donor), colour = 'grey', alpha = 0.8)+
  theme_bw() + 
  xlab('MDS_1') +
  ylab('MDS_2')   +
  facet_grid(.~Group2) +
  theme(aspect.ratio = 1) +
  scale_colour_viridis_c()
ggsave('MDS_Ratio_discharge.pdf', width = 20, height = 8, units = 'cm', useDingbats=F)


ggplot(fit2, aes(x=X1, y=X2, Group = Group2, colour = as.numeric(severity_score))) +
  geom_point(size = 2) +
  geom_line(aes(Group = Donor), colour = 'grey', alpha = 0.8)+
  theme_bw() + 
  xlab('MDS_1') +
  ylab('MDS_2')   +
  facet_grid(.~Group2) +
  theme(aspect.ratio = 1) +
  scale_colour_viridis_c()
ggsave('MDS_Ratio_severity.pdf', width = 20, height = 8, units = 'cm', useDingbats=F)


meds <- read.xlsx('../../medication nasal patients_checked AR_2.xlsx', 1)
meds[,3][is.na(meds[,3])] <- 'No'
meds[,4][is.na(meds[,4])] <- 'No'
fit3 <- merge(fit2, meds, by = 'Donor', all.x=T)

ggplot(fit3, aes(x=X1, y=X2, Group = Group2, colour = Steroids)) +
  geom_point(size = 2) +
  geom_line(aes(Group = Donor), colour = 'grey', alpha = 0.8)+
  theme_bw() + 
  xlab('MDS_1') +
  ylab('MDS_2')   +
  facet_grid(.~Group2) +
  theme(aspect.ratio = 1) 
ggsave('MDS_Ratio_medication.pdf', width = 20, height = 8, units = 'cm', useDingbats=F)


fit3$asthma[fit3$asthma == 1] <- 'Yes'
fit3$asthma[fit3$asthma == 2] <- 'No'
fit3$asthma[fit3$asthma == 3] <- 'Unknown'


fit3$asthma[fit3$Donor == 'BEAT-116'] <- 'No'
ggplot(fit3, aes(x=X1, y=X2, Group = Group2, colour = asthma)) +
  geom_point(size = 2) +
  geom_line(aes(Group = Donor), colour = 'grey', alpha = 0.8)+
  theme_bw() + 
  xlab('MDS_1') +
  ylab('MDS_2')   +
  facet_grid(.~Group2) +
  theme(aspect.ratio = 1) 
ggsave('MDS_Ratio_asthma.pdf', width = 20, height = 8, units = 'cm', useDingbats=F)


fit3$diabetes[fit3$diabetes == 1] <- 'Yes'
fit3$diabetes[fit3$diabetes == 2] <- 'No'
fit3$diabetes[fit3$Donor == 'BEAT-116'] <- 'No'
ggplot(fit3, aes(x=X1, y=X2, Group = Group2, colour = diabetes)) +
  geom_point(size = 2) +
  geom_line(aes(Group = Donor), colour = 'grey', alpha = 0.8)+
  theme_bw() + 
  xlab('MDS_1') +
  ylab('MDS_2')   +
  facet_grid(.~Group2) +
  theme(aspect.ratio = 1) 
ggsave('MDS_Ratio_diabetes.pdf', width = 20, height = 8, units = 'cm', useDingbats=F)



fup <- read.delim('../6W_clinical.txt')
fit4 <- merge( fit3, fup, by.x = 'Donor', by.y = 'Donor_ID', all.x=T)
ggplot(fit4, aes(x=X1, y=X2, Group = Group2, colour = diff_dr)) +
  geom_point(size = 2) +
  geom_line(aes(Group = Donor), colour = 'grey', alpha = 0.8)+
  theme_bw() + 
  xlab('MDS_1') +
  ylab('MDS_2')   +
  facet_grid(.~Group2) +
  theme(aspect.ratio = 1) 



#save.image('Analysis.RData')
#VL data
write.table(Ratio32, 'Ratio32.txt', row.names = F, sep = '\t', quote =F)
Ratio42 <- read.delim('Ratio32_wVL.txt') #manually added in excel


VL <- read.delim('../VLs.txt')
NP <- VL[VL$Sample.type == 'nasopharyngeal swab',]
sputum <- VL[VL$Sample.type == 'sputum',]


minNP <- aggregate(NP[,4]~NP[,1], FUN=min)
minsputum <- aggregate(sputum[,4]~sputum[,1], FUN=min)
colnames(minNP) <- c('Donor', 'MinCtNPS') 
colnames(minsputum) <- c('Donor', 'MinCtSputum') 
minall <- merge(minNP, minsputum, by = 'Donor', all = T)
minall$MinEither <- minall$MinCtNPS
minall$MinEither[is.na(minall$MinCtNPS)] <- minall$MinCtSputum[is.na(minall$MinCtNPS)]
minall$MinEither[which(minall$MinCtSputum < minall$MinCtNPS)] <- minall$MinCtSputum[which(minall$MinCtSputum < minall$MinCtNPS)] 

Ratio42 <- merge(Ratio42, minall, by = 'Donor', all.x = T)
#correlations of data
Acutes <- Ratio42[Ratio42$Timepoint == 'Acute',]
cor_data <- Acutes[,c(60,44,45,46,49,57, 81:117,66)]

for(i in 2:ncol(cor_data)){
  cor_data[,i] <- gsub(' days', '', cor_data[,i])
  cor_data[,i] <- as.numeric(cor_data[,i])
}
cor_data$VL <- 40- cor_data$VL

corr <- rcorr(as.matrix(cor_data[,2:ncol(cor_data)]), 
              type = 'spearman')
pdf('corrplot_inc_VL.pdf', width = 10, height = 10)
corrplot::corrplot(corr$r, 
                   p.mat = corr$P, 
                   sig.level = 0.05, 
                   insig = 'label_sig', 
                   pch.col = 'black',
                   pch.cex = 0.8, 
                   type= 'upper',
                   order = 'hclust')
dev.off()


#correlation without severity score
corr <- rcorr(as.matrix(cor_data[,c(2:5,7:ncol(cor_data))]), 
              type = 'spearman')
pdf('corrplot_inc_VL_no_severity score.pdf', width = 12, height = 12)
corrplot::corrplot(corr$r, 
                   pch.cex = 0.8, 
                   type= 'upper',
                   method = "number",
                   order = 'hclust', 
                   number.cex = 0.5)
dev.off()

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
pdf('corrplot_inc_VL_no_severity score.pdf', width = 10, height = 10)
corrplot::corrplot(corr$r, 
         method = "color", 
         col = rev(col(200)),
         type = "upper", 
         order = "hclust", 
         tl.col = "black", 
         tl.srt = 90,
         p.mat = corr$P, 
         sig.level = c(0.001,0.01,0.05), 
         insig = 'label_sig', 
         pch.col = 'white',
         pch.cex = 0.5, )
dev.off()




pdf('corrplot_inc_VL_no_severity score_wnumbers.pdf', width = 12, height = 12)
corrplot::corrplot(corr$r, 
                   method = "color", 
                   col = rev(col(200)),
                   type = "upper", 
                   order = "hclust", 
                   number.cex = .5,
                   addCoef.col = "black", # Add coefficient of correlation
                   tl.col = "black", 
                   tl.srt = 90)
dev.off()


cs <- data.frame(t(counts))
cs$Mono <- cs$Myeloid_CD206_CD163_Mono + cs$Myeloid_CD163_Mono +
  cs$Myeloid_Classical_Mono 
cs$CD163_PercMono <- cs$Myeloid_CD163_Mono / cs$Mono  * 100
cs$CD163CD206_PercMono <- cs$Myeloid_CD206_CD163_Mono / cs$Mono  * 100

cs$DC <- cs$Myeloid_CD163_DC + cs$Myeloid_CD206_DC
cs$CD206_PercDC <- cs$Myeloid_CD206_DC / cs$DC * 100

cor_data2 <- merge(cor_data, cs[,c(30:31,33)], by.x = 'CyTOF_ID', by.y = 0)



mfi <- read.delim('../tsne and heatmaps/Clustered_data_intensity.txt')
mfi0 <- mfi[mfi$Pop == 'CD8T TRM',]
group <- split(mfi0, f = mfi0$Vol)
quant <- lapply(group, function(x){
  a <- x
  a <- a[,5:ncol(a)]
  apply(a, 2, quantile, .5)
})

quantiles <- do.call(rbind, quant)
colnames(quantiles) <- paste0(colnames(quantiles), '_MFI_CD8TRM')

cor_data3 <- merge(cor_data2, quantiles[,c(25,27,31)], 
                   by.x = 'CyTOF_ID', by.y = 0, all.x = T)


mfi0 <- mfi[mfi$Pop == 'CD4T_EM',]
group <- split(mfi0, f = mfi0$Vol)
quant <- lapply(group, function(x){
  a <- x
  a <- a[,5:ncol(a)]
  apply(a, 2, quantile, .5)
})

quantiles <- do.call(rbind, quant)
colnames(quantiles) <- paste0(colnames(quantiles), '_MFI_CD4TEM')

cor_data4 <- merge(cor_data3, quantiles[,c(22,34)], 
                   by.x = 'CyTOF_ID', by.y = 0, all.x = T)





corr <- rcorr(as.matrix(cor_data4[,2:ncol(cor_data4)]), 
              type = 'spearman')
pdf('corrplot_incactivation.pdf', width = 10, height = 10)
corrplot::corrplot(corr$r, 
                   p.mat = corr$P, 
                   sig.level = 0.05, 
                   insig = 'label_sig', 
                   pch.col = 'black',
                   pch.cex = 0.8, 
                   type= 'upper',
                   order = 'hclust')
dev.off()


SPO <- Ratio32
SPO$spo <- NA
SPO$spo[which(Ratio32$spo2 > 93)] <- '94+'
SPO$spo[which(Ratio32$spo2 < 94)] <- '93-'
ggplot(SPO, aes(x=spo2, y = Granulocytes, colour = Group2)) + geom_point()
ggplot(SPO, aes(x=spo2, y = Gran_CD16Hi, colour = Group2)) + geom_point() +
  stat_smooth(method = 'lm', aes(group=1))

column_dend = hclust(dist(as.matrix(t(cor_data))))
#corr$r[corr$P > 0.05] <- 0
dev.off()
pdf('correlation_heatmap.pdf', width = 8, height = 8)
pheatmap(corr$r, 
         cluster_cols = column_dend,
         cluster_rows = column_dend)
dev.off()



splitted <- split(Ratio32, f = Ratio32$Timepoint)
breaksList = seq(-1, 1, by = 0.01)


ac <- splitted[[1]]
ac <- ac[ac$Timepoint != 'Acute2' & 
           ac$Timepoint != 'Acute3' &
           ac$Timepoint != 'Acute4',]
cor_data <- ac[,15:42]
corr <- rcorr(as.matrix(cor_data))



pheatmap(corr$r, 
         cluster_cols = column_dend,
         cluster_rows = column_dend,
         border_color = NA,
         treeheight_col = 0,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaksList)),
         breaks = breaksList)


cor_data <- splitted[[2]][,15:42]
corr <- rcorr(as.matrix(cor_data))
pheatmap(corr$r, 
         cluster_cols = column_dend,
         cluster_rows = column_dend,
         border_color = NA,
         treeheight_col = 0,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaksList)),
         breaks = breaksList)



cor_data <- splitted[[3]][,15:42]
corr <- rcorr(as.matrix(cor_data))
pheatmap(corr$r, 
         cluster_cols = column_dend,
         cluster_rows = column_dend,
         border_color = NA,
         treeheight_col = 0,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaksList)),
         breaks = breaksList)



acconv <- merge(ac[,c(4,15:42)], splitted[[2]][,c(4,15:42)], by = 'Donor')
corr <- rcorr(as.matrix(acconv[,2:ncol(acconv)]))
corr2 <- corr$r[1:28,29:56]
corp2 <- corr$P[1:28,29:56]
corr2[corp2 > 0.05] <- 0
pheatmap(corr2, 
        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaksList)),
         breaks = breaksList)


scaled <- scale(Ratio32[,15:42])
rownames(scaled) <- Ratio32$CyTOF_ID
ann <- data.frame(Group = Ratio32$Group2)
rownames(ann) <-  Ratio2$CyTOF_ID
pheatmap(scaled, annotation_row = ann)

