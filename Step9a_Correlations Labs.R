#===============================================================================================
# Expression data in rat, mouse and human muscle, L6, C2C12 and human myotubes
#===============================================================================================
# This file assumes that all data has been processed and normalized as in "Step1_Data.Rds"
# and is available as one single file: "GENENAME_norm.Rds" 
library(here)
library(limma)
library(impute)
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(ggpubr)
library(gplots)
library(car)


#=============================================================================================================================
# Diffferences in transcriptomic accross laboratories
#=============================================================================================================================
# Here we want to know how different are samples comming from different laboratories.
muscle <- readRDS(here("Data_processed", "GENENAME_norm.Rds"))
muscle <- muscle[!grepl('HEK', colnames(muscle))] #remove HEK cells
muscle <- muscle[!grepl('HeLa', colnames(muscle))] #remove HeLa cells


C2C12 <- muscle[grepl('MouseC2C12', colnames(muscle))]
colnames(C2C12)
C2C12 <- data.frame(GSE44563=rowMeans(C2C12[,1:3], na.rm=T),
                    GSE11415=rowMeans(C2C12[,4:9], na.rm=T),
                    GSE12993=rowMeans(C2C12[,10:11], na.rm=T),
                    GSE79925=rowMeans(C2C12[,12:15], na.rm=T),
                    GSE1984=rowMeans(C2C12[,16:19], na.rm=T),
                    GSE4330=rowMeans(C2C12[,20:22], na.rm=T))
C2C12 <- na.omit(C2C12)
C2C12.cor <- cor(C2C12, method="spearman")
C2C12.cor[C2C12.cor==1] <- NA
#write.table(C2C12.cor, file="CorLabs_C2C12.txt", sep="\t")


HSMC <- muscle[grepl('HumanCell', colnames(muscle))]
colnames(HSMC)
HSMC <- data.frame(GSE70822=rowMeans(HSMC[,1:6], na.rm=T),
                   GSE55034=rowMeans(HSMC[,7:8], na.rm=T),
                   GSE67326=rowMeans(HSMC[,9:16], na.rm=T),
                   GSE40789=rowMeans(HSMC[,17:20], na.rm=T),
                   GSE44051=rowMeans(HSMC[,21:32], na.rm=T))
HSMC <- na.omit(HSMC)
HSMC.cor <- cor(HSMC)
HSMC.cor[HSMC.cor==1] <- NA



png(filename=here("Figures", "CorrelationsLabs_C2C12.png"), #print graph
    units="cm", width=33, height=33, 
    pointsize=12, res=300)
plot(C2C12)
dev.off()

png(filename=here("Figures", "CorrelationsLabs_HSMC.png"), #print graph
    units="cm", width=33, height=33, 
    pointsize=12, res=300)
plot(HSMC)
dev.off()

