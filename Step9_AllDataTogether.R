setwd("C:/ownCloud/Projects/MuscleModels")
library(readxl)


#Transcriptome
transcriptome <- readRDS("Data/Transcriptomics/GENENAME_norm.Rds")
mRNA_L6 <- transcriptome[grep('RatL6', colnames(transcriptome))]
colnames(mRNA_L6) <- paste('L6_', seq(1,ncol(mRNA_L6),1), sep='')

mRNA_C2 <- transcriptome[grep('MouseC2C12', colnames(transcriptome))] 
mRNA_HSMC <- transcriptome[grep('HumanCell', colnames(transcriptome))] 


#FA oxidation
FAOx <- read_excel("Data/Metabolism/FA oxidation.xlsx", sheet = "Sheet1")

FAOx_L6 <- data.frame(rbind(FAOx_L6_basal=FAOx$L6_basal, FAOx_L6_FC=FAOx$L6_FC))
colnames(FAOx_L6) <- paste('L6_', seq(1,ncol(FAOx_L6),1), sep='')

FAOx_C2 <- data.frame(rbind(FAOx_C2_basal=FAOx$C2C12_basal, FAOx_C2_FC=FAOx$C2C12_FC))
colnames(FAOx_C2) <- paste('C2C12_', seq(1,ncol(FAOx_C2),1), sep='')

FAOx_HSMC <- data.frame(rbind(FAOx_HSMC_basal=FAOx$HSMC_basal, FAOx_HSMC_FC=FAOx$HSMC_FC))
colnames(FAOx_HSMC) <- paste('HSMC_', seq(1,ncol(FAOx_HSMC),1), sep='')


#Glucose oxidation
GOx <- read_excel("Data/Metabolism/Glc oxidation.xlsx", sheet = "Sheet1")

GOx_L6 <- data.frame(rbind(GOx_L6_basal=GOx$L6_basal, GOx_L6_FC=GOx$L6_FC))
colnames(GOx_L6) <- paste('L6_', seq(1,ncol(GOx_L6),1), sep='')

GOx_C2 <- data.frame(rbind(GOx_C2_basal=GOx$C2C12_basal, GOx_C2_FC=GOx$C2C12_FC))
colnames(GOx_C2) <- paste('C2C12_', seq(1,ncol(GOx_C2),1), sep='')

GOx_HSMC <- data.frame(rbind(GOx_HSMC_basal=GOx$HSMC_basal, GOx_HSMC_FC=GOx$HSMC_FC))
colnames(GOx_HSMC) <- paste('HSMC_', seq(1,ncol(GOx_HSMC),1), sep='')


#Glucose uptake
GUp <- read_excel("Data/Metabolism/Glc uptake.xlsx", sheet = "Sheet1")

GUp_L6 <- data.frame(rbind(GUp_L6_basal=GUp$L6_basal, GUp_L6_FC=GUp$L6_FC))
colnames(GUp_L6) <- paste('L6_', seq(1,ncol(GUp_L6),1), sep='')

GUp_C2 <- data.frame(rbind(GUp_C2_basal=GUp$C2C12_basal, GUp_C2_FC=GUp$C2C12_FC))
colnames(GUp_C2) <- paste('C2C12_', seq(1,ncol(GUp_C2),1), sep='')

GUp_HSMC <- data.frame(rbind(GUp_HSMC_basal=GUp$HSMC_basal, GUp_HSMC_FC=GUp$HSMC_FC))
colnames(GUp_HSMC) <- paste('HSMC_', seq(1,ncol(GUp_HSMC),1), sep='')


#Glycogen synthesis
Glyc <- read_excel("Data/Metabolism/Glyc synthesis.xlsx", sheet = "Sheet1")

Glyc_L6 <- data.frame(rbind(Glyc_L6_basal=Glyc$L6_basal, Glyc_L6_FC=Glyc$L6_FC))
colnames(Glyc_L6) <- paste('L6_', seq(1,ncol(Glyc_L6),1), sep='')

Glyc_C2 <- data.frame(rbind(Glyc_C2_basal=Glyc$C2C12_basal, Glyc_C2_FC=Glyc$C2C12_FC))
colnames(Glyc_C2) <- paste('C2C12_', seq(1,ncol(Glyc_C2),1), sep='')

Glyc_HSMC <- data.frame(rbind(Glyc_HSMC_basal=Glyc$HSMC_basal, Glyc_HSMC_FC=Glyc$HSMC_FC))
colnames(Glyc_HSMC) <- paste('HSMC_', seq(1,ncol(Glyc_HSMC),1), sep='')


#Protein content
Prot <- read_excel("Data/Metabolism/Protein content.xlsx", sheet = "Sheet1")

Prot_L6 <- data.frame(t(c(Prot$L6_GS, Prot$L6_FAO)))
colnames(Prot_L6) <- paste('L6_', seq(1,ncol(Prot_L6),1), sep='')

Prot_C2 <- data.frame(t(c(Prot$C2C12_GS, Prot$C2C12_FAO)))
colnames(Prot_C2) <- paste('C2_', seq(1,ncol(Prot_C2),1), sep='')

Prot_HSMC <- data.frame(t(c(Prot$HSMC_GS, Prot$HSMC_FAO)))
colnames(Prot_HSMC) <- paste('HSMC_', seq(1,ncol(Prot_HSMC),1), sep='')


#bind all
library(dplyr)
L6 <- dplyr::bind_rows(mRNA_L6,
                       FAOx_L6,
                       GOx_L6,
                       GUp_L6,
                       Glyc_L6,
                       Prot_L6)

C2 <- dplyr::bind_rows(mRNA_C2,
                       FAOx_C2,
                       GOx_C2,
                       GUp_C2,
                       Glyc_C2,
                       Prot_C2)

HSMC <- dplyr::bind_rows(mRNA_HSMC,
                         FAOx_HSMC,
                         GOx_HSMC,
                         GUp_HSMC,
                         Glyc_HSMC,
                         Prot_HSMC)


AllData <- cbind(L6, C2, HSMC)
rownames(AllData) <- c(rownames(mRNA_L6), 
                       "FAOx_basal", "FAOx_FC",
                       "GOx_basal", "GOx_FC",
                       "GUp_basal", "GUp_FC",
                       "Glyc_basal", "Glyc_FC",
                       "Prot_L6")
AllData <- AllData[,colSums(is.na(AllData))<nrow(AllData)]


