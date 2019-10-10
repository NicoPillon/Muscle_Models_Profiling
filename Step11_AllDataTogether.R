library(readxl)
library(here)

#Transcriptome
transcriptome <- readRDS(here("Data_Processed", "GENENAME_norm.Rds"))
mRNA_L6 <- transcriptome[grep('RatL6', colnames(transcriptome))]
colnames(mRNA_L6) <- paste('L6_', seq(1,ncol(mRNA_L6),1), sep='')

mRNA_C2 <- transcriptome[grep('MouseC2C12', colnames(transcriptome))]
colnames(mRNA_C2) <- paste('C2C12_', seq(1,ncol(mRNA_C2),1), sep='')

mRNA_HSMC <- transcriptome[grep('HumanCell', colnames(transcriptome))] 
colnames(mRNA_HSMC) <- paste('HSMC_', seq(1,ncol(mRNA_HSMC),1), sep='')

#FA oxidation
FAOx <- read_excel(here("Data_Raw", "Metabolism", "FA oxidation.xlsx"), sheet = "Sheet1")

FAOx_L6 <- data.frame(rbind(FAOx_L6_basal=FAOx$L6_basal, FAOx_L6_FC=FAOx$L6_FC))
colnames(FAOx_L6) <- paste('L6_', seq(1,ncol(FAOx_L6),1), sep='')

FAOx_C2 <- data.frame(rbind(FAOx_C2_basal=FAOx$C2C12_basal, FAOx_C2_FC=FAOx$C2C12_FC))
colnames(FAOx_C2) <- paste('C2C12_', seq(1,ncol(FAOx_C2),1), sep='')

FAOx_HSMC <- data.frame(rbind(FAOx_HSMC_basal=FAOx$HSMC_basal, FAOx_HSMC_FC=FAOx$HSMC_FC))
colnames(FAOx_HSMC) <- paste('HSMC_', seq(1,ncol(FAOx_HSMC),1), sep='')


#Glucose oxidation
GOx <- read_excel(here("Data_Raw", "Metabolism", "Glc oxidation.xlsx"), sheet = "Sheet1")

GOx_L6 <- data.frame(rbind(GOx_L6_basal=GOx$L6_basal, GOx_L6_FC=GOx$L6_FC))
colnames(GOx_L6) <- paste('L6_', seq(1,ncol(GOx_L6),1), sep='')

GOx_C2 <- data.frame(rbind(GOx_C2_basal=GOx$C2C12_basal, GOx_C2_FC=GOx$C2C12_FC))
colnames(GOx_C2) <- paste('C2C12_', seq(1,ncol(GOx_C2),1), sep='')

GOx_HSMC <- data.frame(rbind(GOx_HSMC_basal=GOx$HSMC_basal, GOx_HSMC_FC=GOx$HSMC_FC))
colnames(GOx_HSMC) <- paste('HSMC_', seq(1,ncol(GOx_HSMC),1), sep='')


#Glucose uptake
GUp <- read_excel(here("Data_Raw", "Metabolism", "Glc uptake.xlsx"), sheet = "Sheet1")

GUp_L6 <- data.frame(rbind(GUp_L6_basal=GUp$L6_basal, GUp_L6_FC=GUp$L6_FC))
colnames(GUp_L6) <- paste('L6_', seq(1,ncol(GUp_L6),1), sep='')

GUp_C2 <- data.frame(rbind(GUp_C2_basal=GUp$C2C12_basal, GUp_C2_FC=GUp$C2C12_FC))
colnames(GUp_C2) <- paste('C2C12_', seq(1,ncol(GUp_C2),1), sep='')

GUp_HSMC <- data.frame(rbind(GUp_HSMC_basal=GUp$HSMC_basal, GUp_HSMC_FC=GUp$HSMC_FC))
colnames(GUp_HSMC) <- paste('HSMC_', seq(1,ncol(GUp_HSMC),1), sep='')


#Glycogen synthesis
Glyc <- read_excel(here("Data_Raw", "Metabolism", "Glyc synthesis.xlsx"), sheet = "Sheet1")

Glyc_L6 <- data.frame(rbind(Glyc_L6_basal=Glyc$L6_basal, Glyc_L6_FC=Glyc$L6_FC))
colnames(Glyc_L6) <- paste('L6_', seq(1,ncol(Glyc_L6),1), sep='')

Glyc_C2 <- data.frame(rbind(Glyc_C2_basal=Glyc$C2C12_basal, Glyc_C2_FC=Glyc$C2C12_FC))
colnames(Glyc_C2) <- paste('C2C12_', seq(1,ncol(Glyc_C2),1), sep='')

Glyc_HSMC <- data.frame(rbind(Glyc_HSMC_basal=Glyc$HSMC_basal, Glyc_HSMC_FC=Glyc$HSMC_FC))
colnames(Glyc_HSMC) <- paste('HSMC_', seq(1,ncol(Glyc_HSMC),1), sep='')


#Protein content
Prot <- read_excel(here("Data_Raw", "Protein", "Protein_content.xlsx"), sheet = "Sheet1")

Prot_L6 <- data.frame(t(Prot$L6_GS))
colnames(Prot_L6) <- paste('L6_', seq(1,ncol(Prot_L6),1), sep='')

Prot_C2 <- data.frame(t(Prot$C2C12_GS))
colnames(Prot_C2) <- paste('C2_', seq(1,ncol(Prot_C2),1), sep='')

Prot_HSMC <- data.frame(t(Prot$HSMC_GS))
colnames(Prot_HSMC) <- paste('HSMC_', seq(1,ncol(Prot_HSMC),1), sep='')


#bind all
library(dplyr)
L6 <- dplyr::bind_rows(mRNA_L6,
                       FAOx_L6,
                       GOx_L6,
                       GUp_L6,
                       Glyc_L6)

C2 <- dplyr::bind_rows(mRNA_C2,
                       FAOx_C2,
                       GOx_C2,
                       GUp_C2,
                       Glyc_C2)

HSMC <- dplyr::bind_rows(mRNA_HSMC,
                         FAOx_HSMC,
                         GOx_HSMC,
                         GUp_HSMC,
                         Glyc_HSMC)


AllData <- cbind(L6, C2, HSMC)
rownames(AllData) <- c(rownames(mRNA_L6), 
                       "FAOx_basal", "FAOx_FC",
                       "GOx_basal", "GOx_FC",
                       "GUp_basal", "GUp_FC",
                       "Glyc_basal", "Glyc_FC")

#Make table with mRNA data
GeneData <- AllData[1:(nrow(AllData)-8),]
# delete rows that contain too many NAs
GeneData <- GeneData[rowSums(is.na(GeneData)) / ncol(GeneData) < 25/100, ]

#Make table with metabolic assays data
AssayData <- AllData[c("FAOx_basal", "FAOx_FC",
                       "GOx_basal", "GOx_FC",
                       "GUp_basal", "GUp_FC",
                       "Glyc_basal", "Glyc_FC"),]


#===================================================================================
Spearman.estimate <- function(x) cor.test(x, y, method="spearman", exact=F)$estimate
Spearman.p.value  <- function(x) cor.test(x, y, method="spearman", exact=F)$p.value

estimates <- data.frame(rownames(GeneData))
pvalues <- data.frame(rownames(GeneData))
for (i in 1:nrow(GeneData)){
  y <- as.numeric(AssayData[i,])
  estimates[,i] <-  apply(GeneData, 1, Spearman.estimate)
  pvalues[,i] <-  apply(GeneData, 1, Spearman.p.value)
  print(i)
}
rownames(estimates) <- rownames(GeneData)
rownames(pvalues) <- rownames(GeneData)
colnames(estimates) <- rownames(AssayData)
colnames(pvalues) <- rownames(AssayData)



#===================================================================================
#Keep genes with at least one significant parameter
sign <- pvalues[rowSums(abs(estimates) > 0.9) >= 1, ]
estimates <- estimates[rownames(sign),]

#===================================================================================
library(gplots)
cormatrix <- t(estimates)
#png(filename="Figure_CorrelationMatrix.png", #print graph
#    units="cm", width=40, height=10, 
#    pointsize=12, res=1200)
heatmap.2(cormatrix, scale="none",
          margins =c(round(max(nchar(colnames(cormatrix)))/2),  # adjust margins to col names
                     round(max(nchar(rownames(cormatrix)))/2)), # adjust margins to row names
          key=F, keysize=1,                                  # remove color key legend
          density.info="none",                  # turns off density plot inside color legend
          trace="none",
          col=redgreen,
          cexRow=0.8, cexCol=0.8,
          main="")
dev.off()

sign['SLC7A1',]
rownames(sign)
