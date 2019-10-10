#===============================================================================================
# Expression data in rat, mouse and human mRNA, L6, C2C12 and human myotubes
#===============================================================================================
# This file assumes that all data has been processed and normalized as in "Step1_Data.Rds"
# and is available as one single file: "GENENAME_norm.Rds" 
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(gplots)
library(stringr)
library(grid)
library(gridExtra)
library(here)
cbPalette <- c("#56B4E9", "#D3C839", "#CC79A7", "#0072B2", "#E69F00", "#D55E00") #color palette for colorblind people

#Transcriptome
mRNA <- readRDS(here("Data_processed", "GENENAME_norm.Rds"))
mRNA <- mRNA[!grepl('HEK', colnames(mRNA))] #remove HEK cells
mRNA <- mRNA[!grepl('HeLa', colnames(mRNA))] #remove HeLa cells

#Proteome Rat L6
RatL6_PXD001543 <- read_excel("Data_Raw/Proteomics/RatL6_PXD001543.xlsx")
RatL6_PXD001543 <- RatL6_PXD001543[,c(2,45:53)]
SYMBOL <- gsub("_RAT.*", "", RatL6_PXD001543$`Majority protein IDs`)
RatL6_PXD001543 <- data.frame(SYMBOL=gsub(".*\\|", "", SYMBOL),
                              RatL6=log10(RatL6_PXD001543[,2:10]))
RatL6_PXD001543[RatL6_PXD001543=='-Inf'] <- NA
RatL6_PXD001543 <- RatL6_PXD001543[!is.na(RatL6_PXD001543$SYMBOL),]


#Proteome Rat Tissue
RatTissue_PXD007182 <- read_excel("Data_Raw/Proteomics/RatTissue_PXD007182.xlsx")
RatTissue_PXD007182 <- RatTissue_PXD007182[,c(7,52:57)]
RatTissue_PXD007182 <- data.frame(SYMBOL=toupper(gsub(";.*", "", RatTissue_PXD007182$`Gene names`)),
                                  RatTissue=log10(RatTissue_PXD007182[,2:7]))
RatTissue_PXD007182[RatTissue_PXD007182=='-Inf'] <- NA
RatTissue_PXD007182 <- RatTissue_PXD007182[!is.na(RatTissue_PXD007182$SYMBOL),]


#Proteome Mouse C2C12
MouseC2C12_PXD000288 <- read_excel("Data_Raw/Proteomics/MouseC2C12_PXD000288.xlsx")
MouseC2C12_PXD000288 <- data.frame(SYMBOL=toupper(gsub(";.*", "", MouseC2C12_PXD000288$`Gene names`)),
                                   MouseC2C12=log10(MouseC2C12_PXD000288[,11:16]))
MouseC2C12_PXD000288 <- MouseC2C12_PXD000288[!is.na(MouseC2C12_PXD000288$SYMBOL),]


#Proteome Mouse Tissue
MouseTissue_PXD000288 <- read_excel("Data_Raw/Proteomics/MouseTissue_PXD000288.xlsx")
MouseTissue_PXD000288 <- data.frame(SYMBOL=toupper(gsub(";.*", "", MouseTissue_PXD000288$`Gene names`)),
                                    MouseTissue=log10(MouseTissue_PXD000288[,11:16]))
MouseTissue_PXD000288 <- MouseTissue_PXD000288[!is.na(MouseTissue_PXD000288$SYMBOL),]


# merge all
library(dplyr)
list.data <- list(RatL6_PXD001543, RatTissue_PXD007182,
                  MouseC2C12_PXD000288, MouseTissue_PXD000288)
All  <- data.frame(list.data[[1]]$SYMBOL)
colnames(All) <- 'SYMBOL'
for ( df in list.data ) {
    asd <- df
    All <- full_join(All, asd)
}
All <- aggregate(All[,2:ncol(All)], by=list(All$SYMBOL), FUN=mean, na.rm=T)
All[is.na(All)] <- NA
All <- data.frame(All, row.names = 1)
rm(list.data, asd, df, i, list.filenames)

#normalize
Norm <- normalizeBetweenArrays(All, method='quantile')
Norm <- data.frame(Norm - median(Norm, na.rm=T))

#plot
res <- Norm
Sample <- c(rep("RatL6",  length(grep('RatL6', colnames(res)))),
            rep("RatTissue", length(grep('RatTissue',  colnames(res)))),
            rep("MouseC2C12", length(grep('MouseC2C12',  colnames(res)))),
            rep("MouseTissue",    length(grep('MouseTissue',       colnames(res)))))
PlotFunction <- function(genename) {
    data   <- data.frame()                                #create an empty dataframe to collect data
    for( i in 1:length(genename)) { 
        y     <- as.numeric(res[genename[i],])              #collect data for gene name i
        datay <- cbind.data.frame(Sample, y, rep(genename[i]))   #create table with x="sample type", y="data", "gene name"
        colnames(datay) <- c("x","y","Gene")                #rename column names to make it possible to rbind later
        data  <- rbind.data.frame(data, datay)              #bind the new gene data at the bottom of the previous one
    }
    data$x <- factor(data$x, levels=c("RatL6", "RatTissue", "MouseC2C12", "MouseTissue")) #for a box plot, x should be a factor
    ggplot(data, aes(x=x, y=y, fill=x)) +  
        geom_boxplot() +
        labs(x="",
             y=paste(genename)) + 
        theme_bw() + 
        theme(plot.title  = element_text(face="bold", color="black", size=7, angle=0),
              axis.text.x = element_text(color="black", size=6, angle=45, hjust=1),
              axis.text.y = element_text(color="black", size=6, angle=0),
              axis.title  = element_text(face="bold", color="black", size=7, angle=0),
              legend.text = element_text(face="bold", color="black", size=6, angle=0),
              legend.position="none", legend.title = element_blank()) + 
        scale_fill_manual(values=cbPalette) +
        scale_color_manual(values=cbPalette)
}

PlotFunction('MYH1')















#=============================================================================================================================
# Differences in transcriptomic and proteomics in C2C12
#=============================================================================================================================
C2C12_Prot <- MouseC2C12_PXD000288
C2C12_mRNA <- mRNA[grepl('MouseC2C12', colnames(mRNA))]
C2C12_mRNA <- data.frame(rownames(C2C12_mRNA), rowMeans(C2C12_mRNA))

C2C12 <- merge(C2C12_Prot, C2C12_mRNA, by=1)
colnames(C2C12) <- c('SYMBOL', 'Proteome', 'Transcriptome')
C2C12 <- na.omit(data.frame(aggregate(C2C12[,2:3], b=list(C2C12$SYMBOL), FUN=mean), row.names = 1))
C2C12 <- normalizeBetweenArrays(C2C12, method="quantile")
C2C12 <- C2C12 - median(C2C12)

plot(C2C12)
abline(lm(C2C12[,2]~C2C12[,1]))

C2C12_stats <- cor.test(C2C12[,1], C2C12[,2], method='spearman', exact=F)
C2C12_stats


#=============================================================================================================================
# Differences in transcriptomic and proteomics in mouse Tissue
#=============================================================================================================================
MouseTissue_Prot <- MouseTissue_PXD000288
MouseTissue_mRNA <- mRNA[grepl('MouseTissue', colnames(mRNA))]
MouseTissue_mRNA <- data.frame(rownames(MouseTissue_mRNA), rowMeans(MouseTissue_mRNA))

MouseTissue <- merge(MouseTissue_Prot, MouseTissue_mRNA, by=1)
colnames(MouseTissue) <- c('SYMBOL', 'Proteome', 'Transcriptome')  
MouseTissue <- na.omit(data.frame(aggregate(MouseTissue[,2:3], b=list(MouseTissue$SYMBOL), FUN=mean), row.names = 1))
MouseTissue <- normalizeBetweenArrays(MouseTissue, method="quantile")
MouseTissue <- MouseTissue - median(MouseTissue)

plot(MouseTissue)
abline(lm(MouseTissue[,2]~MouseTissue[,1]))

MouseTissue_stats <- cor.test(MouseTissue[,1], MouseTissue[,2], method='spearman', exact=F)
MouseTissue_stats


#=============================================================================================================================
# Differences in transcriptomic and proteomics in rat L6
#=============================================================================================================================
RatL6_Prot <- RatL6_PXD001543
RatL6_mRNA <- mRNA[grepl('RatL6', colnames(mRNA))]
RatL6_mRNA <- data.frame(rownames(RatL6_mRNA), rowMeans(RatL6_mRNA))

RatL6 <- merge(RatL6_Prot, RatL6_mRNA, by=1)
colnames(RatL6) <- c('SYMBOL', 'Proteome', 'Transcriptome')  
RatL6 <- na.omit(data.frame(aggregate(RatL6[,2:3], b=list(RatL6$SYMBOL), FUN=mean), row.names = 1))
RatL6 <- normalizeBetweenArrays(RatL6, method="quantile")
RatL6 <- RatL6 - median(RatL6)

plot(RatL6)
abline(lm(RatL6[,2]~RatL6[,1]))

RatL6_stats <- cor.test(RatL6[,1], RatL6[,2], method='spearman', exact=F)
RatL6_stats


#=============================================================================================================================
# Differences in transcriptomic and proteomics in rat tissue
#=============================================================================================================================
RatTissue_Prot <- RatTissue_PXD007182
RatTissue_mRNA <- mRNA[grepl('RatTissue', colnames(mRNA))]
RatTissue_mRNA <- data.frame(rownames(RatTissue_mRNA), rowMeans(RatTissue_mRNA))

RatTissue <- merge(RatTissue_Prot, RatTissue_mRNA, by=1)
colnames(RatTissue) <- c('SYMBOL', 'Proteome', 'Transcriptome')  
RatTissue <- na.omit(data.frame(aggregate(RatTissue[,2:3], b=list(RatTissue$SYMBOL), FUN=mean), row.names = 1))
RatTissue <- normalizeBetweenArrays(RatTissue, method="quantile")
RatTissue <- data.frame(RatTissue - median(RatTissue))

plot(RatTissue)
abline(lm(RatTissue[,2]~RatTissue[,1]))

RatTissue_stats <- cor.test(RatTissue[,1], RatTissue[,2], method='spearman', exact=F)
RatTissue_stats



#=============================================================================================================================
# Diffferences in transcriptomic and proteomics in mouse Tissue
#=============================================================================================================================
png(filename=here("Figures", "Proteome_transcriptome_Mouse.png"), #print graph
    units="cm", width=24, height=12, 
    pointsize=12, res=300)
par(mfrow=c(1,2))

plot(C2C12, main="Mouse C2C12",
     sub=paste(nrow(C2C12), " proteins, ", "r=", signif(C2C12_stats$estimate, 2), ", p=", signif(C2C12_stats$p.value, 2), sep=""))
abline(lm(C2C12[,2]~C2C12[,1]), lwd=4, col="#D3C839")
abline(h=0, lty=3)
abline(v=0, lty=3)

plot(MouseTissue, main='Mouse Tissue',
     sub=paste(nrow(MouseTissue), " proteins, ", "r=", signif(MouseTissue_stats$estimate, 2), ", p=", signif(MouseTissue_stats$p.value, 2), sep=""))
abline(lm(MouseTissue[,2]~MouseTissue[,1]), lwd=4, col="#E69F00")
abline(h=0, lty=3)
abline(v=0, lty=3)
dev.off()




#=============================================================================================================================
# Diffferences in transcriptomic and proteomics in mouse and rat
#=============================================================================================================================
png(filename=here("Figures", "Proteome_transcriptome.png"), #print graph
    units="cm", width=24, height=24, 
    pointsize=12, res=300)
par(mfrow=c(2,2))

plot(C2C12, main="Mouse C2C12",
     sub=paste(nrow(C2C12), " proteins, ", "r=", signif(C2C12_stats$estimate, 2), ", p=", signif(C2C12_stats$p.value, 2), sep=""))
abline(lm(C2C12[,2]~C2C12[,1]), lwd=4, col="#D3C839")
abline(h=0, lty=3)
abline(v=0, lty=3)

plot(MouseTissue, main='Mouse Tissue',
     sub=paste(nrow(MouseTissue), " proteins, ", "r=", signif(MouseTissue_stats$estimate, 2), ", p=", signif(MouseTissue_stats$p.value, 2), sep=""))
abline(lm(MouseTissue[,2]~MouseTissue[,1]), lwd=4, col="#E69F00")
abline(h=0, lty=3)
abline(v=0, lty=3)

plot(RatL6, main='Rat L6',
     sub=paste(nrow(RatL6), " proteins, ", "r=", signif(RatL6_stats$estimate, 2), ", p=", signif(RatL6_stats$p.value, 2), sep=""))
abline(lm(RatL6[,2]~RatL6[,1]), lwd=4, col="#CC79A7")
abline(h=0, lty=3)
abline(v=0, lty=3)

plot(RatTissue, main='Rat Tissue',
     sub=paste(nrow(RatTissue), " proteins, ", "r=", signif(RatTissue_stats$estimate, 2), ", p=", signif(RatTissue_stats$p.value, 2), sep=""))
abline(lm(RatTissue[,2]~RatTissue[,1]), lwd=4, col="#D55E00")
abline(h=0, lty=3)
abline(v=0, lty=3)
dev.off()


#=============================================================================================================================
# Find "outliers": proteins and mRNA dispcrepencies
#=============================================================================================================================
C2C12_o <- C2C12[,1] / C2C12[,2]
C2C12_o <- C2C12_o[C2C12_o<0]
C2C12_o <- data.frame(C2C12[names(C2C12_o),])
C2C12_mRNAhigh <- C2C12_o[C2C12_o$Proteome<2 & C2C12_o$Transcriptome>2,]
write.table(C2C12_mRNAhigh, file="C2C12_mRNAhigh.txt", sep="\t")
C2C12_mRNAlow <- C2C12_o[C2C12_o$Proteome>2 & C2C12_o$Transcriptome<2,]
write.table(C2C12_mRNAlow, file="C2C12_mRNAlow.txt", sep="\t")
