setwd("C:/ownCloud/Projects/MuscleModels/Data/Transcriptomics")
rawdata <- readRDS("GENENAME_norm.Rds")
rawdata <- cbind(rawdata[grepl('HumanCell', colnames(rawdata))],
                 rawdata[grepl('MouseC2C12', colnames(rawdata))],
                 rawdata[grepl('RatL6', colnames(rawdata))])

#Filter out unwanted genes: microRNA, long non coding RNA...
rawdata <- rawdata[!grepl("AC[0-9][0-9][0-9].*.",    rownames(rawdata)),] 
rawdata <- rawdata[!grepl("AABR",    rownames(rawdata)),] 



#=============================================================================================================================
# Calculate model-specific genes
#=============================================================================================================================
library(limma)
design <- model.matrix(~ 0+factor(c(rep(1 ,length(grep('HumanCell', colnames(rawdata)))),
                                    rep(3 ,length(grep('MouseC2C12', colnames(rawdata)))),
                                    rep(5 ,length(grep('RatL6', colnames(rawdata)))))))
colnames(design) <- c("HumanCell", "MouseC2C12", "RatL6")

fit <- lmFit(rawdata, design)

# Function for threshold selection
Threshold <- function(x) {
  x1 <- data.frame(topTable(fit2, coef=1, adjust.method="fdr",n=Inf, confint=T, sort.by="none")) # A list of differential expressed in group2 versus group1
  x2 <- data.frame(topTable(fit2, coef=2, adjust.method="fdr",n=Inf, confint=T, sort.by="none")) # A list of differential expressed in group2 versus group1
  x <- data.frame(x1$AveExpr, 
                  x1$logFC, x1$adj.P.Val,
                  x2$logFC, x2$adj.P.Val,
                  LogFC_mean=(x1$logFC + x2$logFC)/2,
                  row.names = rownames(x1))
  x <- x[x1$AveExpr > 1 &
         x$x1.adj.P.Val < 0.01 & 
         x$x2.adj.P.Val < 0.01 & 
         x$x1.logFC*x$x2.logFC > 10 ,] #Select genes where both groups go in the same direction
  x <- na.omit(x)
  return(x)
}



#Genes enriched in Human Myotubes
contrast.matrix <- makeContrasts(HumanCell-MouseC2C12,
                                 HumanCell-RatL6,
                                 levels=design)
fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))
HumanCell <- Threshold(x)
colnames(HumanCell) <- c('AveExpr',
                         'LogFC_MouseC2C12', 'FDR_MouseC2C12',
                         'LogFC_RatL6', 'FDR_RatL6', 'LogFC_mean')
HumanCell_up <- HumanCell[order(HumanCell$LogFC_mean, decreasing=T),][1:1,c(1,6)]
HumanCell_down <- HumanCell[order(HumanCell$LogFC_mean, decreasing=F),][1:1,c(1,6)]


#Genes enriched in L6
contrast.matrix <- makeContrasts(RatL6-MouseC2C12,
                                 RatL6-HumanCell,
                                 levels=design)
fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))
RatL6 <- Threshold(x)
colnames(RatL6) <- c('AveExpr',
                     'LogFC_MouseC2C12', 'FDR_MouseC2C12',
                     'LogFC_HumanCell', 'FDR_HumanCell', 'LogFC_mean')
RatL6_up <- RatL6[order(RatL6$LogFC_mean, decreasing=T),][1:1,c(1,6)]
RatL6_down <- RatL6[order(RatL6$LogFC_mean, decreasing=F),][1:1,c(1,6)]



#Genes enriched in C2C12
contrast.matrix <- makeContrasts(MouseC2C12-RatL6,
                                 MouseC2C12-HumanCell,
                                 levels=design)
fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))
MouseC2C12 <- Threshold(x)
colnames(MouseC2C12) <- c('AveExpr',
                          'LogFC_RatL6', 'FDR_RatL6',
                          'LogFC_HumanCell', 'FDR_HumanCell', 'LogFC_mean')
MouseC2C12_up <- MouseC2C12[order(MouseC2C12$LogFC_mean, decreasing=T),][1:1,c(1,6)]
MouseC2C12_down <- MouseC2C12[order(MouseC2C12$LogFC_mean, decreasing=F),][1:1,c(1,6)]


TestPCR <- rbind(HumanCell_up, HumanCell_down,
                 RatL6_up, RatL6_down,
                 MouseC2C12_up, MouseC2C12_down)

TestPCR <- c(rownames(TestPCR),
             'MKI67', 'MYH7', 'MYH4',
             'SLC2A1', 'SLC2A4',
             'GYS1', 'GYS2', 'GSK3A', 'GSK3B',
             'CPT1A', 'CPT1B', 'CPT2',
             'TBP', 'HPRT1', 'B2M')
TestPCR
