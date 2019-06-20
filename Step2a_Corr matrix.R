#===============================================================================================
# Expression data in rat, mouse and human muscle, L6, C2C12 and human myotubes
#===============================================================================================
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(gplots)
library(stringr)
library(grid)
library(gridExtra)
cbPalette <- c("#E69F00", "#0072B2", "#CC79A7", "#009E73", "#D3C839", "#BC5300", "#84C4E8", "#000000") #color palette for colorblind people
cbShapes  <- c(   21    ,    21    ,    24    ,    24    ,    22    ,    22    ,    23    ,     23   )
cbLines   <- c(   'a'   ,    'b'   ,    'c'   ,    'd'   ,    'e'   ,    'f'   ,    'g'   ,     'h'  )
setwd("C:/ownCloud/Projects/MuscleModels/Data/Transcriptomics")

#load data and define samples
res <- readRDS("GENENAME_norm.Rds")

# Impute missing values with KNN imputation
library(impute)
noNA <- res[rowSums(is.na(res)) / ncol(res) < 10/100, ]
noNA <- (as.matrix((noNA)))
noNA <- impute.knn(noNA, k=100, rowmax=0.8, colmax=1, maxp=nrow(noNA))
noNA  <- (data.frame(noNA$data))

# calculate mean for each group
noNAmean <- data.frame(HumanCell=rowMeans(noNA[grepl('HumanCell', colnames(noNA))]),
                       HumanTissue=rowMeans(noNA[grepl('HumanTissue', colnames(noNA))]),
                       MouseC2C12=rowMeans(noNA[grepl('MouseC2C12', colnames(noNA))]),
                       MouseTissue=rowMeans(noNA[grepl('MouseTissue', colnames(noNA))]),
                       RatL6=rowMeans(noNA[grepl('RatL6', colnames(noNA))]),
                       RatTissue=rowMeans(noNA[grepl('RatTissue', colnames(noNA))]))
colnames(noNAmean) <- c('Human Myotube', 'Human Tissue', 'Mouse C2C12', 'Mouse Tissue', 'Rat L6', 'Rat Tissue')

# make correlation matrix
cormatrix <- cor(noNAmean, method="spearman")
cormatrix[cormatrix==1]

# plot
library(RColorBrewer)
Colors=c("white", "red")
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(99)

png(filename="../../Figures/Figure_CorrelationMatrix.png", #print graph
    units="cm", width=15, height=15, 
    pointsize=12, res=1200)
heatmap.2(cormatrix, scale="none",
          margins =c(11,11), # adjust margins to row names
          key=T, keysize=2,                                  # remove color key legend
          density.info="none",                  # turns off density plot inside color legend
          trace="none",
          col=Colors,
          cexRow=1.5, cexCol=1.5,
          breaks = seq(0.6, 0.9, length.out = 100))
dev.off()

