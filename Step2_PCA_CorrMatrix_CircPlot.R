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


#--------------------------------------------------------------------------------------------------
# Figure: Primary Component Analysis
#--------------------------------------------------------------------------------------------------
noNA <- na.omit(res) # remove missing values for PCA
groups <- t(noNA) #transpose to compare groups, not genes
groups <- cbind(groups, Sample=Sample)

PCA_data = prcomp(t(noNA)) #transpose to compare groups, not genes
PCAbar <- plot(PCA_data)
pca <- autoplot(PCA_data, data=groups,
                colour='Sample', shape='Sample', fill='Sample',
                frame=TRUE, frame.type='convex') +
  theme_bw() +
  theme(plot.title  = element_text(face="bold", color="black", size=7, angle=0),
        axis.text.x = element_text(color="black", size=6, angle=0),
        axis.text.y = element_text(color="black", size=6, angle=0),
        axis.title  = element_text(face="bold", color="black", size=7, angle=0),
        legend.text = element_text(face="bold", color="black", size=6, angle=0),
        legend.position="right", legend.title = element_blank()) +
  scale_fill_manual(breaks=c("Human Myotube", "Mouse C2C12", "Rat L6", 
                             "Human Tissue", "Mouse Tissue", "Rat Tissue"),
                    values=c("#56B4E9", "#0072B2", "#D3C839", "#E69F00", "#CC79A7", "#D55E00")) +
  scale_color_manual(breaks=c("Human Myotube", "Mouse C2C12", "Rat L6", 
                              "Human Tissue", "Mouse Tissue", "Rat Tissue"),
                     values=c("#56B4E9", "#0072B2", "#D3C839", "#E69F00", "#CC79A7", "#D55E00")) +
  scale_shape_manual(breaks=c("Human Myotube", "Mouse C2C12", "Rat L6", 
                              "Human Tissue", "Mouse Tissue", "Rat Tissue"),
                     values=cbShapes)
pca


png(filename="../../Figures/PCA.png", 
    units="cm", width=9, height=5, 
    pointsize=12, res=300)
pca
dev.off()


png(filename="../../Figures/Annotated_PCA.png", 
    units="cm", width=19, height=16, 
    pointsize=12, res=300)
pca +  geom_text_repel(aes(label=colnames(noNA)), size=1, colour='black')
dev.off()


#--------------------------------------------------------------------------------------------------
# Figure: Correlation matrix
#--------------------------------------------------------------------------------------------------
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



#--------------------------------------------------------------------------------------------------
# Figure: Circular plot
#--------------------------------------------------------------------------------------------------
# load data and define samples
res <- readRDS("GENENAME_norm.Rds")
res_means <- data.frame(GENENAME=rownames(res),
                        HumanCell=rowMeans(res[grepl('HumanCell', colnames(res))]),
                        MouseC2C12=rowMeans(res[grepl('MouseC2C12', colnames(res))]),
                        RatL6=rowMeans(res[grepl('RatL6', colnames(res))]),
                        HumanTissue=rowMeans(res[grepl('HumanTissue', colnames(res))]),
                        MouseTissue=rowMeans(res[grepl('MouseTissue', colnames(res))]),
                        RatTissue=rowMeans(res[grepl('RatTissue', colnames(res))]))

# find top 100 genes for each study
signENSEMBL <- list(HumanCell <- (res_means[order(res_means$HumanCell, decreasing=T),][1:100,1]),
                    MouseC2C12 <- (res_means[order(res_means$MouseC2C12, decreasing=T),][1:100,1]),
                    RatL6 <- (res_means[order(res_means$RatL6, decreasing=T),][1:100,1]),
                    HumanTissue <- (res_means[order(res_means$HumanTissue, decreasing=T),][1:100,1]),
                    MouseTissue <- (res_means[order(res_means$MouseTissue, decreasing=T),][1:100,1]),
                    RatTissue <- (res_means[order(res_means$RatTissue, decreasing=T),][1:100,1])
)

# make a matrix of the genes common between groups
dat <- data.frame(row.names=c('Human Myotube', 'Mouse C2C12', 'Rat L6', 'Human Tissue', 'Mouse Tissue', 'Rat Tissue'))
for (y in 1:length(signENSEMBL)){
  study1 <- numeric()
  for (i in 1:length(signENSEMBL)){
    x <- c(signENSEMBL[[y]], signENSEMBL[[i]])
    x <- length(x[duplicated(x)])
    study1[i] <- x
  }
  dat[,y] <- study1
}
colnames(dat) <- rownames(dat)

# Set all self directions to 0
for (i in 1:length(signENSEMBL)){
  dat[i,i] <- 0
}
dat <- as.matrix(dat)

# plot
library(circlize)

png(filename="../../Figures/ChordPlot.png", #print graph
    units="cm", width=11, height=11, 
    pointsize=12, res=1200)

circos.clear()
circos.par(gap.after = 10)
chordDiagram(dat, symmetric=T, grid.col=cbPalette[1:6])

dev.off()

# re-organize to have human muscle on top
library(circlize)
dat2 <- dat[c(2,5,3,6,1,4),c(2,5,3,6,1,4)]
png(filename="../../Figures/ChordPlot2.png", #print graph
    units="cm", width=11, height=11, 
    pointsize=12, res=1200)

circos.clear()
circos.par(gap.after = 10)
chordDiagram(dat2, symmetric=T, grid.col=cbPalette[c(2,5,3,6,1,4)])

dev.off()




#--------------------------------------------------------------------------------------------------
# Figure: correlation contraction/proliferation
#--------------------------------------------------------------------------------------------------
Gene1 <- 'CKM'
Gene2 <- 'MKI67'
data  <- data.frame(as.numeric(res[Gene1,] ), as.numeric(res[Gene2,]))
data$Sample <- Sample
colnames(data) <- c("Gene1", "Gene2", "Sample")

Stats <- paste("Spearman", 
               "\nr =", signif(cor.test(data$Gene1, data$Gene2, method="spearman", exact=F)$estimate, digits=2),
               "\np =", signif(cor.test(data$Gene1, data$Gene2, method="spearman", exact=F)$p.value, digits=2), sep=" ")

#getting the convex hull of each unique point set
library(plyr)
df <- na.omit(data)
find_hull <- function(df) df[chull(df$Gene1, df$Gene2), ]
hulls <- ddply(df, "Sample", find_hull)

#plot
correlation <- ggplot(data, aes(x=Gene1, y=Gene2)) +
  geom_point(data=data, aes(x=Gene1, y=Gene2, colour=Sample, shape=Sample, fill=Sample)) +
  geom_polygon(data=hulls, aes(fill=Sample, colour=Sample), alpha=0.2) +
  labs(x=paste(Gene1, ", log2(expression)", sep=""),
       y=paste(Gene2, ", log2(expression)", sep=""),
       title="Proliferation vs Differentiation") +
  theme_bw() +
  theme(plot.title  = element_text(face="bold", color="black", size=7, angle=0),
        axis.text.x = element_text(size=6, angle=0, hjust = 1),
        axis.text.y = element_text(size=6, angle=0, hjust = 1),
        axis.title  = element_text(face="bold", color="black", size=7, angle=0),
        legend.text = element_text(face="bold", color="black", size=6, angle=0),
        legend.title = element_blank(),
        legend.position="right") +
  geom_smooth(method="lm", color= "black", formula= y ~ x, se = T) +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(values=cbPalette) +
  scale_shape_manual(values=cbShapes)
correlation

