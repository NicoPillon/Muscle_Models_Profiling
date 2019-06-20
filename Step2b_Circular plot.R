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
cbPalette <- c("#56B4E9", "#D3C839", "#CC79A7", "#0072B2", "#E69F00", "#D55E00") #color palette for colorblind people
#             light blue    yellow     pink     dark blue    orange      red
cbShapes  <- c(   21    ,    21    ,    24    ,    24    ,    22    ,    22    ,    23    ,     23   )
cbLines   <- c(   'a'   ,    'b'   ,    'c'   ,    'd'   ,    'e'   ,    'f'   ,    'g'   ,     'h'  )
setwd("C:/ownCloud/Projects/MuscleModels/Data/Transcriptomics")

# load data and define samples
res <- readRDS("GENENAME_norm.Rds")
res <- data.frame(GENENAME=rownames(res),
                  HumanCell=rowMeans(res[grepl('HumanCell', colnames(res))]),
                  MouseC2C12=rowMeans(res[grepl('MouseC2C12', colnames(res))]),
                  RatL6=rowMeans(res[grepl('RatL6', colnames(res))]),
                  HumanTissue=rowMeans(res[grepl('HumanTissue', colnames(res))]),
                  MouseTissue=rowMeans(res[grepl('MouseTissue', colnames(res))]),
                  RatTissue=rowMeans(res[grepl('RatTissue', colnames(res))]))

# find top 100 genes for each study
signENSEMBL <- list(HumanCell <- (res[order(res$HumanCell, decreasing=T),][1:100,1]),
                    MouseC2C12 <- (res[order(res$MouseC2C12, decreasing=T),][1:100,1]),
                    RatL6 <- (res[order(res$RatL6, decreasing=T),][1:100,1]),
                    HumanTissue <- (res[order(res$HumanTissue, decreasing=T),][1:100,1]),
                    MouseTissue <- (res[order(res$MouseTissue, decreasing=T),][1:100,1]),
                    RatTissue <- (res[order(res$RatTissue, decreasing=T),][1:100,1])
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

