setwd("C:/ownCloud/Projects/MuscleModels/Data/Transcriptomics")
rawdata <- readRDS("GENENAME_norm.Rds")
rawdata <- cbind(rawdata[grepl('HumanCell', colnames(rawdata))],
                 rawdata[grepl('MouseC2C12', colnames(rawdata))],
                 rawdata[grepl('RatL6', colnames(rawdata))])


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
                  row.names = rownames(x1))
  x <- x[x$x1.adj.P.Val < 0.05 & 
         x$x2.adj.P.Val < 0.05 & 
         x$x1.logFC*x$x2.logFC > 2 ,] #Select genes where both groups go in the same direction
  x <- na.omit(x)
  return(x)
}

#Genes enriched in Human Myotubes
contrast.matrix <- makeContrasts(HumanCell-MouseC2C12,
                                 HumanCell-RatL6,
                                 levels=design)
fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))
HumanCell <- Threshold(x)
HumanCell <- HumanCell[order(HumanCell$x1.AveExpr, decreasing=T),]#[1:50,]
colnames(HumanCell) <- c('AveExpr',
                         'LogFC_MouseC2C12', 'FDR_MouseC2C12',
                         'LogFC_RatL6', 'FDR_RatL6')
HumanCell$LogFC_mean <- HumanCell$LogFC_MouseC2C12 + HumanCell$LogFC_RatL6

#Genes enriched in L6
contrast.matrix <- makeContrasts(RatL6-MouseC2C12,
                                 RatL6-HumanCell,
                                 levels=design)
fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))
RatL6 <- Threshold(x)
RatL6 <- RatL6[order(RatL6$x1.AveExpr, decreasing=T),]#[1:50,]
colnames(RatL6) <- c('AveExpr',
                         'LogFC_MouseC2C12', 'FDR_MouseC2C12',
                         'LogFC_HumanCell', 'FDR_HumanCell')
RatL6$LogFC_mean <- RatL6$LogFC_MouseC2C12 + RatL6$LogFC_HumanCell

#Genes enriched in C2C12
contrast.matrix <- makeContrasts(MouseC2C12-RatL6,
                                 MouseC2C12-HumanCell,
                                 levels=design)
fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))
MouseC2C12 <- Threshold(x)
MouseC2C12 <- MouseC2C12[order(MouseC2C12$x1.AveExpr, decreasing=T),]#[1:50,]
colnames(MouseC2C12) <- c('AveExpr',
                     'LogFC_RatL6', 'FDR_RatL6',
                     'LogFC_HumanCell', 'FDR_HumanCell')
MouseC2C12$LogFC_mean <- MouseC2C12$LogFC_HumanCell + MouseC2C12$LogFC_RatL6


write.table(HumanCell, file="../Stats/HumanCell_Specific.txt", sep="\t")
write.table(RatL6, file="../Stats/RatL6_Specific.txt", sep="\t")
write.table(MouseC2C12, file="../Stats/MouseC2C12_Specific.txt", sep="\t")



#======================================================================================================
# Figure: boxplots for top specific genes
#======================================================================================================
setwd("C:/ownCloud/Projects/MuscleModels/Data/Transcriptomics")
MouseC2C12 <- read.table("Stats/MouseC2C12_Specific.txt")
RatL6 <- read.table("Stats/RatL6_Specific.txt")
HumanCell <- read.table("Stats/HumanCell_Specific.txt")

#genes that should be excluded
RatL6 <- RatL6[!(rownames(RatL6) %in% 'GM24698'),]
RatL6 <- RatL6[!(rownames(RatL6) %in% 'AC111885.1'),]
RatL6 <- RatL6[!(rownames(RatL6) %in% 'AC104389.4'),]
MouseC2C12 <- MouseC2C12[!(rownames(MouseC2C12) %in% 'AABR07016578.1'),]

#select top 10 genes
TopHuman <- rownames(HumanCell[order(HumanCell$LogFC_mean, decreasing=T),])[1:10]
TopMouse <- rownames(MouseC2C12[order(MouseC2C12$LogFC_mean, decreasing=T),])[1:10]
TopRat   <- rownames(RatL6[order(RatL6$LogFC_mean, decreasing=T),])[1:10]

TopHuman
TopMouse
TopRat

res <- readRDS("GENENAME_norm.Rds")
res <- cbind(res[grep('HumanCell', colnames(res))],
             res[grep('HumanTissue', colnames(res))],
             res[grep('MouseC2C12', colnames(res))],
             res[grep('MouseTissue', colnames(res))],
             res[grep('RatL6', colnames(res))],
             res[grep('RatTissue', colnames(res))])
Sample <- c(rep("Human Primary",  length(grep('HumanCell', colnames(res)))),
            rep("Human Tissue", length(grep('HumanTissue',   colnames(res)))),  #list of sample types
            rep("Mouse C2C12",   length(grep('MouseC2C12',  colnames(res)))),
            rep("Mouse Tissue",   length(grep('MouseTissue',  colnames(res)))),
            rep("Rat L6",        length(grep('RatL6',       colnames(res)))),
            rep("Rat Tissue",   length(grep('RatTissue',  colnames(res)))))

# function to plot boxplots
library(ggplot2)
PlotFunction <- function(genename) {
  data   <- data.frame()                                #create an empty dataframe to collect data
  for( i in 1:length(genename)) { 
    y     <- as.numeric(res[genename[i],])              #collect data for gene name i
    datay <- cbind.data.frame(Sample, y, rep(genename[i]))   #create table with x="sample type", y="data", "gene name"
    colnames(datay) <- c("x","y","Gene")                #rename column names to make it possible to rbind later
    data  <- rbind.data.frame(data, datay)              #bind the new gene data at the bottom of the previous one
  }
  data$x <- factor(data$x, levels=c("Human Primary", "Mouse C2C12",  "Rat L6", 
                                    "Human Tissue",  "Mouse Tissue", "Rat Tissue")) #for a box plot, x should be a factor
  ggplot(data, aes(x=x, y=y, fill=x)) +  
    geom_boxplot() +
    labs(x="",
         y=paste(genename, "expression, log2")) +
    theme(plot.title  = element_text(face="bold", color="black", size=7, angle=0),
          axis.text.x = element_text(color="black", size=6, angle=45, hjust=1),
          axis.text.y = element_text(color="black", size=6, angle=0),
          axis.title  = element_text(face="bold", color="black", size=7, angle=0),
          legend.text = element_text(face="bold", color="black", size=6, angle=0),
          legend.position="none", legend.title = element_blank()) +
    scale_fill_manual(values=c("#E69F00", "#CC79A7", "#D3C839", 
                               "#0072B2", "#009E73", "#BC5300")) +
    scale_color_manual(values=c("#E69F00", "#CC79A7", "#D3C839", 
                                "#0072B2", "#009E73", "#BC5300"))
}


library(grid)
library(gridExtra)
png(filename="../Figures/TopGenes.png", #print graph
    units="cm", width=18, height=18, 
    pointsize=12, res=300)
matrix <- rbind(c(1,2, 3, 4),
                c(5,6, 7, 8),
                c(9,10,11,12))
grid.arrange(PlotFunction(TopHuman[1]), PlotFunction(TopHuman[2]),
             PlotFunction(TopHuman[3]), PlotFunction(TopHuman[4]), 
             PlotFunction(TopMouse[1]), PlotFunction(TopMouse[2]), 
             PlotFunction(TopMouse[3]), PlotFunction(TopMouse[4]), 
             PlotFunction(TopRat[1]),   PlotFunction(TopRat[2]),   
             PlotFunction(TopRat[3]),   PlotFunction(TopRat[4]),  
             layout_matrix=matrix)
dev.off()
