##############################################################################################################################
#      Expression data in rat, mouse and human muscle, L6, C2C12 and human myotubes                                  
##############################################################################################################################
# This file aims at finding genes that thare significantly enriched in on the the cell model compared to the other two
# For instance genes that would be enriched in L6 compared to C2C12 and HSMC at FDR<0.05
library(here)
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(gplots)
library(stringr)
library(grid)
library(gridExtra)
#cbPalette <- c("#E69F00", "#0072B2", "#CC79A7", "#009E73", "#D3C839", "#BC5300", "#84C4E8", "#000000") #color palette for colorblind people
cbPalette <- c("#56B4E9", "#D3C839", "#CC79A7", "#0072B2", "#E69F00", "#D55E00") #color palette for colorblind people
#             light blue    yellow     pink     dark blue    orange      red
cbShapes  <- c(   21    ,    21    ,    24    ,    24    ,    22    ,    22    ,    23    ,     23   )
cbLines   <- c(   'a'   ,    'b'   ,    'c'   ,    'd'   ,    'e'   ,    'f'   ,    'g'   ,     'h'  )

res <- readRDS(here("Data_Processed", "GENENAME_norm.Rds"))
res <- cbind(res[grepl('HumanCell', colnames(res))],
                 res[grepl('MouseC2C12', colnames(res))],
                 res[grepl('RatL6', colnames(res))])
Sample <- c(rep("Human Myotube",  length(grep('HumanCell', colnames(res)))),
            rep("Mouse C2C12",   length(grep('MouseC2C12',  colnames(res)))),
            rep("Rat L6",        length(grep('RatL6',       colnames(res)))))

#=============================================================================================================================
# Calculate model-specific genes
#=============================================================================================================================
library(limma)
design <- model.matrix(~ 0+factor(c(rep(1 ,length(grep('HumanCell', colnames(res)))),
                                    rep(3 ,length(grep('MouseC2C12', colnames(res)))),
                                    rep(5 ,length(grep('RatL6', colnames(res)))))))
colnames(design) <- c("HumanCell", "MouseC2C12", "RatL6")

fit <- lmFit(res, design)

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


write.table(HumanCell, file=here("Stats", "HumanCell_Specific.txt"))
write.table(RatL6, file=here("Stats", "RatL6_Specific.txt"))
write.table(MouseC2C12, file=here("Stats", "MouseC2C12_Specific.txt"))



#======================================================================================================
# Figure: boxplots for top specific genes
#======================================================================================================
MouseC2C12 <- read.table(here("Stats", "MouseC2C12_Specific.txt"))
RatL6 <- read.table(here("Stats", "RatL6_Specific.txt"))
HumanCell <- read.table(here("Stats", "HumanCell_Specific.txt"))

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


# function to plot boxplots
PlotFunction <- function(genename) {
  data   <- data.frame()                                #create an empty dataframe to collect data
  for( i in 1:length(genename)) { 
    y     <- as.numeric(res[genename[i],])              #collect data for gene name i
    datay <- cbind.data.frame(Sample, y, rep(genename[i]))   #create table with x="sample type", y="data", "gene name"
    colnames(datay) <- c("x","y","Gene")                #rename column names to make it possible to rbind later
    data  <- rbind.data.frame(data, datay)              #bind the new gene data at the bottom of the previous one
  }
  data$x <- factor(data$x, levels=c("Human Myotube", "Mouse C2C12",  "Rat L6", 
                                    "Human Tissue",  "Mouse Tissue", "Rat Tissue")) #for a box plot, x should be a factor
  ggplot(data, aes(x=x, y=y, fill=x)) +  
    geom_boxplot() +
    labs(x="",
         y=paste(genename, "expression, log2")) +
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


png(filename=here("Figures", "TopGenes.png"), #print graph
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



##############################################################################################################################
# Gene Ontology                               
##############################################################################################################################
# Gene Ontology using ClusterProfiler https://github.com/GuangchuangYu/clusterProfiler/issues/32
# using the genes enriches in each model, we run gene enrichment to look at pathways specific to each cell model
library(here)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
keytypes(org.Hs.eg.db) # use keytypes to list all supporting types
MouseC2C12 <- read.table(here("Stats", "MouseC2C12_Specific.txt"))
RatL6 <- read.table(here("Stats", "RatL6_Specific.txt"))
HumanCell <- read.table(here("Stats", "HumanCell_Specific.txt"))

#genes have to be ENTREZID!
Annotation <- bitr(c(rownames(MouseC2C12), rownames(RatL6), rownames(HumanCell)),
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

MouseC2C12 <- merge(Annotation, MouseC2C12, by.x=1, by.y=0, all=F)
RatL6 <- merge(Annotation, RatL6, by.x=1, by.y=0, all=F)
HumanCell <- merge(Annotation, HumanCell, by.x=1, by.y=0, all=F)

#Make matrix
HumanCell$entrez <- HumanCell$ENTREZID
HumanCell$group <- "Human Primary"
HumanCell$direction <- "up"
HumanCell$direction[HumanCell$LogFC_mean <0] <- "down"
MouseC2C12$entrez <- MouseC2C12$ENTREZID
MouseC2C12$group <- "Mouse C2C12"
MouseC2C12$direction <- "up"
MouseC2C12$direction[MouseC2C12$LogFC_mean <0] <- "down"
RatL6$entrez <- RatL6$ENTREZID
RatL6$group <- "Rat L6"
RatL6$direction <- "up"
RatL6$direction[RatL6$LogFC_mean <0] <- "down"
mydf <- rbind(HumanCell[,9:11], MouseC2C12[,9:11], RatL6[,9:11])
rm(Annotation, HumanCell, MouseC2C12, RatL6)

#Gene Ontology
bp <- compareCluster(entrez~group+direction, data=mydf,
                     fun="enrichGO", ont = "BP",
                     OrgDb = org.Hs.eg.db, pvalueCutoff=0.05)
saveRDS(bp, here("Stats", "enrichGO.Rds"))


#=============================================================================================================================
#Gene Ontology plot and simplify
#=============================================================================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
keytypes(org.Hs.eg.db) # use keytypes to list all supporting types
bp <- readRDS(here("Stats", "enrichGO.Rds"))
bp.result <- bp@compareClusterResult
dotplot(bp, x=~direction) + ggplot2::facet_grid(~group)

#save original GO
png(filename=here("Figures", "enrichGO.png"), #print graph
    units="cm", width=24, height=18, 
    pointsize=12, res=600)
dotplot(bp, x=~direction, showCategory = 5, font.size=9) + ggplot2::facet_grid(~group)
dev.off()

#manually drop redundant GO to simplify
bp2 <- dropGO(bp, level=NULL, term=c('GO:0009167', 'GO:0009126', 'GO:0009205', 'GO:0009199', #purine metabolism
                                     'GO:0000398', 'GO:0000377', 'GO:0000375', #RNA splicing
                                     'GO:0033141', #STAT
                                     'GO:0003206', 'GO:0003007',#cardiac
                                     'GO:0000070', 'GO:0140014', #mitotic
                                     'GO:0060541') #respiratory
)
dotplot(bp2, x=~direction, show=3) + ggplot2::facet_grid(~group)

#save cleaned GO
png(filename=here("Figures", "enrichGO_simplified.png"), #print graph
    units="cm", width=24, height=12, 
    pointsize=12, res=600)
dotplot(bp2, x=~direction, showCategory = 3, font.size=9) + ggplot2::facet_grid(~group)
dev.off()


