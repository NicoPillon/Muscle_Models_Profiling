setwd("C:/Dropbox/NICO/R/Across_Species/Data/Transcriptomics")
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
keytypes(org.Hs.eg.db) # use keytypes to list all supporting types
setwd("C:/Dropbox/NICO/R/Across_Species/Data/Transcriptomics")
rawdata <- readRDS("GENENAME_norm.Rds")
rawdata <- cbind(rawdata[grepl('HumanCell', colnames(rawdata))],
                 rawdata[grepl('HumanTissue', colnames(rawdata))],
                 rawdata[grepl('MouseC2C12', colnames(rawdata))],
                 rawdata[grepl('MouseTissue', colnames(rawdata))],
                 rawdata[grepl('RatL6', colnames(rawdata))],
                 rawdata[grepl('RatTissue', colnames(rawdata))])

#Annotate with ENTREZID
Annotation <- readRDS("C:/Dropbox/NICO/R/Annotation_All.Rds") # This is my homemade annotation file
rawdata <- merge(Annotation, rawdata, by.x=3, by.y=0, all=F) #convert to ENTREZID
rawdata <- rawdata[!duplicated(rawdata$ENTREZID),] #remove duplicated ENTREZID
rawdata <- rawdata[!is.na(rawdata$ENTREZID),] #remove duplicated ENTREZID
rawdata <- data.frame(rawdata[,6:ncol(rawdata)], row.names=rawdata$ENTREZID)

#=============================================================================================================================
# Calculate model-specific genes
#=============================================================================================================================
library(limma)
design <- model.matrix(~ 0+factor(c(rep(1 ,length(grep('HumanCell', colnames(rawdata)))),
                                    rep(2 ,length(grep('HumanTissue', colnames(rawdata)))),
                                    rep(3 ,length(grep('MouseC2C12', colnames(rawdata)))),
                                    rep(4 ,length(grep('MouseTissue', colnames(rawdata)))),
                                    rep(5 ,length(grep('RatL6', colnames(rawdata)))),
                                    rep(6 ,length(grep('RatTissue', colnames(rawdata)))))))
colnames(design) <- c("HumanCell", "HumanTissue","MouseC2C12", "MouseTissue", "RatL6", "RatTissue")

fit <- lmFit(rawdata, design)


# Function for threshold selection
Threshold <- function(x) {
  x1 <- data.frame(topTable(fit2, coef=1, adjust.method="fdr",n=Inf, confint=T, sort.by="none")) # A list of differential expressed in group2 versus group1
  x2 <- data.frame(topTable(fit2, coef=2, adjust.method="fdr",n=Inf, confint=T, sort.by="none")) # A list of differential expressed in group3 versus group1
  x3 <- data.frame(topTable(fit2, coef=3, adjust.method="fdr",n=Inf, confint=T, sort.by="none")) # A list of differential expressed in group4 versus group1
  x4 <- data.frame(topTable(fit2, coef=4, adjust.method="fdr",n=Inf, confint=T, sort.by="none")) # A list of differential expressed in group5 versus group1
  x5 <- data.frame(topTable(fit2, coef=5, adjust.method="fdr",n=Inf, confint=T, sort.by="none")) # A list of differential expressed in group6 versus group1
  x <- data.frame(x1$AveExpr,
                  x1$adj.P.Val,
                  x2$adj.P.Val,
                  x3$adj.P.Val,
                  x4$adj.P.Val,
                  x5$adj.P.Val,
                  row.names = rownames(x1))
  x <- x[x$x1.adj.P.Val < 0.001 &
         x$x2.adj.P.Val < 0.001 &
         x$x3.adj.P.Val < 0.001 &
         x$x4.adj.P.Val < 0.001 &
         x$x5.adj.P.Val < 0.001 ,]
  return(x)
}



#Genes enriched in Human Myotubes
contrast.matrix <- makeContrasts(HumanCell-MouseC2C12,
                                 HumanCell-MouseTissue,
                                 HumanCell-RatL6,
                                 HumanCell-RatTissue,
                                 HumanCell-MouseTissue,
                                 levels=design)
fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))
HumanCell <- Threshold(x)


#Genes enriched in L6
contrast.matrix <- makeContrasts(RatL6-MouseC2C12,
                                 RatL6-MouseTissue,
                                 RatL6-HumanCell,
                                 RatL6-HumanTissue,
                                 RatL6-RatTissue,
                                 levels=design)
fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))
RatL6 <- Threshold(x)


#Genes enriched in C2C12
contrast.matrix <- makeContrasts(MouseC2C12-RatL6,
                                 MouseC2C12-RatTissue,
                                 MouseC2C12-HumanCell,
                                 MouseC2C12-HumanTissue,
                                 MouseC2C12-MouseTissue,
                                 levels=design)
fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))
MouseC2C12 <- Threshold(x)


#Genes enriched in human tissue
contrast.matrix <- makeContrasts(HumanTissue-RatL6,
                                 HumanTissue-RatTissue,
                                 HumanTissue-HumanCell,
                                 HumanTissue-MouseC2C12,
                                 HumanTissue-MouseTissue,
                                 levels=design)
fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))
HumanTissue <- Threshold(x)


#Genes enriched in rat tissue
contrast.matrix <- makeContrasts(RatTissue-RatL6,
                                 RatTissue-HumanTissue,
                                 RatTissue-HumanCell,
                                 RatTissue-MouseC2C12,
                                 RatTissue-MouseTissue,
                                 levels=design)
fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))
RatTissue <- Threshold(x)


#Genes enriched in mouse tissue
contrast.matrix <- makeContrasts(MouseTissue-RatL6,
                                 MouseTissue-HumanTissue,
                                 MouseTissue-HumanCell,
                                 MouseTissue-MouseC2C12,
                                 MouseTissue-RatTissue,
                                 levels=design)
fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))
MouseTissue <- Threshold(x)


#Make matrix
HumanCell$entrez <- rownames(HumanCell)
HumanCell$group <- "Human Myotube"
HumanCell$direction <- "up"
HumanCell$direction[HumanCell$LogFC_mean <0] <- "down"

MouseC2C12$entrez <- rownames(MouseC2C12)
MouseC2C12$group <- "Mouse C2C12"
MouseC2C12$direction <- "up"
MouseC2C12$direction[MouseC2C12$LogFC_mean <0] <- "down"

RatL6$entrez <- rownames(RatL6)
RatL6$group <- "Rat L6"
RatL6$direction <- "up"
RatL6$direction[RatL6$LogFC_mean <0] <- "down"

HumanTissue$entrez <- rownames(HumanTissue)
HumanTissue$group <- "Human Tissue"
HumanTissue$direction <- "up"
HumanTissue$direction[HumanTissue$LogFC_mean <0] <- "down"

MouseTissue$entrez <- rownames(MouseTissue)
MouseTissue$group <- "Mouse Tissue"
MouseTissue$direction <- "up"
MouseTissue$direction[MouseTissue$LogFC_mean <0] <- "down"

RatTissue$entrez <- rownames(RatTissue)
RatTissue$group <- "Rat Tissue"
RatTissue$direction <- "up"
RatTissue$direction[RatTissue$LogFC_mean <0] <- "down"

mydf <- rbind(HumanCell[,7:9],
              MouseC2C12[,7:9],
              RatL6[,7:9],
              HumanTissue[,7:9],
              MouseTissue[,7:9],
              RatTissue[,7:9])


#Gene Ontology
bp <- compareCluster(entrez~group, data=mydf,
                     fun="enrichGO", ont = "BP",
                     OrgDb = org.Hs.eg.db, pvalueCutoff=0.05)
dotplot(bp, show=5)








#=============================================================================================================================
# Gene Ontology using ClusterProfiler https://github.com/GuangchuangYu/clusterProfiler/issues/32
#=============================================================================================================================
setwd("C:/Dropbox/NICO/R/Across_Species/Data/Transcriptomics")
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
keytypes(org.Hs.eg.db) # use keytypes to list all supporting types

res <- readRDS("GENENAME_norm.Rds")
res <- data.frame(GENENAME=rownames(res),
                  HumanCell=rowMeans(res[grepl('HumanCell', colnames(res))]),
                  MouseC2C12=rowMeans(res[grepl('MouseC2C12', colnames(res))]),
                  RatL6=rowMeans(res[grepl('RatL6', colnames(res))]),
                  HumanTissue=rowMeans(res[grepl('HumanTissue', colnames(res))]),
                  MouseTissue=rowMeans(res[grepl('MouseTissue', colnames(res))]),
                  RatTissue=rowMeans(res[grepl('RatTissue', colnames(res))]))


#Annotate with ENTREZID
Annotation <- read.delim("C:/Dropbox/NICO/R/Annotation_All.txt", sep="") # This is my homemade annotation file
res <- merge(Annotation, res, by.x=3, by.y=0, all=F) #convert to ENTREZID
res <- res[!duplicated(res$ENTREZID),] #remove duplicated ENTREZID


#find top genes for each study
HumanCell <- (res[order(res$HumanCell, decreasing=T),][1:100,])
MouseC2C12 <- (res[order(res$MouseC2C12, decreasing=T),][1:100,])
RatL6 <- (res[order(res$RatL6, decreasing=T),][1:100,])
HumanTissue <- (res[order(res$HumanTissue, decreasing=T),][1:100,])
MouseTissue <- (res[order(res$MouseTissue, decreasing=T),][1:100,])
RatTissue <- (res[order(res$RatTissue, decreasing=T),][1:100,])


#Make matrix
HumanCell$entrez <- HumanCell$ENTREZID
HumanCell$group <- "Human Myotube"
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

HumanTissue$entrez <- HumanTissue$ENTREZID
HumanTissue$group <- "Human Tissue"
HumanTissue$direction <- "up"
HumanTissue$direction[HumanTissue$LogFC_mean <0] <- "down"

MouseTissue$entrez <- MouseTissue$ENTREZID
MouseTissue$group <- "Mouse Tissue"
MouseTissue$direction <- "up"
MouseTissue$direction[MouseTissue$LogFC_mean <0] <- "down"

RatTissue$entrez <- RatTissue$ENTREZID
RatTissue$group <- "Rat Tissue"
RatTissue$direction <- "up"
RatTissue$direction[RatTissue$LogFC_mean <0] <- "down"

mydf <- rbind(HumanCell[,13:15],
              MouseC2C12[,13:15],
              RatL6[,13:15],
              HumanTissue[,13:15],
              MouseTissue[,13:15],
              RatTissue[,13:15])


#Gene Ontology
bp <- compareCluster(entrez~group, data=mydf,
                              fun="enrichGO", ont = "BP",
                              OrgDb = org.Hs.eg.db, pvalueCutoff=0.05)
dotplot(bp, show=3)
saveRDS(bp, "../Stats/enrichGO_cells+Tissues.Rds")


#=============================================================================================================================
#Gene Ontology plot and simplify
#=============================================================================================================================
setwd("C:/Dropbox/NICO/R/Across_Species/Data/Transcriptomics")
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
keytypes(org.Hs.eg.db) # use keytypes to list all supporting types
bp <- readRDS("../../Stats/enrichGO.Rds")
bp.result <- bp@compareClusterResult
dotplot(bp, x=~direction) + ggplot2::facet_grid(~group)

#save original GO
png(filename="../../Figures/enrichGO.png", #print graph
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
png(filename="../../Figures/enrichGO_simplified.png", #print graph
    units="cm", width=24, height=12, 
    pointsize=12, res=600)
dotplot(bp2, x=~direction, showCategory = 3, font.size=9) + ggplot2::facet_grid(~group)
dev.off()

