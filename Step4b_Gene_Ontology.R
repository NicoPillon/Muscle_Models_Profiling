##############################################################################################################################
#      Expression data in rat, mouse and human muscle, L6, C2C12 and human myotubes                                  
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

