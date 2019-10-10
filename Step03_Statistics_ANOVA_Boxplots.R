#===============================================================================================
# Expression data in rat, mouse and human muscle, L6, C2C12 and human myotubes
#===============================================================================================
# This file assumes that all data has been processed and normalized as in "Step1_Data.Rds"
# and is available as one single file: "GENENAME_norm.Rds" 
library(here)
library(limma)
library(impute)
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(ggpubr)
library(gplots)
library(car)


#=============================================================================================================================
# Genes different within species
#=============================================================================================================================
# Here we want to know how many genes are different in each cell model compared
# to muscle tissue. So we perform intra-species comparisons: L6 vs Rat tissue, C2C12
# vs mouse tissue and human primary vs human tissue. Simple t-test is used here with limma.
library(limma)
muscle <- readRDS(here("Data_processed", "GENENAME_norm.Rds"))
muscle <- muscle[!grepl('HEK', colnames(muscle))] #remove HEK cells
muscle <- muscle[!grepl('HeLa', colnames(muscle))] #remove HeLa cells

design <- model.matrix(~ 0+factor(c(rep(1, length(grep('HumanCell', colnames(muscle)))),
                                    rep(2, length(grep('HumanTissue', colnames(muscle)))),
                                    rep(3, length(grep('MouseC2C12', colnames(muscle)))),
                                    rep(4, length(grep('MouseTissue', colnames(muscle)))),
                                    rep(5, length(grep('RatL6', colnames(muscle)))),
                                    rep(6, length(grep('RatTissue', colnames(muscle)))))))
colnames(design) <- c("HumanCell", "HumanTissue", "MouseC2C12","MouseTissue", "RatL6", "RatTissue")
fit <- lmFit(muscle, design)

contrast.matrix <- makeContrasts(HumanCell-HumanTissue,  # 1 - comparison human
                                 MouseC2C12-MouseTissue, # 2 - comparison mouse
                                 RatL6-RatTissue,        # 3 - comparison rat
                                 levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)

Human <- data.frame(topTable(fit2, coef=1, adjust.method="fdr", n=Inf)) # 1 - comparison human
Mouse <- data.frame(topTable(fit2, coef=2, adjust.method="fdr", n=Inf)) # 2 - comparison mouse
Rat   <- data.frame(topTable(fit2, coef=3, adjust.method="fdr", n=Inf)) # 3 - comparison rat

#Select significant genes at FDR<0.05
Human <- Human[Human$adj.P.Val < 0.01,]
Mouse <- Mouse[Mouse$adj.P.Val < 0.01,]
Rat   <- Rat[Rat$adj.P.Val < 0.01,]

#How many genes are significant?
nrow(Human)*100/nrow(muscle) # percentage for human
nrow(Mouse)*100/nrow(muscle) # percentage for mouse
nrow(Rat)*100/nrow(muscle)   # percentage for rat

#Venn Diagram
require(VennDiagram)
venn.diagram(x=list(A=rownames(Human), B=rownames(Mouse), C=rownames(Rat)), filename=here("Figures", "VennDiagram.tiff"),
             units="cm", width=9, height=5, res=1200, margin=0.05,
             col = "black", lwd=0.5,
             fill = c("#0072B2", "#009E73", "#D55E00"), 
             alpha = 0.50,
             fontfamily="sans", cex = 0.6,
             category.names=c('Human', 'Mouse', 'Rat'),
             cat.col = c("#0072B2", "#009E73", "#D55E00"), cat.cex=0.7, cat.fontface="bold", cat.fontfamily="sans",
             cat.pos = c(0, 0, 180),
             cat.dist= c(0.05, 0.05, 0.05))


#=============================================================================================================================
# Calculate 2-way ANOVA
#=============================================================================================================================
# This is to try to find profiles of genes where there is an effect of species but not samples,
# an effect of cell vs tissue or interaction. Note that the 2-way ANOVA here is not valid for 
# the genes with non-uniform variance and can only be used as an approximate tool to find profiles.
library(multcomp)
library(car)
muscle <- readRDS(here("Data_Processed", "GENENAME_batch.Rds"))
muscle <- muscle[!grepl('HEK', colnames(muscle))] #remove HEK cells
muscle <- muscle[!grepl('HeLa', colnames(muscle))] #remove HeLa cells

#list of sample types
design <- c(rep(1, length(grep('HumanCell',   colnames(muscle)))), 
            rep(2, length(grep('HumanTissue', colnames(muscle)))),
            rep(3, length(grep('MouseC2C12',  colnames(muscle)))),
            rep(4, length(grep('MouseTissue', colnames(muscle)))),
            rep(5, length(grep('RatL6',       colnames(muscle)))),
            rep(6, length(grep('RatTissue',   colnames(muscle)))))
i <- 1
#prepare data in proper table
stats <- data.frame(Sample=numeric(), Species=numeric(), Interaction=numeric(), Residuals=numeric())
for (i in 1:nrow(muscle)) {
  testdata <- data.frame(t(muscle[i,]))
  colnames(testdata) <- 'Gene'
  #Make a new column "Species" with only species
  Species <- rownames(testdata)
  Species<-gsub("C.*","",Species)
  Species<-gsub("HeLa.*","",Species)
  Species<-gsub("HEK.*","",Species)
  Species<-gsub("Tissue.*","",Species)
  Species<-gsub("L6.*","",Species)
  testdata$Species <- as.factor(Species)
  rm(Species)
  #Make a new column "Sample" with either Cell or issue
  Sample <- rownames(testdata)
  Sample<-gsub("Human","",Sample)
  Sample<-gsub("Mouse","",Sample)
  Sample<-gsub("Rat","",Sample)
  Sample<-gsub("_.*","",Sample)
  Sample<-gsub("L6","Cell",Sample)
  Sample<-gsub("C2C12","Cell",Sample)
  testdata$Sample <- as.factor(Sample)
  rm(Sample)
  #two-way ANOVA test in R for unbalanced designs
  my_anova <- aov(Gene ~ Sample*Species, data=testdata)
  my_anova <- Anova(my_anova, type = "II")
  pval <- my_anova$`Pr(>F)`
  stats <- rbind(stats, pval)
  colnames(stats) <- c('Sample.pval', 'Species.pval', 'Interaction.pval', 'Residuals')
}
rownames(stats) <- rownames(muscle)

#calculate Mean, Sd and n
coeff <- data.frame(HumanCell_Mean = rowMeans(muscle[,design==1], na.rm=T),
                    HumanCell_Sd = rowSds(as.matrix(muscle[,design==1]), na.rm=T),
                    HumanCell_n = rowSums(!is.na(muscle[,design==1])),
                    HumanTissue_Mean = rowMeans(muscle[,design==2], na.rm=T),
                    HumanTissue_Sd   = rowSds(as.matrix(muscle[,design==2]), na.rm=T),
                    HumanTissue_n    = rowSums(!is.na(muscle[,design==2])),
                    MouseC2C12_Mean = rowMeans(muscle[,design==3], na.rm=T),
                    MouseC2C12_Sd   = rowSds(as.matrix(muscle[,design==3]), na.rm=T),
                    MouseC2C12_n    = rowSums(!is.na(muscle[,design==3])),
                    MouseTissue_Mean = rowMeans(muscle[,design==4], na.rm=T),
                    MouseTissue_Sd   = rowSds(as.matrix(muscle[,design==4]), na.rm=T),
                    MouseTissue_n    = rowSums(!is.na(muscle[,design==4])),
                    RatL6_Mean = rowMeans(muscle[,design==5], na.rm=T),
                    RatL6_Sd   = rowSds(as.matrix(muscle[,design==5]), na.rm=T),
                    RatL6_n    = rowSums(!is.na(muscle[,design==5])),
                    RatTissue_Mean = rowMeans(muscle[,design==6], na.rm=T),
                    RatTissue_Sd   = rowSds(as.matrix(muscle[,design==6]), na.rm=T),
                    RatTissue_n    = rowSums(!is.na(muscle[,design==6])))

#merge all stats
stats_all <- cbind(coeff,
                   stats[,1:3],
                   Sample.FDR=p.adjust(stats$Sample, method="fdr"),
                   Species.FDR=p.adjust(stats$Species, method="fdr"),
                   Interaction.FDR=p.adjust(stats$Interaction, method="fdr"))
  
write.table(stats_all, file=here("Stats", "2-way_ANOVA.txt"))



#------------------------------------------------------------------------------------------------------------------
# Check individual genes
#------------------------------------------------------------------------------------------------------------------
# This code for a single gene name will plot the data and allows the calculation of the homogeneity
# of variances to check if the use of the 2-way ANOVA is valid for that specific gene.
genename <- 'CR1L'

{testdata <- data.frame(t(muscle[genename,]))
  colnames(testdata) <- 'Gene'
  t1 <- rownames(testdata)
  t1<-gsub("C.*","",t1)
  t1<-gsub("HeLa.*","",t1)
  t1<-gsub("HEK.*","",t1)
  t1<-gsub("Tissue.*","",t1)
  t1<-gsub("L6.*","",t1)
  testdata$Species <- as.factor(t1)
  rm(t1)
  t2 <- rownames(testdata)
  t2<-gsub("Human","",t2)
  t2<-gsub("Mouse","",t2)
  t2<-gsub("Rat","",t2)
  t2<-gsub("_.*","",t2)
  t2<-gsub("L6","Cell",t2)
  t2<-gsub("C2C12","Cell",t2)
  testdata$Sample <- as.factor(t2)
  rm(t2)
}

#plot data
ggboxplot(testdata, x="Species", y="Gene", color="Sample", xlab='', ylab=genename)

#2-way ANOVA
res.aov <- aov(Gene ~ Species * Sample, data=testdata)
summary(res.aov)

#check residuals and normality - detect outliers
par(mfrow=c(1,2))
plot(res.aov, 1)
plot(res.aov, 2)

#check homogeneity of variances: if p<0.05 -> violation of asumption = variances not homogeneous
leveneTest(Gene ~ Sample*Species, data=testdata)




##############################################################################################################################
# Boxplots for genes with typical profiles
##############################################################################################################################
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


#load data and define samples
res <- readRDS(here("Data_Processed", "GENENAME_norm.Rds"))
res <- cbind(res[grep('HumanCell', colnames(res))],
             res[grep('MouseC2C12', colnames(res))],
             res[grep('RatL6', colnames(res))],
             res[grep('HumanTissue', colnames(res))],
             res[grep('MouseTissue', colnames(res))],
             res[grep('RatTissue', colnames(res))])
Sample <- c(rep("HSMC",  length(grep('HumanCell', colnames(res)))),
            rep("C2C12",   length(grep('MouseC2C12',  colnames(res)))),
            rep("L6",        length(grep('RatL6',       colnames(res)))),
            rep("Human Tissue", length(grep('HumanTissue',   colnames(res)))),  #list of sample types
            rep("Mouse Tissue",   length(grep('MouseTissue',  colnames(res)))),
            rep("Rat Tissue",   length(grep('RatTissue',  colnames(res)))))


#======================================================================================================
# Function to make boxplots
#======================================================================================================
PlotFunction <- function(genename) {
  data   <- data.frame()                                #create an empty dataframe to collect data
  for( i in 1:length(genename)) { 
    y     <- as.numeric(res[genename[i],])              #collect data for gene name i
    datay <- cbind.data.frame(Sample, y, rep(genename[i]))   #create table with x="sample type", y="data", "gene name"
    colnames(datay) <- c("x","y","Gene")                #rename column names to make it possible to rbind later
    data  <- rbind.data.frame(data, datay)              #bind the new gene data at the bottom of the previous one
  }
  data$x <- factor(data$x, levels=c("HSMC", "C2C12",  "L6", 
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


#======================================================================================================
# Boxplots of typical gene profiles (based on the file "Stats_ANOVA")
#======================================================================================================
# Representative genes selected in Stats for the effect of species, model or interation.
stats <- read.delim(here("Stats", "2-way_ANOVA.txt"))[,22:24]


#typical species effect
species <- stats[stats$Species.FDR<0.01 & stats$Sample.FDR>0.5 & stats$Interaction.FDR>0.5,]
species <- species[order(species$Species.FDR),]
t1 <- PlotFunction(rownames(species[1,]))
t2 <- PlotFunction(rownames(species[2,]))
t1

#typical cell vs tissue effect
sample <- stats[stats$Species.FDR>0.5 & stats$Sample.FDR<0.01 & stats$Interaction.FDR>0.5,]
sample <- sample[order(sample$Sample.FDR),]
t3 <- PlotFunction(rownames(sample[1,]))
t4 <- PlotFunction(rownames(sample[1,]))
t3

#typical interaction
interaction <- stats[stats$Species.FDR<.001 & stats$Sample.FDR<0.001 & stats$Interaction.FDR<0.001,]
interaction <- interaction[order(interaction$Interaction.FDR),]
t5 <- PlotFunction(rownames(interaction[1,]))
t6 <- PlotFunction(rownames(interaction[1,]))
t5

# Print figure
tiff(filename=here("Figures", "Profiles.tiff"), #print graph
    units="cm", width=17.2, height=12, 
    pointsize=12, res=1200)
matrix <- rbind(c(1,2,3), c(4,5,6))
grid.arrange(t1, t3, t5,
             t2, t4, t6,
             layout_matrix=matrix)
dev.off()


#======================================================================================================
# Gene markers of glycolysis, contraction and beta-oxidation
#======================================================================================================
# Glycolysis: phosphofructokinase (PFKM) step is the rate-limiting step
# Beta-oxidation: Carnitine Palmitoyltransferase (CPT1A) is the rate limiting step
# Contraction response: myosin heavy chain (MYH1) is specific to adult striated muscle
# Insulin response: glucose transporter system (GLUT4) is rate limiting
library(grid)
library(gridExtra)
tiff(filename=here("Figures", "Markers.tiff"), #print graph
    units="cm", width=18, height=12, 
    pointsize=12, res=300)
matrix <- rbind(c(1,2,3,4), c(5,6,7,8))
grid.arrange(PlotFunction('PFKM') + labs(title="Glycolysis"), 
             PlotFunction('CPT1A') + labs(title="Beta-oxidation"), 
             PlotFunction('MYH1') + labs(title="Contractile apparatus"),
             PlotFunction('SLC2A4') + labs(title="Glucose transport"),
             PlotFunction('LDHB') + labs(title="Glycolysis"), 
             PlotFunction('PPARG') + labs(title="Beta-oxidation"), 
             PlotFunction('RYR1') + labs(title="Contractile apparatus"),
             PlotFunction('SLC2A1') + labs(title="Glucose transport"), 
             layout_matrix=matrix)
dev.off()


#======================================================================================================
# Gene markers of differentiation
#======================================================================================================
library(grid)
library(gridExtra)
tiff(filename=here("Figures", "Markers_Differentiation.tiff"), #print graph
    units="cm", width=18, height=12, 
    pointsize=12, res=300)
matrix <- rbind(c(1,2,3,4), c(5,6,7,8))
grid.arrange(PlotFunction('DES'),
             PlotFunction('MYOG'),
             PlotFunction('DMD'),
             PlotFunction('CKM'),
             PlotFunction('MYOD1'),
             PlotFunction('MSTN'),
             PlotFunction('MYF5'),
             PlotFunction('PAX3'),
             layout_matrix=matrix)
dev.off()
