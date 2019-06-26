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
venn.diagram(x=list(A=rownames(Human), B=rownames(Mouse), C=rownames(Rat)), filename=here("Figures", "VennDiagram.png"),
             units="cm", width=9, height=5, res=300, margin=0.05,
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

# Extract the residuals and run Shapiro-Wilk test
aov_residuals <- residuals(object = res.aov)
shapiro.test(x = aov_residuals )


pairwise.t.test(testdata$Gene, testdata$Species,
                p.adjust.method = "BH")

library(multcomp)
compstats <- glht(res.aov, linfct = mcp(Species = "Tukey"))
summary(compstats)
