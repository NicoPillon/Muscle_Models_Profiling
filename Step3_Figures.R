##############################################################################################################################
#      Expression data in rat, mouse and human muscle, L6, C2C12 and human myotubes                                   ========
##############################################################################################################################
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
setwd("C:/ownCloud/Projects/MuscleModels/Data/Transcriptomics")

#load data and define samples
res <- readRDS("GENENAME_norm.Rds")
res <- cbind(res[grep('HumanCell', colnames(res))],
             res[grep('MouseC2C12', colnames(res))],
             res[grep('RatL6', colnames(res))],
             res[grep('HumanTissue', colnames(res))],
             res[grep('MouseTissue', colnames(res))],
             res[grep('RatTissue', colnames(res))])
Sample <- c(rep("Human Myotube",  length(grep('HumanCell', colnames(res)))),
            rep("Mouse C2C12",   length(grep('MouseC2C12',  colnames(res)))),
            rep("Rat L6",        length(grep('RatL6',       colnames(res)))),
            rep("Human Tissue", length(grep('HumanTissue',   colnames(res)))),  #list of sample types
            rep("Mouse Tissue",   length(grep('MouseTissue',  colnames(res)))),
            rep("Rat Tissue",   length(grep('RatTissue',  colnames(res)))))



#======================================================================================================
# Figure: Primary Component Analysis
#======================================================================================================
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



#======================================================================================================
# Figure: correlation contraction/proliferation
#======================================================================================================
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


#======================================================================================================
# Print Figure: PCA + correlation contraction/proliferation in one figure
#======================================================================================================
library(grid)
library(gridExtra)
png(filename="../../Figures/PCA+Corr.png", #print graph
    units="cm", width=18, height=8, 
    pointsize=12, res=300)
matrix <- rbind(c(1,1,1,1,1,1,2,2,2,2,2,2,2,2))
grid.arrange(pca, correlation,
             layout_matrix=matrix)
dev.off()



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


#======================================================================================================
# Boxplots of typical gene profiles (based on Stats)
#======================================================================================================
# Representative genes selected in Stats for the effect of species, model or interation.
stats <- read.delim("Stats/2-way_ANOVA.txt")[,22:24]




#typical species effect
species <- stats[stats$Species.FDR<0.01 & stats$Sample.FDR>0.5 & stats$Interaction.FDR>0.5,]
species <- species[order(species$Species.FDR),]
t1 <- PlotFunction("GRIN2B")
t2 <- PlotFunction("IFNA14")
t1

#typical cell vs tissue effect
sample <- stats[stats$Species.FDR>0.5 & stats$Sample.FDR<0.01 & stats$Interaction.FDR>0.5,]
sample <- sample[order(sample$Sample.FDR),]
t3 <- PlotFunction("LDHD")
t4 <- PlotFunction("SLC38A3")

#typical interaction
interaction <- stats[stats$Species.FDR<.001 & stats$Sample.FDR<0.001 & stats$Interaction.FDR<0.001,]
interaction <- interaction[order(interaction$Interaction.FDR),]
t5 <- PlotFunction("PLN")
t6 <- PlotFunction("CADM1")


# Print figure
png(filename="../../Figures/Profiles.png", #print graph
    units="cm", width=17.2, height=12, 
    pointsize=12, res=300)
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
png(filename="../../Figures/Markers.png", #print graph
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