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
library(here)
cbPalette <- c("#56B4E9", "#D3C839", "#CC79A7", "#0072B2", "#E69F00", "#D55E00") #color palette for colorblind people
#             light blue    yellow     pink     dark blue    orange      red
cbShapes  <- c(   21    ,    21    ,    24    ,    24    ,    22    ,    22    ,    23    ,     23   )
cbLines   <- c(   'a'   ,    'b'   ,    'c'   ,    'd'   ,    'e'   ,    'f'   ,    'g'   ,     'h'  )


#load data and define samples
res <-  readRDS(here("Data_Processed", "GENENAME_norm.Rds"))
res <- cbind(res[grep('HumanCell', colnames(res))],
             res[grep('MouseC2C12', colnames(res))],
             res[grep('RatL6', colnames(res))])
Sample <- c(rep("HSMC",  length(grep('HumanCell', colnames(res)))),
            rep("C2C12", length(grep('MouseC2C12',  colnames(res)))),
            rep("L6",    length(grep('RatL6',       colnames(res)))))


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
  data$x <- factor(data$x, levels=c("HSMC", "C2C12",  "L6")) #for a box plot, x should be a factor
  ggplot(data, aes(x=x, y=y, fill=x)) +  
    geom_boxplot() +
    labs(x="",
         y=paste(genename)) + 
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

PlotFunction('RAC1')
PlotFunction('CAMK2A')
PlotFunction('MYOC')


#======================================================================================================
# Gene markers of glycolysis, contraction and beta-oxidation
#======================================================================================================
# Glycolysis: phosphofructokinase (PFKM) step is the rate-limiting step
# Beta-oxidation: Carnitine Palmitoyltransferase (CPT1A) is the rate limiting step
# Contraction response: myosin heavy chain (MYH1) is specific to adult striated muscle
# Insulin response: glucose transporter system (GLUT4) is rate limiting
library(grid)
library(gridExtra)



tiff(filename=here("Figures", "Contraction.tiff"), #print graph
    units="cm", width=12, height=10, 
    pointsize=12, res=1200)
matrix <- rbind(c(1,2,3), c(4,5,6))
grid.arrange(PlotFunction('MYH1') + ylim(-4, 7.5),
             PlotFunction('MYH3') + ylim(-4, 7.5),
             PlotFunction('MYH4') + ylim(-4, 7.5),
             PlotFunction('MYH6') + ylim(-4, 7.5),
             PlotFunction('MYH7')+ ylim(-4, 7.5),
             PlotFunction('MYH9')+ ylim(-4, 7.5),
             layout_matrix=matrix)
dev.off()


tiff(filename=here("Figures", "Glc_Uptake_basal.tiff"), #print graph
    units="cm", width=(3*1.6), height=4, 
    pointsize=12, res=1200)
matrix <- rbind(c(1,2))
grid.arrange(PlotFunction('SLC2A1'),
             PlotFunction('SLC2A3'),
             layout_matrix=matrix)
dev.off()

tiff(filename=here("Figures", "Glc_Uptake_insulin.tiff"), #print graph
    units="cm", width=(3*1.6), height=4, 
    pointsize=12, res=1200)
matrix <- rbind(c(1,2))
grid.arrange(PlotFunction('SLC2A4'),
             PlotFunction('PIK3CD'),
             layout_matrix=matrix)
dev.off()

tiff(filename=here("Figures", "Glyg_Synthesis_basal.tiff"), #print graph
    units="cm", width=(3*1.6), height=4, 
    pointsize=12, res=1200)
matrix <- rbind(c(1,2))
grid.arrange(PlotFunction('GYS1'),
             PlotFunction('GYS2'),
             layout_matrix=matrix)
dev.off()

tiff(filename=here("Figures", "Glyg_Synthesis_insulin.tiff"), #print graph
    units="cm", width=(3*1.6), height=4, 
    pointsize=12, res=1200)
matrix <- rbind(c(1,2))
grid.arrange(PlotFunction('GSK3A'),
             PlotFunction('GSK3B'),
             layout_matrix=matrix)
dev.off()

tiff(filename=here("Figures", "FAOx.tiff"), #print graph
    units="cm", width=6.5, height=4.3, 
    pointsize=12, res=1200)
matrix <- rbind(c(1,2))
grid.arrange(PlotFunction('CPT1B'),
             PlotFunction('CPT2'),
             layout_matrix=matrix)
dev.off()

tiff(filename=here("Figures", "GlcOx.tiff"), #print graph
    units="cm", width=6.5, height=4.3, 
    pointsize=12, res=1200)
matrix <- rbind(c(1,2))
grid.arrange(PlotFunction('PDHA1'),
             PlotFunction('PDHB'),
             layout_matrix=matrix)
dev.off()

tiff(filename=here("Figures", "Glycolysis.tiff"), #print graph
    units="cm", width=(3*2), height=6, 
    pointsize=12, res=1200)
matrix <- rbind(c(1,2))
grid.arrange(PlotFunction('PFKM'),
             PlotFunction('LDHA'),
             layout_matrix=matrix)
dev.off()


tiff(filename=here("Figures", "Proliferation.tiff"), #print graph
    units="cm", width=9.5, height=5, 
    pointsize=12, res=1200)
matrix <- rbind(c(1,2,3))
grid.arrange(PlotFunction('MKI67'),
             PlotFunction('PLK1'),
             PlotFunction('CCNB1'),
             layout_matrix=matrix)
dev.off()


tiff(filename=here("Figures", "Lactate.tiff"), #print graph
    units="cm", width=6, height=4.5, 
    pointsize=12, res=1200)
matrix <- rbind(c(1,2))
grid.arrange(PlotFunction('LDHA') + ylim(-1.5,7),
             PlotFunction('LDHB') + ylim(-1.5,7),
             layout_matrix=matrix)
dev.off()


tiff(filename=here("Figures", "CS.tiff"), #print graph
    units="cm", width=5, height=4.5, 
    pointsize=12, res=1200)
matrix <- rbind(c(1,2))
grid.arrange(PlotFunction('CS') + ylim(-0.5,5),
             PlotFunction('CS') + ylim(-0.5,5),
             layout_matrix=matrix)
dev.off()


tiff(filename=here("Figures", "EPS_GlcUptake.tiff"), #print graph
    units="cm", width=7.5, height=4.25, 
    pointsize=12, res=1200)
matrix <- rbind(c(1,2,3))
grid.arrange(PlotFunction('CAMK2A'),
             PlotFunction('TBC1D1'),
             PlotFunction('RAC1'),
             layout_matrix=matrix)
dev.off()

