############################################################################################
#=== Expression data in rat muscle, L6 cells, human myotubes and human muscle ==============
############################################################################################
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(gplots)
library(stringr)
library(grid)
library(gridExtra)
library(here)
cbPalette <- c("#E69F00", "#0072B2", "#CC79A7", "#009E73", "#D3C839", "#BC5300", "#84C4E8", "#000000") #color palette for colorblind people
cbShapes  <- c(   21    ,    21    ,    24    ,    24    ,    22    ,    22    ,    23    )
cbLines   <- c(   'a'   ,    'b'   ,    'c'   ,    'd'   ,    'e'   ,    'f'   ,    'g'   )


rawdata <- readRDS(here("Data_Processed", "GENENAME_batch.Rds"))
rawdata <- rawdata[!grepl('HEK', colnames(rawdata))] #remove HEK cells
rawdata <- rawdata[!grepl('HeLa', colnames(rawdata))] #remove HeLa cells

All <- data.frame(HC=rowMeans(rawdata[grep('HumanCell', colnames(rawdata))], na.rm=T),
                  MC=rowMeans(rawdata[grep('MouseC2C12', colnames(rawdata))], na.rm=T),
                  RT=rowMeans(rawdata[grep('RatL6', colnames(rawdata))], na.rm=T))
All <- 2^All

#theme for figures
library(ggplot2)
theme <- theme(plot.title  = element_text(face="bold", color="black", size=7, angle=0),
               axis.text.x = element_text(color="black", size=6, angle=45, hjust = 1),
               axis.text.y = element_text(color="black", size=6, angle=0, hjust = 0.5),
               axis.title  = element_text(face="bold", color="black", size=7, angle=0),
               legend.text = element_text(color="black", size=4.5, angle=0, hjust = 0),
               legend.key.size = unit(0.3, "cm"))


#=====================================================================================================================
# Oxidative metabolism
Complex1 <- c('ND1',     'ND2',     'ND3',     'ND4',    'ND4L',   'ND5',     'ND6',    'NDUFS1', 'NDUFS2',
              'NDUFS3',  'NDUFS4',  'NDUFS5',  'NDUFS6', 'NDUFS7', 'NDUFS8',  'NDUFV1', 'NDUFV2', 'NDUFV3',
              'NDUFAB1', 'NDUFA1',  'NDUFA2',  'NDUFA3', 'NDUFA4', 'NDUFA5',  'NDUFA6', 'NDUFA7', 'NDUFA8',
              'NDUFA9',  'NDUFA10', 'NDUFA11', 'NDUFA12','NDUFA13','NDUFB1',  'NDUFB2', 'NDUFB3', 'NDUFB4',
              'NDUFB5',  'NDUFB6',  'NDUFB7',  'NDUFB8', 'NDUFB9', 'NDUFB10', 'NDUFB11','NDUFC1', 'NDUFC2')
Complex2 <- c('SDHA',    'SDHB',    'SDHC',    'SDHD')
Complex3 <- c('CYC1',    'CYTB',    'UQCRB',   'UQCRC1', 'UQCRC2', 'UQCRFS1', 'UQCRH',  'UQCRQ',  'UQCR10',
              'UQCR11',  'BCS1L' )
Complex4 <- c('COX1',    'COX2',    'COX3',    'COX4I1', 'COX4I2', 'COX5A',   'COX5B',  'COX6A1', 
              'COX6A2',  'COX6B1',  'COX6B2',  'COX6C',  'COX7A1', 'COX7A2',  'COX7A2L','COX7B',
              'COX7B2',  'COX7C',   'COX8A',   'COX8C')
Complex5 <- c('ATP5A1',  'ATP5B',   'ATP5C1',  'ATP5D',  'ATP5E',  'ATP5G1',  'ATP5G2', 'ATP5G3',
              'USMG5',   'ATP5I',   'ATP5J2',  'ATP5L',  'C14orf2','ATP6',    'ATP8',   'ATP5F1',
              'ATP5H',   'ATP5J',   'ATP5O',   'ATPIF1', 'OXA1L')
Complex1.mean <- colMeans(All[Complex1, ], na.rm=T)
Complex2.mean <- colMeans(All[Complex2, ], na.rm=T)
Complex3.mean <- colMeans(All[Complex3, ], na.rm=T)
Complex4.mean <- colMeans(All[Complex4, ], na.rm=T)
Complex5.mean <- colMeans(All[Complex5, ], na.rm=T)
Complex1.sd <- colMeans(All[Complex1, ], na.rm=T)
Complex2.sd <- colMeans(All[Complex2, ], na.rm=T)
Complex3.sd <- colMeans(All[Complex3, ], na.rm=T)
Complex4.sd <- colMeans(All[Complex4, ], na.rm=T)
Complex5.sd <- colMeans(All[Complex5, ], na.rm=T)
values <- as.numeric(c(Complex1.mean, Complex2.mean, Complex3.mean, Complex4.mean, Complex5.mean))
sd <- as.numeric(c(Complex1.sd, Complex2.sd, Complex3.sd, Complex4.sd, Complex5.sd))
Gene <-c(rep("Complex I", 3), rep("Complex II", 3), rep("Complex III", 3),
         rep("Complex IV", 3), rep("Complex V", 3))
Exercise <-factor(rep(c("HSMC", "C2C12", "L6"), length(unique(Gene))), levels=c('HSMC', 'C2C12', 'L6'))
mydata <-data.frame(Exercise, values, sd, Gene)

Respiration <- ggplot(mydata, aes(Exercise, values, fill=Gene)) +
  geom_bar(stat="identity", colour="black", size=0.2) +
  labs(x= "",
       y= "OXPHOS",
       fill="",
       label="")  + theme_bw() + theme +
  scale_fill_brewer(palette="Set3")
Respiration

png(filename=here("Figures", "Mitochondria.png"), #print graph
    units="cm", width=5, height=5, 
    pointsize=12, res=600)
Respiration
dev.off()
Respiration




#=====================================================================================================================
# Figure for Contractile Proteins 
#MYH2, MYH7B, MYH8, MYH15, MYH16 and MYL5 are not present in the dataset
Genelist <-c("ACTA1", "ACTA2", "ACTC1",
             "MYH1", "MYH3", "MYH4", "MYH6", "MYH7", "MYH9", "MYH10", "MYH11", "MYH13", "MYH14",
             "MYL1", "MYL2", "MYL3", "MYL4", "MYL6", "MYL6B", "MYL7", "MYL9", "MYL10", "MYL12A", "MYL12B", "MYLPF")
Values <- numeric()
for (i in 1:length(Genelist)){   Values <- c(Values, as.numeric(All[Genelist[i],])) }
Samples <-factor(rep(c("HSMC", "C2C12", "L6"), length(unique(Genelist))), levels=c('HSMC', 'C2C12', 'L6'))
Genes <- character()
for (i in 1:length(Genelist)){   Genes <- c(Genes, rep(Genelist[i], 3)) }
mydata <-data.frame(Samples, Values, Genes)
mydata$Genes <- factor(mydata$Genes, levels=Genelist)
Contraction <- ggplot(mydata, aes(Samples, Values, fill=Genes, cutoff = factor(0) )) + theme_bw() +
  geom_bar(stat="identity", colour="black", size=0.2) +
  labs(x= "", y= "Relative abundance",
       fill="") +  theme 
Contraction

png(filename=here("Figures", "Contraction.png"), #print graph
    units="cm", width=6, height=7, 
    pointsize=12, res=600)
Contraction
dev.off()



#=====================================================================================================================
# Figure for Insulin signalling
INSR <- paste(All['INSR',])
IRS1 <- paste(All['IRS1',])
IRS2 <- paste(All['IRS2',])
PI3K  <- paste((All['PIK3R1', ]+All['PIK3R2', ]+All['PIK3CA', ] +All['PIK3CB',]  +All['PIK3CD', ])/5)
PDK1 <- paste(All['PDK1',])
AKT1 <- paste(All['AKT1',])
AKT2 <- paste(All['AKT2',])
TBC1D4 <- paste(All['TBC1D4',])
GLUT4 <- paste(All['SLC2A4',])
Gene <-c(rep("INSR", 3), rep("IRS1", 3), rep("IRS2", 3),  rep("PI3K", 3), rep("PDK1", 3),
         rep("AKT1", 3), rep("AKT2", 3), rep("TBC1D4", 3), rep("GLUT4", 3))
values <-as.numeric(c(INSR, IRS1, IRS2, PI3K, PDK1, AKT1, AKT2, TBC1D4, GLUT4))
Exercise <-rep(c("Human Myotube",
                 "Mouse C2C12",
                 "Rat L6"), 9)
mydata <-data.frame(Exercise, values, Gene)
mydata$Gene <- factor(mydata$Gene, levels=c('INSR','IRS1','IRS2','PI3K','PDK1','AKT1','AKT2','TBC1D4','GLUT4')) #for a box plot, x should be a factor

InsulinSignalling <- ggplot(mydata, aes(Exercise, values, fill=Gene)) +
  geom_bar(stat="identity", colour="black", size=0.2) +

  labs(title="Insulin signalling",
       x= "", y= "Relative abundance",
       fill="") + theme +
  scale_fill_brewer(palette="Set3")
InsulinSignalling


#=====================================================================================================================
# Fatty Acid Metabolism
ACAT <- paste((All['ACAA1', ]+All['ACAA2', ]+All['ACAT1', ])/3)        #Acetyl-CoA Transferases
ACAD <- paste((All['ACAD9', ]+All['ACAD10',]                           #Acyl-CoA Dehydrogenases
             +All['ACAD11',]+All['ACADL', ]
             +All['ACADM', ]+All['ACADS', ]
             +All['ACADSB',]+All['ACADVL',]
             +All['EHHADH',]+All['GCDH',  ])/10)
ACOX  <- paste((All['ACOX1',]+All['ACOX2', ]+All['ACOX3',])/3) #Acyl-CoA Oxidases
ACS   <- paste((All['ACSBG1',]                          #Acyl-CoA Synthetases
              +All['ACSL1', ]+All['ACSL3',]
              +All['ACSL4', ]+All['ACSL5',]
              +All['ACSL6', ]+All['ACSM3', ])/10)
ACOT   <- paste((All['ACOT1',]+All['ACOT2', ]                          #Acyl-CoA Thioesterases
               +All['ACOT7',]
               +All['ACOT8', ]+All['ACOT9',]
               +All['ACOT12', ])/10)
CPT   <- paste((All['CPT1A',]+All['CPT1B', ]                          #Carnitine Transferases
              +All['CPT1C', ]+All['CPT2',]
              +All['CRAT', ]+All['CROT',])/9)
FABP  <- paste((All['FABP1',]+All['FABP2', ]                          #Fatty Acid Transport
              +All['FABP3', ]+All['FABP4',]
              +All['FABP5', ]+All['FABP6',])/9)
AMPK   <- paste((All['PRKAA1',]+All['PRKAA2', ]                        #Fatty Acid Biosynthesis Regulation (AMPK)
                +All['PRKAB1', ]+All['PRKAB2',]
                +All['PRKACA', ]+All['PRKACB',]
                +All['PRKAG2', ]+All['PRKAG3',])/9)

#ALDH2, DECR1, DECR2ECHS1, HADHA, MCEE, MUT, ECI2, PECR, PPA1 #Others
#SLC27A1, SLC27A2, SLC27A3, SLC27A4, SLC27A5, SLC27A6 #Fatty Acid Transport
#BDH1, BDH2, HMGCL, HMGCS1, HMGCS2, OXCT2 #Ketogenesis & Ketone Body Metabolism
#GK, GK2, GPD1, GPD2, LIPE, LPL #Triacylglycerol Metabolism

values <- as.numeric(c(ACAT, ACAD, ACOX, ACS, ACOT, CPT, FABP, AMPK))
Gene <-c(rep("ACAT", 3), rep("ACAD", 3), rep("ACOX", 3),  rep("ACS", 3),
         rep("ACOT", 3), rep("CPT", 3), rep("FABP", 3),  rep("AMPK", 3))
Exercise <-rep(c("Human Myotube",
                 "Mouse C2C12",
                 "Rat L6"), 8)
mydata <-data.frame(Exercise, values, Gene)
LipidMetab <- ggplot(mydata, aes(Exercise, values, fill=Gene)) +
  geom_bar(stat="identity", colour="black", size=0.2) +

  labs(title="Lipid Metabolism",
       x= "",  y= "Relative abundance",
       fill="") + theme +
  scale_fill_brewer(palette="Set3")
LipidMetab


#=====================================================================================================================
# Glycolysis
HK    <- paste(All['HK2',  ])                                   #Hexokinase
GPI   <- paste(All['GPI',  ])                                   #Phosphoglucose isomerase
PFKL  <- paste(All['PFKL', ])                                   #Phosphofructokinase
ALDO  <- paste((All['ALDOA',]+All['ALDOB',]+All['ALDOC',])/3)   #Aldolase
TPI   <- paste(All['TPI1', ])                                   #triosephosphate isomerase
GAPDH <- paste(All['GAPDH',])                                   #Glyceraldehyde-3-phosphate dehydrogenase
PGK   <- paste(All['PGK1', ])                                   #Phosphoglycerate kinase
PGAM  <- paste((All['PGM1', ]+All['PGM2', ]+All['PGM3', ]
                +All['PGAM1',]+All['PGAM2', ])/5)               #Phosphoglycerate mutase
ENO   <- paste((All['ENO1', ]+All['ENO2', ]+All['ENO3', ])/3)   #Phosphopyruvate hydratase (Enolase)
PKM   <- paste(All['PKM',  ])                                   #Pyruvate kinase muscle
LDH   <- paste((All['LDHA',  ]+All['LDHB',  ])/2)               #lactate dehydrogenase
values <- as.numeric(c(HK, GPI, PFKL, ALDO, TPI, GAPDH, PGK, PGAM, ENO, PKM, LDH))
Gene <-c(rep("HK", 3), rep("GPI", 3), rep("PFKL", 3), rep("ALDO", 3),
         rep("TPI", 3), rep("GAPDH", 3), rep("PGK", 3), rep("PGAM", 3),
         rep("ENO", 3), rep("PKM", 3), rep("LDH", 3))
Exercise <-factor(rep(c("HSMC", "C2C12", "L6"), length(unique(Gene))), levels=c('HSMC', 'C2C12', 'L6'))
mydata <-data.frame(Exercise, values, Gene)

Glycolysis <- ggplot(mydata, aes(Exercise, values, fill=Gene)) +
  geom_bar(stat="identity", colour="black", size=0.2) +
  labs(x= "", y= "Relative abundance",
       fill="") + theme_bw() + theme +
  scale_fill_brewer(palette="Set3")
Glycolysis


#=====================================================================================================================
# AMPK
PRKAA1 <- paste(All['PRKAA1',])
PRKAA2 <- paste(All['PRKAA2',])
PRKAB1 <- paste(All['PRKAB1',])
PRKAB2 <- paste(All['PRKAB2',])
PRKAG1 <- paste(All['PRKAG1',]) #not expressed
PRKAG2 <- paste(All['PRKAG2',])
PRKAG3 <- paste(All['PRKAG3',])
values <- as.numeric(c(PRKAA1, PRKAA2, PRKAB1, PRKAB2, PRKAG1, PRKAG2, PRKAG3))
Gene <-c(rep("PRKAA1", 3), rep("PRKAA2", 3),
         rep("PRKAB1", 3),  rep("PRKAB2", 3),
         rep("PRKAG1", 3),  rep("PRKAG2", 3),  rep("PRKAG3", 3))
Exercise <-factor(rep(c("HSMC", "C2C12", "L6"), length(unique(Gene))), levels=c('HSMC', 'C2C12', 'L6'))
mydata <-data.frame(Exercise, values, Gene)

AMPKisoforms <- ggplot(mydata, aes(Exercise, values, fill=Gene)) +
  geom_bar(stat="identity", colour="black", size=0.2) +
  
  labs(title="Lipid Metabolism",
       x= "",  y= "Relative abundance",
       fill="") + theme +
  scale_fill_brewer(palette="Set3")
AMPKisoforms
