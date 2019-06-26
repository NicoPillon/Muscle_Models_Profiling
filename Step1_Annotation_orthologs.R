#===============================================================================================
# Expression data in rat, mouse and human muscle, L6, C2C12 and human myotubes
#===============================================================================================
library(here)
# This file requires that:
# - all the CEL files are downloaded and placed in placed in /Data/Transcriptomics/
# - all the CEL files are annotated by sample type (HumanCell, MouseTissue, etc)
# - all the CEL files are sorted by platform in specific subfolders (GPL81, GPL6244, etc)
# - all corresponding libraries are installed:
BiocManager::install('pd.mg.u74av2', ask=F) #Affymetrix Murine Genome U74A Version 2 Array
BiocManager::install('pd.hg.u133.plus.2', ask=F) #Affymetrix Human Genome U133 Plus 2.0 Array
BiocManager::install('pd.mouse430.2', ask=F) #Affymetrix Mouse Expression 430 2.0 Array
BiocManager::install('pd.rat230.2', ask=F) #Affymetrix Rat Genome 230 2.0 Array
BiocManager::install('pd.hugene.1.0.st.v1', ask=F) #Affymetrix Human Gene 1.0 ST Array
BiocManager::install('pd.mogene.2.1.st', ask=F) #Affymetrix Mouse Gene 2.1 ST Array
BiocManager::install('pd.hta.2.0', ask=F) #Affymetrix Human Transcriptome Array 2.0
# - all necessary libraries are installed:
BiocManager::install("oligo")
BiocManager::install("oligoClasses")
BiocManager::install("Rcpp")
BiocManager::install('biomaRt', ask=F)


#===============================================================================================
# Make annotation file with BiomaRt with human, rat and mouse
# On each plateform, probes are annotated differently so we need a master file
# with the proper ENSEMBL annotation for all the arrays used.
#===============================================================================================
require(biomaRt)
ensembl <- useMart("ensembl")
mart <- listMarts(ensembl)

# rat annotation
ensembl <- useDataset("rnorvegicus_gene_ensembl",mart=ensembl)
list <- listAttributes(ensembl)
annot1  <- getBM(attributes = c("affy_rat230_2", "ensembl_gene_id","external_gene_name"), mart = ensembl)
annot2  <- getBM(attributes = c("agilent_sureprint_g3_ge_8x60k", "ensembl_gene_id","external_gene_name"), mart = ensembl)
annot4  <- getBM(attributes = c("agilent_wholegenome_4x44k_v1", "ensembl_gene_id","external_gene_name"), mart = ensembl)
annot5  <- getBM(attributes = c("agilent_wholegenome_4x44k_v3", "ensembl_gene_id","external_gene_name"), mart = ensembl)

# human annotation
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
list <- listAttributes(ensembl)
annot6 <- getBM(attributes = c("affy_hg_u133_plus_2", "ensembl_gene_id","external_gene_name"), mart = ensembl)
annot7 <- getBM(attributes = c("affy_hugene_1_0_st_v1", "ensembl_gene_id","external_gene_name"), mart = ensembl)
annot8 <- getBM(attributes = c("affy_hta_2_0", "ensembl_gene_id","external_gene_name"), mart = ensembl)
annot9 <- annot8
annot9$affy_hta_2_0 <- paste(annot8$affy,".1", sep="")

# mouse annotation
ensembl <- useDataset("mmusculus_gene_ensembl",mart=ensembl)
list <- listAttributes(ensembl)
annot10 <- getBM(attributes = c("affy_mg_u74av2", "ensembl_gene_id","external_gene_name"), mart = ensembl)
annot11 <- getBM(attributes = c("affy_mogene_2_1_st_v1", "ensembl_gene_id","external_gene_name"), mart = ensembl)
annot12 <- getBM(attributes = c("affy_mouse430_2", "ensembl_gene_id","external_gene_name"), mart = ensembl)

# merge all annotation
colnames(annot1) <- c('affy', 'ENSEMBL', 'GENENAME')
colnames(annot2) <- c('affy', 'ENSEMBL', 'GENENAME')
colnames(annot4) <- c('affy', 'ENSEMBL', 'GENENAME')
colnames(annot5) <- c('affy', 'ENSEMBL', 'GENENAME')
colnames(annot6) <- c('affy', 'ENSEMBL', 'GENENAME')
colnames(annot7) <- c('affy', 'ENSEMBL', 'GENENAME')
colnames(annot8) <- c('affy', 'ENSEMBL', 'GENENAME')
colnames(annot9) <- c('affy', 'ENSEMBL', 'GENENAME')
colnames(annot10) <- c('affy', 'ENSEMBL', 'GENENAME')
colnames(annot11) <- c('affy', 'ENSEMBL', 'GENENAME')
colnames(annot12) <- c('affy', 'ENSEMBL', 'GENENAME')
annot <- rbind(annot1, annot2, annot4, annot5, annot6, 
               annot7, annot8, annot9, annot10, annot11, annot12)
annot[annot == "" ] = NA
annot <- annot[!is.na(annot$affy),]
annot$GENENAME <- NULL
saveRDS(annot, file=here("Data", "Transcriptomics", "Annotation_ENSEMBL.Rds"))


#=============================================================================================================================
# Batch normalize and annotate each platform with ENSEMBL
#=============================================================================================================================
# This section will look for CEL files, apply RMA normalization and annotate with ENSEMBL.
# The final file for each platform is in each folder in the format "GPLxxx_ENSEMBL.Rds".
annot <- readRDS(here("Data_Processed", "Annotation_ENSEMBL.Rds"))
GEOs <- c('GPL81',   #Affymetrix Murine Genome U74A Version 2 Array
          'GPL570',  #Affymetrix Human Genome U133 Plus 2.0 Array
          'GPL1261', #Affymetrix Mouse Genome 430 2.0 Array
          'GPL1355', #Affymetrix Rat Genome 230 2.0 Array
          'GPL6244', #Affymetrix Human Gene 1.0 ST Array
          'GPL17400',#Affymetrix Mouse Gene 2.1 ST Array
          'GPL17586' #Affymetrix Human Transcriptome Array 2.0
) 

for (i in 1:length(GEOs)) {
  library(oligo)
  celFiles <- list.celfiles(here("Data_Raw", "Transcriptomics", GEOs[i]), full=T, listGzipped=T)
  data <- read.celfiles(celFiles, checkType=F)
  eset <- rma(data)
  eset  <- data.frame(exprs(eset))
  annotdata <- merge(annot, eset, by.x=1, by.y=0, all=F)
  annotdata <- data.frame(aggregate(annotdata[,3:ncol(annotdata)],by = list(annotdata$ENSEMBL),FUN = mean), row.names = 1)
  colnames(annotdata) <- paste(GEOs[i], "_", colnames(annotdata), sep="")
  saveRDS(annotdata, file=here("Data_Processed", paste(GEOs[i], "_ENSEMBL.Rds", sep=""))) # export file
}


#===============================================================================================
# Find ortholog genes with BiomaRt 
# Because of the different annotation between mouse, human and rat, we need to use the 
# corresponding orthologs for each species.
#===============================================================================================
require(biomaRt)

# orthologs human
ensMart<-useMart("ensembl")
human = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart <- listMarts(human)
list <- listAttributes(human)
attributes = c("ensembl_gene_id",
               "rnorvegicus_homolog_ensembl_gene",
               "mmusculus_homolog_ensembl_gene",
               "external_gene_name")
human_orthologs = getBM(attributes, values=TRUE, mart=human, uniqueRows=T)
colnames(human_orthologs) <- c('human', 'rat', 'mouse', 'genename')
human_orthologs[human_orthologs==''] <- NA

# orthologs mouse
ensMart<-useMart("ensembl")
mouse = useMart("ensembl", dataset="mmusculus_gene_ensembl")
mart <- listMarts(mouse)
list <- listAttributes(mouse)
attributes = c("ensembl_gene_id",
               "rnorvegicus_homolog_ensembl_gene",
               "hsapiens_homolog_ensembl_gene",
               "external_gene_name")
mouse_orthologs = getBM(attributes, values=TRUE, mart=mouse, uniqueRows=T)
mouse_orthologs$external_gene_name <- toupper(mouse_orthologs$external_gene_name)
colnames(mouse_orthologs) <- c('mouse', 'rat', 'human', 'genename')
mouse_orthologs[mouse_orthologs==''] <- NA

# orthologs rat
ensMart<-useMart("ensembl")
rat = useMart("ensembl", dataset="rnorvegicus_gene_ensembl")
mart <- listMarts(rat)
list <- listAttributes(rat)
attributes = c("ensembl_gene_id",
               "hsapiens_homolog_ensembl_gene",
               "mmusculus_homolog_ensembl_gene",
               "external_gene_name")
rat_orthologs = getBM(attributes, values=TRUE, mart=rat, uniqueRows=T)
rat_orthologs$external_gene_name <- toupper(rat_orthologs$external_gene_name)
colnames(rat_orthologs) <- c('rat', 'human', 'mouse', 'genename')
rat_orthologs[rat_orthologs==''] <- NA

# merge by ENSEMBL for the genes for which ENSEMBL codes of orthologs are known
library(dplyr)
ENSEMBL <- full_join(human_orthologs, mouse_orthologs)
ENSEMBL <- full_join(ENSEMBL, rat_orthologs)
ENSEMBL[ENSEMBL==''] <- NA
ENSEMBL <- ENSEMBL[!is.na(ENSEMBL$human),] # not in human
ENSEMBL <- ENSEMBL[rowSums(is.na(ENSEMBL)) < 2, ] #remove genes with only one species present

# Merge by GENENAME for the genes for which ENSEMBL annotation is missing
# but have the same gene symbol in all three species
GENENAME <- merge(human_orthologs[,c(1,4)], mouse_orthologs[,c(1,4)], by=2, all=F)
GENENAME <- unique(GENENAME)
GENENAME <- merge(GENENAME, rat_orthologs[,c(1,4)], by.x=1, by.y=2, all=F)
GENENAME <- unique(GENENAME)

# full merge
All <- full_join(ENSEMBL, GENENAME)

# Clean from pseudogenes and unwanted genes
All <- All
All <- All[!is.na(All$human),] # not in human
All <- All[!grepl("MIR", All$genename),] # microRNA
All <- All[!grepl("\\bLOC", All$genename),] # unknown locus
All <- All[!grepl("orf", All$genename),] # orfans
All <- All[!grepl("RPL", All$genename),] # ribosomal proteins
All <- All[!grepl("RPS", All$genename),] # ribosomal proteins
All <- All[!grepl("RNU", All$genename),] # ribosomal proteins
All <- All[!grepl("SNOR", All$genename),] # snoRNA
All <- All[!grepl("RF0", All$genename),] # undefined RNA
All <- All[!grepl("RIK", All$genename),] # Riken annotation

saveRDS(All, file=here("Data_Processed", "Annotation_Orthologs.Rds"))


#=============================================================================================================================
# Merge data from all platforms and normalize
# This section will merge the data from all plateforms, annotate with orthologs and normalize. 
#=============================================================================================================================
# create a list of all ENSEMBL files
list.filenames <- list.files(here("Data_Processed"), pattern="ENSEMBL.Rds", recursive=T)
list.data <- list()
for (i in 1:length(list.filenames)) {
  list.data[[i]]<-readRDS(here("Data_Processed", list.filenames[i]))
}
names(list.data) <- list.filenames

# merge all by ENSEMBL name
library(dplyr)
All  <- data.frame(rownames(list.data[[1]]))
colnames(All) <- 'rownames'
for ( df in list.data ) {
  asd <- df
  asd$rownames <- rownames(asd)
  All <- full_join(All, asd)
}
All <- data.frame(All, row.names = 1)
rm(list.data, asd, df, i, list.filenames)

# merge by orthologs human using custom made file (see Across_Species_Orthologs.R)
orthologs <- readRDS(here("Data_Processed", "Annotation_Orthologs.Rds"))
human <- merge(orthologs, All[,grep('Human', colnames(All))], by.x=1, by.y=0, all=T)
rm(orthologs) #needed to free some memory on slow machines
rat   <- merge(human,     All[,grep('Rat', colnames(All))], by.x=2, by.y=0, all=T)
rm(human) #needed to free some memory on slow machines
mouse <- merge(rat,       All[,grep('Mouse', colnames(All))], by.x=3, by.y=0, all=T)
rm(rat, All) #needed to free some memory on slow machines
Final <- mouse[!is.na(mouse$genename),]
rm(mouse) #needed to free some memory on slow machines

# Aggregate by gene symbol
Final <- data.frame(aggregate(Final[,5:ncol(Final)], by=list(Final$genename), FUN=mean, na.rm=TRUE), row.names = 1)
Final[Final == 'NaN'] <- NA

# delete genes that contain too many NAs
noNA <- Final[rowSums(is.na(Final)) / ncol(Final) < 20/100, ]

# Normalize between arrays
library(limma)
Norm <- normalizeBetweenArrays(noNA, method="quantile")
Norm <- data.frame(Norm - median(Norm, na.rm = T))
boxplot(Norm)
Norm <- Norm[order(colnames(Norm))]

saveRDS(Norm, here("Data_Processed", "GENENAME_norm.Rds"))


#=============================================================================================================================
# Batch correction (Not used in manuscript)
# I tried to to a batch correction using the different plateforms as batches. However, this erases differences
# between species because the gene arrays are different in human/rat/mouse. In the end, I decided to not use
# any batch correction. The downside of this choice is that I am unable to seprate differences dues to species 
# from differences due to the different arrays used.
#=============================================================================================================================
res <- readRDS(here("Data_Processed", "GENENAME_batch.Rds"))

library(stringr)
batch <- gsub("_.*","",colnames(res) ) # Extract the name of the array used from the column names
unique(batch) # there are 8 different types of arrays used in the analysis

library(limma)
res <- removeBatchEffect(res, batch)
res <- normalizeBetweenArrays(res, method="quantile")
res <- data.frame(res)
boxplot(res)

saveRDS(res, here("Data_Processed", "GENENAME_batch.Rds"))
