library(here)
muscle  <- readRDS(here("Data_Processed", "GENENAME_norm.Rds"))
muscle <- cbind(muscle[grepl('HumanCell', colnames(muscle))],
                muscle[grepl('MouseC2C12', colnames(muscle))],
                muscle[grepl('RatL6', colnames(muscle))],
                muscle[grepl('HumanTissue', colnames(muscle))],
                muscle[grepl('MouseTissue', colnames(muscle))],
                muscle[grepl('RatTissue', colnames(muscle))])
colnames(muscle)
min(muscle, na.rm=T)
max(muscle, na.rm=T)


#list of gene names as RDS
list_genes <- rownames(muscle)
saveRDS(list_genes, here("Shiny_app", "data", "Muscle_Models_Profiling_genelist.Rds"))


#Make a list of data for ggplot
res <- muscle
x      <- c(rep('A1', length(grep('HumanCell',         colnames(muscle)))),  #list of sample types
            rep('A2', length(grep('MouseC2C12',   colnames(muscle)))),
            rep('A3', length(grep('RatL6',    colnames(muscle)))),
            rep('A4', length(grep('HumanTissue',    colnames(muscle)))),
            rep('A5', length(grep('MouseTissue', colnames(muscle)))),
            rep('A6', length(grep('RatTissue', colnames(muscle)))))
datalist <- vector("list", nrow(res))
names(datalist) <- rownames(res)
for (i in 1:nrow(res)){
  data   <- data.frame(x=factor(), y=numeric(), Gene=character(), stringsAsFactors=FALSE) #empty dataframe to collect data
  y     <- as.numeric(res[i,])            #collect data for gene name i
  data <- data.frame(x, y, rep(rownames(res[i,]))) #create table with x="sample type", y="data", "gene name"
  colnames(data) <- c("x","y","Gene")              #rename column names to make it possible to rbind later
  datalist[[i]] <- data
}

datalist[['SLC2A4']] #check example
saveRDS(datalist, here("Shiny_app", "data", "Muscle_Models_Profiling_data.Rds"))


#stats for all genes
res <- muscle
library(matrixStats)
statslist <- vector("list", nrow(res))
names(statslist) <- rownames(res)

for (i in 1:nrow(res)){
mean <- cbind(
  rowMeans(res[i, grepl('HumanCell',   colnames(res))], na.rm=T),
  rowMeans(res[i, grepl('MouseC2C12',        colnames(res))], na.rm=T),
  rowMeans(res[i, grepl('RatL6',  colnames(res))], na.rm=T),
  rowMeans(res[i, grepl('HumanTissue',colnames(res))], na.rm=T),
  rowMeans(res[i, grepl('MouseTissue',   colnames(res))], na.rm=T),
  rowMeans(res[i, grepl('RatTissue',colnames(res))], na.rm=T))
Sd <- cbind(
  rowSds(as.matrix(res[i, grepl('HumanCell',   colnames(res))]), na.rm=T),
  rowSds(as.matrix(res[i, grepl('MouseC2C12',        colnames(res))]), na.rm=T),
  rowSds(as.matrix(res[i, grepl('RatL6',  colnames(res))]), na.rm=T),
  rowSds(as.matrix(res[i, grepl('HumanTissue',colnames(res))]), na.rm=T),
  rowSds(as.matrix(res[i, grepl('MouseTissue',   colnames(res))]), na.rm=T),
  rowSds(as.matrix(res[i, grepl('RatTissue',colnames(res))]), na.rm=T))
nsize <- cbind(
  rowSums(!is.na(res[i, grepl('HumanCell',   colnames(res))])),
  rowSums(!is.na(res[i, grepl('MouseC2C12',        colnames(res))])),
  rowSums(!is.na(res[i, grepl('RatL6',  colnames(res))])),
  rowSums(!is.na(res[i, grepl('HumanTissue',colnames(res))])),
  rowSums(!is.na(res[i, grepl('MouseTissue',   colnames(res))])),
  rowSums(!is.na(res[i, grepl('RatTissue',colnames(res))])))
stats <- data.frame(t(mean), t(Sd), t(nsize))
stats[,3] <- as.factor(stats[,3])
colnames(stats) <- c('Mean', 'Sd', 'n')
rownames(stats) <- c("HumanCell", "MouseC2C12", "RatL6", 
                     "HumanTissue", "MouseTissue", "RatTissue") 
statslist[[i]] <- stats
}

statslist[['SLC2A4']] #check example
saveRDS(statslist, here("Shiny_app", "data", "Muscle_Models_Profiling_statslist.Rds"))



###################################################################################
# Deployment on shinyapps.io
#install.packages('rsconnect')
library(here)
library(rsconnect)
rsconnect::setAccountInfo(name='nicopillon',
                          token='359847F440F90E4AFB8A08FAD87DF504',
                          secret='9Tik3y8B9W26RJzkfLxqPbHGQzeltWye3I14pqD4')

rsconnect::deployApp(here("Shiny_app"), appName="Muscle_Models_Profiling")

