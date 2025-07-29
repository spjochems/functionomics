rm(list=ls())

#packages that are needed
library(reshape2) 
library(ggplot2)
library(flowCore)
library(cowplot)
library(RColorBrewer)
library(plyr)
library(umap)
library(Rtsne)
library(stringi)
library(pheatmap)
library(flowCore)


####################################################################
#set the working directory
setwd('V:\\para-yazdanbakhsh\\JochemsLab\\Simon\\LUF aurora panel/2022_01_17 new comp and gating checked/final comp & files/exported populations/fcs/surface/')
dir.create('V:/para-yazdanbakhsh/JochemsLab/Simon/LUF aurora panel/2022_01_17 new comp and gating checked/analysis/Output FCS/Counts/total')

#read in data and fix column names
tubes <- list.files()
samples <- do.call(c, lapply(strsplit(tubes, split = '-'), function(i) i[2]))
samples <- samples[!duplicated(samples)]




for(tube in tubes){

setwd('V:\\para-yazdanbakhsh\\JochemsLab\\Simon\\LUF aurora panel/2022_01_17 new comp and gating checked/final comp & files/exported populations/fcs/surface/')
setwd(tube)

files <- list.files()

if(length(files) == 17){
  print(paste0('Correct number of surface pops for ', tube))
} else {
  print(paste0('Wrong number of surface pops for ', tube, '!!!!'))
}

#read in the first file
data <- read.FCS(files[1], transformation  = F, truncate_max_range =F)
dat2 <- data.frame(data@exprs)
desc <- data@parameters@data$desc
desc[is.na(desc)] <- data@parameters@data$name[is.na(data@parameters@data$desc)]
colnames(dat2) <- desc

#creat a column based on time and a marker that serves as a unique. 
#Normally there is an event parameter in these files, but not in this file so need to do this
len <- length(strsplit(files[1], '_')[[1]])
dat2$Population <- strsplit(files[1], '_')[[1]][len]

##########Add in the populations
###########################################################3
#run throught all the files
for(i in 2:length(files)){
  
  #get the names of the population and get rid of the extra stuff in it
  data <- read.FCS(files[i], transformation  = F, truncate_max_range =F)
  len <- length(strsplit(files[i], '_')[[1]])
  
  if(nrow(data)>0){
    dat <- data.frame(data@exprs)
    desc <- data@parameters@data$desc
    desc[is.na(desc)] <- data@parameters@data$name[is.na(data@parameters@data$desc)]
    colnames(dat) <- desc
    dat$Population <- strsplit(files[i], '_')[[1]][len]
  
  dat2 <- rbind(dat2, dat)
  }
}

#fix name
dat2$Population <- gsub('.fcs', '', dat2$Population)
table(dat2$Population) #make a table of all the populations

print(paste0('all pops are read:', length(table(dat2$Population)) == length(files))) #check that all the populations are read in, if this is FALSE,
#then perhaps some of the populations are overwritten by later ones that are higher in the gating layout 
#or different files with name in the folder
print(paste0('number of surface cells: ', nrow(dat2))) #check that all the barcodes are read in

#add ID
dat2$ID <- paste0(dat2$Time, '_', dat2$CD4, '_', dat2$CD3, '_', dat2$CD8,
                   '_', dat2$CD14, '_', dat2$CD19, '_', dat2$CD88, '_', dat2$CCL2)

print(paste0('number of duplicates surface: ', sum(duplicated(dat2$ID)))) # check duplicates

#remove duplicates
dat2 <- dat2[!duplicated(dat2$ID),]



##########Add in the barcodes
###########################################################3
#create an empty column that we can use to add the population info to
dat2$Barcode <- '0.NotAssigned'


setwd('V:\\para-yazdanbakhsh\\JochemsLab\\Simon\\LUF aurora panel/2022_01_17 new comp and gating checked/final comp & files/exported populations/fcs/barcodes/')
setwd(tube)
#load in the  fcs files for each of the populations and add the names to the cells in the large dataframe
files <- list.files()
print(paste0('number of stimulation files for ', tube, ' = ', length(files)))

#run through all the files
for(i in 1:length(files)){
  
  #get the names of the population and get rid of the extra stuff in it
  len <- length(strsplit(files[i], '_')[[1]])
  name <-  strsplit(files[i], '_')[[1]][len]
  name <-  gsub('.fcs', '', name)

  pop.data <- read.FCS(files[i], transformation  = F, truncate_max_range =F)
  pop.dat <- data.frame(pop.data@exprs)
  desc <- pop.data@parameters@data$desc
  desc[is.na(desc)] <- pop.data@parameters@data$name[is.na(pop.data@parameters@data$desc)]
  colnames(pop.dat) <- desc
  
  
  pop.dat$ID <- paste0(pop.dat$Time, '_', pop.dat$CD4, '_', pop.dat$CD3, '_', pop.dat$CD8,
                    '_', pop.dat$CD14, '_', pop.dat$CD19, '_', pop.dat$CD88, '_', pop.dat$CCL2)
  
    #add barcodes to cells
  dat2$Barcode[dat2$ID %in% pop.dat$ID] <- name
}


table(dat2$Barcode) #make a table of all the barcodes
print(paste0('number of loaded barcodes = ', length(table(dat2$Barcode))-1))

#change names
dat2$Barcode <- gsub('IL18', 'IL12_IL18', dat2$Barcode)
dat2$Barcode <- gsub('IL17-TSLP', 'IL33_IL25_TSLP', dat2$Barcode)
dat2$Barcode[dat2$Barcode == 'C'] <- 'PolyIC'
dat2$Barcode <- gsub('dC', 'dGdC', dat2$Barcode)


#######################################################################
#check the populations per barcode to make all populations are present equally
############################################################################
Combos <- data.frame(table(dat2$Population, dat2$Barcode)) #make a table
colnames(Combos) <- c('Population', 'Barcode', 'Count') #fix column names


Combos2 <- Combos[Combos$Barcode != '0.NotAssigned' ,] #remove the unassigned barcodes and populations

write.table(Combos2, paste0('V:/para-yazdanbakhsh/JochemsLab/Simon/LUF aurora panel/2022_01_17 new comp and gating checked/analysis/Output FCS/Counts/total/', tube, '.txt'), 
                                   sep = '\t', row.names = F,  quote = F)


}


