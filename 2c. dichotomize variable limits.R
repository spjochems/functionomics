rm(list = ls())

library(reshape2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
library(openCyto)
library(xlsx)

setwd('C:\\Users\\spjochems\\Desktop\\Current\\ICS\\Output FCS/Annotated')

#function
arcsinhn <- function(x, n){
  return(asinh(x/n))  
}



#get the cytokine names
markers <- read.delim('R:\\Para-CIH\\CIH Group member folders\\Simon\\Results\\flow\\LUF Aurora ICS\\2022_01_17 18 samples new\\panel.txt')
cytokines <- markers$Maker[markers$Type == 'Cytokine']

#get the cofactors
arcsin <- read.xlsx('R:\\Para-CIH\\CIH Group member folders\\Simon\\Results\\flow\\LUF Aurora ICS\\2022_01_17 18 samples new\\arcsin_cofactor.xlsx', 1)

#load in text files of mapped cells
files <- list.files(pattern = 'txt')

dir.create('../Cytokine gates')

for(i in 1:18){
  
  setwd('C:\\Users\\spjochems\\Desktop\\Current\\ICS\\Output FCS/Annotated')
  a <-  read.delim(files[i])
  b <-  read.delim(files[i+18])
  
  #get samples id from file name
  sampleid <- strsplit(files[i], '-')[[1]][2]
  sampleid <- gsub('.txt', '', sampleid) #cut out the .txt at end
  sampleid <- gsub('3173', '3713', sampleid) #fix misnamed file 3173
  
  
  setwd('../Cytokine gates')

  responses <- data.frame(Sample = NA,
                          Stimulation = NA, 
                          Population = NA,
                          Cytokine = NA,
                          Percentage = NA)
  
  #add duration
  a$Duration <- 'ON'
  b$Duration <- '4h'
  
  
  #fix names
  a$Population <- gsub('B cells', 'Bcells', a$Population)
  a$Population <- gsub('CD16 NK', 'CD16NK', a$Population)
  a$Population <- gsub('CD4[+]', 'CD4T', a$Population)
  a$Population <- gsub('CD56Hi NK', 'CD56HiNK', a$Population)
  a$Population <- gsub('CD56Mid NK', 'CD56MidNK', a$Population)
  a$Population <- gsub('CD8T conventional', 'CD8T', a$Population)
  a$Population[a$Population == 'Mono'] <- 'MonoIntClassic'
  a$Population <- gsub('DN T', 'DNT', a$Population)
  a$Population <- gsub('DP', 'DPT', a$Population)
  a$Population <- gsub('NKT cells', 'NKT', a$Population)
  a$Population <- gsub('Non-classical Mono', 'MonoNonClassic', a$Population)
  
  
  b$Population <- gsub('B cells', 'Bcells', b$Population)
  b$Population <- gsub('CD16 NK', 'CD16NK', b$Population)
  b$Population <- gsub('CD4[+]', 'CD4T', b$Population)
  b$Population <- gsub('CD56Hi NK', 'CD56HiNK', b$Population)
  b$Population <- gsub('CD56Mid NK', 'CD56MidNK', b$Population)
  b$Population <- gsub('CD8T conventional', 'CD8T', b$Population)
  b$Population[b$Population == 'Mono'] <- 'MonoIntClassic'
  b$Population <- gsub('DN T', 'DNT', b$Population)
  b$Population <- gsub('DP', 'DPT', b$Population)
  b$Population <- gsub('NKT cells', 'NKT', b$Population)
  b$Population <- gsub('Non-classical Mono', 'MonoNonClassic', b$Population)
  
  
  #check oke
  table(a$Population)
  table(b$Population)
  
  #add total monocytes and NK cells as extra groups
  a.monocytes <- a[grepl('Mono', a$Population),]
  a.monocytes$Population <- 'Monocytes'
  a <- rbind(a, a.monocytes)  
  a.nk <- a[grepl('CD16NK', a$Population) | 
              grepl('CD56', a$Population)  ,]
  a.nk$Population <- 'NK'
  a <- rbind(a, a.nk)  
  
  
  b.monocytes <- b[grepl('Mono', b$Population),]
  b.monocytes$Population <- 'Monocytes'
  b <- rbind(b, b.monocytes)  
  b.nk <- b[grepl('CD16NK', b$Population) | 
              grepl('CD56', b$Population)  ,]
  b.nk$Population <- 'NK'
  b <- rbind(b, b.nk)  
  
  #change name
  colnames(a)[colnames(a) == 'Barcode'] <- 'Stimulus'
  colnames(b)[colnames(b) == 'Barcode'] <- 'Stimulus'
  
  #set medium name
  a$Stimulus <- gsub('Medium', '0.ON_Medium', a$Stimulus) #set medium to be first on plot
  b$Stimulus <- gsub('Medium', '0.4h_Medium', b$Stimulus) #set medium to be first on plot
  
  #check oke
  table(a$Stimulus)
  table(b$Stimulus)
  
  
  #select cytokines, population and mulation, if put cytokines between '', the word cytokines is pasted in the dataframe, and it doesn't read cytokines files
  keep <- c('Population', 'Stimulus', 'Duration', cytokines)
  a1.keep <- a[,keep]
  b1.keep <- b[,keep]
  
  for(cyto in 4:ncol(a1.keep)){
    a1.keep[,cyto] <- arcsinhn(a1.keep[,cyto], arcsin[cyto-3,2])
    b1.keep[,cyto] <- arcsinhn(b1.keep[,cyto], arcsin[cyto-3,2])
  }
  
  #split per cell population
  a1.pops <- split(a1.keep, f = a1.keep$Population)
  b1.pops <- split(b1.keep, f = b1.keep$Population)
  
  #make for CCL2 and IL-1b the cutoff based on cd4 T cells
  a1.cd4n <- a1.keep[a1.keep$Population == 'CD4T' & grepl('Medium', a1.keep$Stimulus),]
  a1.ccl2_neg <- quantile(a1.cd4n[,'CCL2'], probs = .95)*2.2 #make sure to write correct CCL2
  a1.il1b_neg <- quantile(a1.cd4n[,'IL_1b'], probs = .95)*2.2 #make sure to write correct IL-1b
  
  b1.cd4n <- b1.keep[b1.keep$Population == 'CD4T' & grepl('Medium', b1.keep$Stimulus),]
  b1.ccl2_neg <- quantile(b1.cd4n[,'CCL2'], probs = .95)*2.2 #make sure to write correct CCL2_MCP_1 and GMCSF
  b1.il1b_neg <- quantile(b1.cd4n[,'IL_1b'], probs = .95)*2.2 #make sure to write correct CCL2_MCP_1 and GMCSF
  

    #dot plots and frequencies
for(m in 1:length(a1.pops)){ #cycle through the cell pops
  
  #select populations
  a.pop <- a1.pops[[m]] #select one cell population
  
  if(nrow(a.pop)>100){  #only work with populations that have at least 100 cells
  name <- names(a1.pops)[m] #make name
  
  b.pop <- b1.pops[[name]] #select one cell population
  
  
  #set negative gates
  medium.a <- a.pop[grepl('Medium', a.pop$Stimulus),] #select medium
  a1.neg <- data.frame(ncol=2, nrow=length(cytokines))
  for(j in 4:ncol(medium.a)){
    a1.neg[j-3,1] <- colnames(medium.a)[j]
    a1.neg[j-3,2] <- quantile(medium.a[,j], probs = .95)*2 #set all cutoffs at 1.2x
    
    if(colnames(medium.a)[j] == "CCL2"){a1.neg[j-3,2] <- a1.ccl2_neg} #for CCL2 add the values from CD4 negative
    if(colnames(medium.a)[j] == "IL_1b"){a1.neg[j-3,2] <- a1.il1b_neg} #for IL1b add the values from CD4 negative
    
    if(colnames(medium.a)[j] == "IFN_a"){a1.neg[j-3,2] <- 2} #for IFNa and TNF set at 2
    if(colnames(medium.a)[j] == "TNFa"){a1.neg[j-3,2] <- quantile(medium.a[,j], probs = .95)*2} #for IFNa and TNF set at 2
    
    if(colnames(medium.a)[j] == "IL_12"){a1.neg[j-3,2] <- 
      quantile(medium.a[,j], probs = .95)*2.5} 
    if(colnames(medium.a)[j] == "IL_10"){a1.neg[j-3,2] <- 
      quantile(medium.a[,j], probs = .95)*2.5} 
    if(colnames(medium.a)[j] == "IL_6"){a1.neg[j-3,2] <- 
      quantile(medium.a[,j], probs = .95)*2.7} 
    if(colnames(medium.a)[j] == "IL4_5_13"){a1.neg[j-3,2] <- 
      quantile(medium.a[,j], probs = .95)*2.5} 
    
 }
  colnames(a1.neg) <- c('Cytokine', 'Limit')
  
  #set max for IL8 and IL12
  if(a1.neg[2,2] > 2){a1.neg[2,2] <- 2}
  if(a1.neg[3,2] > 2.5){a1.neg[2,2] <- 2.5}
  
  #do the same for the short term stimulations
  medium.b <- b.pop[grepl('Medium', b.pop$Stimulus),] #select medium
  b1.neg <- data.frame(ncol=2, nrow=length(cytokines))
  for(j in 4:ncol(medium.b)){
    b1.neg[j-3,1] <- colnames(medium.b)[j]
    b1.neg[j-3,2] <- quantile(medium.b[,j], probs = .95)*2 #set all cutoffs at 1.2x
    
    if(colnames(medium.b)[j] == "CCL2"){b1.neg[j-3,2] <- b1.ccl2_neg} #for CCL2 add the values from CD4 negative
    if(colnames(medium.b)[j] == "IL_1b"){b1.neg[j-3,2] <- b1.il1b_neg} #for IL1b add the values from CD4 negative
    
    if(colnames(medium.b)[j] == "IFN_a"){b1.neg[j-3,2] <- 2} #for IFNa and TNF set at 2
    if(colnames(medium.b)[j] == "TNFa"){b1.neg[j-3,2] <- quantile(medium.a[,j], probs = .95)*2} #for IFNa and TNF set at 2
    
    if(colnames(medium.b)[j] == "IL_12"){b1.neg[j-3,2] <- 
      quantile(medium.b[,j], probs = .95)*2.5} 
    if(colnames(medium.b)[j] == "IL_10"){b1.neg[j-3,2] <- 
      quantile(medium.b[,j], probs = .95)*2.5} 
    if(colnames(medium.b)[j] == "IL_6"){b1.neg[j-3,2] <- 
      quantile(medium.b[,j], probs = .95)*2.7} 
    if(colnames(medium.b)[j] == "IL4_5_13"){b1.neg[j-3,2] <- 
      quantile(medium.b[,j], probs = .95)*2.5} 
  }
  colnames(b1.neg) <- c('Cytokine', 'Limit')
  
  #set max for IL8 and IL12
  if(b1.neg[2,2] > 2){b1.neg[2,2] <- 2}
  if(b1.neg[3,2] > 2.5){b1.neg[2,2] <- 2.5}
  

  dich.a <- a.pop
  for(j in 4:ncol(dich.a)){
    dich.a[j][dich.a[j] < a1.neg[j-3,2]] <- 0
    dich.a[j][dich.a[j] >= a1.neg[j-3,2]] <- 1
  }
  
  dich.b <- b.pop
  for(j in 4:ncol(dich.b)){
    dich.b[j][dich.b[j] < b1.neg[j-3,2]] <- 0
    dich.b[j][dich.b[j] >= b1.neg[j-3,2]] <- 1
  }
  
  dich <- rbind(dich.a, dich.b)
  dich <- dich[rowSums(dich[,4:14]) > 0,]  
  
  if(m == 1){
    dich2 <- dich
  } else {
    dich2 <- rbind(dich2, dich)
  }
  
  }
}
  
  
write.table(dich2, paste0('Dichotomous_', sampleid, '.txt'), 
            sep = '\t', row.names = F)


}
    


