rm(list=ls())

#packages that are needed
library(reshape2) 
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(plyr)

#function
arcsinhn <- function(x, n){
  return(asinh(x/n))  
}

####################################################################
#set the working directory
setwd('C:\\Users\\spjochems\\Desktop\\Current\\ICS\\Output FCS/Annotated')

#read in data and fix column names
files <- list.files()

#read in the first file with 100,000 cells
n=100000
a <-  read.delim(files[1], nrows = n)

for(i in 2:length(files)){
  b <-  read.delim(files[i], nrows = n)
  a <- rbind(a,b)
}


a$Sample <- rep(files, each=n)


a.sub <- a[sample(1:nrow(a), 1000000),]
#get the cytokine names
markers <- read.delim('R:\\Para-CIH\\CIH Group member folders\\Simon\\Results\\flow\\LUF Aurora ICS\\2022_01_17 18 samples new\\panel.txt')
cytokines <- markers$Maker[markers$Type == 'Cytokine']


a.sub <- a.sub[,c('CD14', cytokines)]

a150 <- a.sub
a150[,2:12] <- arcsinhn(a150[,2:12], 150)
a150[,1] <- arcsinhn(a150[,1], 1500)

a500 <- a.sub
a500[,2:12] <- arcsinhn(a500[,2:12], 500)  
a500[,1] <- arcsinhn(a500[,1], 1500)

a1500 <- a.sub
a1500[,2:12] <- arcsinhn(a1500[,2:12], 1500)  
a1500[,1] <- arcsinhn(a1500[,1], 1500)  

a5000 <- a.sub
a5000[,2:12] <- arcsinhn(a5000[,2:12], 5000)  
a5000[,1] <- arcsinhn(a5000[,1], 1500)  

a15000 <- a.sub
a15000[,2:12] <- arcsinhn(a15000[,2:12], 15000)  
a15000[,1] <- arcsinhn(a15000[,1], 1500)  


pdf('arcsin transforms.pdf', onefile = T, width = 10, height = 10)

for(i in 2:ncol(a.sub)){
  
name <- colnames(a.sub)[i]  
colnames(a150)[i] <- colnames(a500)[i] <- colnames(a1500)[i] <- 
  colnames(a5000)[i] <- colnames(a15000)[i] <- 'Marker' 

p1 <- ggplot(a150, aes(x=CD14, y=Marker)) + geom_hex(bins = 128)+
  scale_fill_viridis_c() + ggtitle('150') + ylab(name) + 
  theme(legend.position = 'none')
p2 <- ggplot(a500, aes(x=CD14, y=Marker)) + geom_hex(bins = 128)+
  scale_fill_viridis_c() + ggtitle('500') + ylab(name)+ 
  theme(legend.position = 'none')
p3 <- ggplot(a1500, aes(x=CD14, y=Marker)) + geom_hex(bins = 128)+
  scale_fill_viridis_c() + ggtitle('1500') + ylab(name)+ 
  theme(legend.position = 'none')
p4 <- ggplot(a5000, aes(x=CD14, y=Marker)) + geom_hex(bins = 128)+
  scale_fill_viridis_c() + ggtitle('5000') + ylab(name)+ 
  theme(legend.position = 'none')
p5 <- ggplot(a15000, aes(x=CD14, y=Marker)) + geom_hex(bins = 128)+
  scale_fill_viridis_c() + ggtitle('15000') + ylab(name)+ 
  theme(legend.position = 'none')

print(plot_grid(p1,p2,p3, p4, p5, ncol = 3))

colnames(a150)[i] <- colnames(a500)[i] <- colnames(a1500)[i] <- 
  colnames(a5000)[i] <- colnames(a15000)[i] <- name 

}

dev.off()

