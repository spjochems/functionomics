rm(list = ls())

library(reshape2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
library(pheatmap)
library(ggrepel)
library('limma')
library(ggtern)
library(writexl)

setwd('V:\\para-yazdanbakhsh\\JochemsLab\\Simon\\LUF aurora panel\\2022_01_17 new comp and gating checked\\analysis\\Output FCS\\Cytokine gates\\')

#load in data
files <- list.files(pattern = 'Freqs')
files.read <- lapply(files, read.delim)

#get for one sample the id
freqs <- files.read[[1]]
colnames(freqs)[5] <- freqs$Sample[1]
freqs$ID <- paste(freqs$Population, freqs$Stimulation, freqs$Cytokine, sep = '_')
freqs <- freqs[,-1]




#freqs <- freqs[freqs$Stimulation != '0.Medium',]

#merge all files together
for(i in 2:length(files)){
    freqs2 <- files.read[[i]]
    colnames(freqs2)[5] <- freqs2$Sample[1]
    
    #freqs2 <- freqs2[freqs2$Stimulation != '0.Medium',]
    
    freqs2$ID <- paste(freqs2$Population, freqs2$Stimulation, freqs2$Cytokine, sep = '_')
    
    freqs <- merge(freqs, freqs2[,5:6], by = 'ID', all = T)
}


#remove cDC1 cells with too few cells to be able to interpret
freqs <- freqs[!is.na(freqs$Population) & 
                 freqs$Population != 'Monocytes' &
                 freqs$Population != 'NK',]


write.table(freqs, 'frequencies_ICS_all.txt', sep='\t', row.names=T)
write_xlsx(freqs, 'frequencies_ICS_all.xlsx')

dir.create('Analysis results_20240801_co0_5')
setwd('Analysis results_20240801_co0_5')

#order donors
ord <- c('449','7335','6106','55','3460','3713','D5','D17','D18','D19','D20',
         'Angga','D4','SJ','Wouter','JS','D7','EUR6')
freqs3 <- freqs[,append(colnames(freqs)[1:4], ord)]


freqs3 <- freqs3[rowSums(is.na(freqs3[,5:22])) <10,] #remove the ones with too many missing stuff


#get average
freqs3$Average <- apply(freqs3[,5:22], 1, mean, na.rm=T)
freqs3$AverageRural <- apply(freqs3[,5:10], 1, mean, na.rm=T)
freqs3$AverageUrban <- apply(freqs3[,11:16], 1, mean, na.rm=T)
freqs3$AverageEur <- apply(freqs3[,17:22], 1, mean, na.rm=T)

#add filter
freqs3$Sufficient <- 'No'
freqs3$Sufficient[freqs3$Average > 0.5] <- 'Yes'
sum(freqs3$Sufficient == 'Yes')/nrow(freqs3)*100

#order by freq and plot
freqs3$ID <- factor(freqs3$ID, levels = freqs3$ID[order(freqs3$Average)])
ggplot(freqs3, aes(y=ID, x=Average, colour = Sufficient))+ 
    geom_bar(stat='identity') + 
    theme(axis.text.y=element_blank()) + 
    xlab('Average expression') + 
    ggtitle(paste0(sum(freqs3$Sufficient == 'Yes'), ' features kept at 0.5% out of ',
                   nrow(freqs3)))
ggsave('Freqs_to_keep_co0,5.pdf', width = 4, height = 5)


#calculate differences
freqs3$EU_Rur_Diff <- freqs3$AverageEur / freqs3$AverageRural
freqs3$EU_Urb_Diff <- freqs3$AverageEur / freqs3$AverageUrban
freqs3$Urb_Rur_Diff <- freqs3$AverageUrban / freqs3$AverageRural


#get pvalues
freqs3$EU_Rur_pval <- NA
freqs3$EU_Urb_pval <- NA
freqs3$Urb_Rur_pval <- NA

#select significant ones
freqs4 <- freqs3[freqs3$Sufficient == 'Yes',]


#use limma analysis
freqs.limma <- freqs4[,5:22]
rownames(freqs.limma) <- freqs4$ID
design <- model.matrix(~ 0+factor(c(rep(1,6), rep(2,6), rep(3,6))))
colnames(design) <- c("Rural", "Urban", "European")
fit <- lmFit(freqs.limma, design)
contrast.matrix <- makeContrasts(Rural-Urban, Rural-European, Urban-European, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef=1, adjust="fdr")
topTable(fit2, coef=2, adjust="fdr")
topTable(fit2, coef=3, adjust="fdr")

#g <- decideTests(fit2,method="global",adjust.method="none", lfc = 0)

#add the data to the freqs4. Multiple by 3 for 3 comparisons
freqs4$Urb_Rur_pval <- topTable(fit2, coef=1, adjust="fdr",sort.by = 'none', number = nrow(freqs4))$P.Value*3
freqs4$EU_Rur_pval <- topTable(fit2, coef=2, adjust="fdr",sort.by = 'none', number = nrow(freqs4))$P.Value*3
freqs4$EU_Urb_pval <- topTable(fit2, coef=3, adjust="fdr",sort.by = 'none', number = nrow(freqs4))$P.Value*3

freqs4$Urb_Rur_pval[freqs4$Urb_Rur_pval > 1] <- 1
freqs4$EU_Rur_pval[freqs4$EU_Rur_pval > 1] <- 1
freqs4$EU_Urb_pval[freqs4$EU_Urb_pval > 1] <- 1

#add the BH corrected pvalue per population
freqs.split <- split(freqs4, f=freqs4$Population)
freqs.list <- lapply(freqs.split, function(i){
  i$Urb_Rur_padj <- p.adjust(i$Urb_Rur_pval, method = 'BH') 
  i$EU_Rur_padj <- p.adjust(i$EU_Rur_pval, method = 'BH') 
  i$EU_Urb_padj <- p.adjust(i$EU_Urb_pval, method = 'BH') 
  i <- i
})
freqs5 <- do.call(rbind, freqs.list)

write.table(freqs5, 'frequencies_ICS.txt', sep = '\t', row.names = F)
write_xlsx(freqs5[,-(34:36)], 'frequencies_ICS.xlsx')





freqs.sig <- freqs4[(freqs4$EU_Rur_pval < 0.05 | 
                       freqs4$EU_Urb_pval < 0.05 |
                       freqs4$Urb_Rur_pval < 0.05) ,]
freqs.sig <- freqs.sig[!is.na(freqs.sig$ID),]


#plot of these ones
freqs.melt <- melt(freqs.sig[,1:22], id = colnames(freqs.sig)[1:4])

meta <- readxl::read_xlsx('V:\\para-yazdanbakhsh\\JochemsLab\\Simon\\R drive copied 2022\\LUF Aurora ICS\\2022_01_17 18 samples new\\ICS-metadata.xlsx', 1)


freqs.melt <- merge(meta[,1:7], freqs.melt, by.x='ID', by.y='variable')
pops.ord <- freqs.melt$ID.y[order(freqs.melt$Population)]
freqs.melt$ID.y <- factor(freqs.melt$ID.y, 
                          levels = pops.ord[!duplicated(pops.ord)])
ggplot(freqs.melt, aes(x=group, y=value, colour = group)) + 
    geom_boxplot(outlier.colour = NA) + 
    geom_jitter(width = 0.1) +
    theme_bw() + 
    facet_wrap(.~ID.y, scales = 'free') + 
    theme(strip.text.x = element_text(size = 8))
ggsave('sig_boxplots.pdf', width = 12, height = 12)


#heatmap all
hm <- t(scale(t(freqs.sig[,5:22]))) 
hm[hm < -3] <- -3
hm[hm > 2] <- 2
rownames(hm) <- freqs.sig$ID
dev.off()
pdf('heatmap_sig_diffs.pdf', width = 8, height = 7)
pheatmap(hm)
dev.off()

#volcano
volc <- freqs4
volc$EU_Rur_sig <- 'No'
volc$EU_Rur_sig[volc$EU_Rur_pval<0.05] <- 'Yes'
volc$EU_Rur_lab <- volc$ID
volc$EU_Rur_lab[!volc$EU_Rur_pval<0.05] <- NA
pos <- sum(volc$EU_Rur_sig == 'Yes' & log2(volc$EU_Rur_Diff)>0)
neg <- sum(volc$EU_Rur_sig == 'Yes' & log2(volc$EU_Rur_Diff)<0)
ggplot(volc, aes(x = log2(EU_Rur_Diff), y = -log10(EU_Rur_pval))) + 
    geom_point(aes(colour = EU_Rur_sig)) + theme_bw() + 
    geom_text_repel(aes(label = EU_Rur_lab), size = 3) + 
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
    ggtitle('Rural vs European') + theme(aspect.ratio = 1)+
    scale_colour_manual(values = c('grey', 'blue')) + 
    theme(legend.position = 'none') +
    annotate('text', x = Inf, y = Inf, hjust = 2, vjust = 2, label = pos, colour = 'red') +
    annotate('text', x = -Inf, y = Inf, hjust = -2, vjust = 2, label = neg, colour = 'red')
ggsave('volcano_EU_Rur.pdf', width = 10, height = 10)



p1<-ggplot(volc, aes(x = log2(EU_Rur_Diff), y = -log10(EU_Rur_pval))) + 
  geom_point(aes(colour = EU_Rur_sig, size = EU_Rur_sig)) + theme_bw() + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  ggtitle('Rural vs European') + theme(aspect.ratio = 1)+
  scale_colour_manual(values = c('grey', 'blue')) + 
  theme(legend.position = 'none') +
  ylim(c(0,3)) +
  xlim(c(-5,5))+
  annotate('text', x = Inf, y = Inf, hjust = 2, vjust = 2, label = pos, colour = 'red') +
  annotate('text', x = -Inf, y = Inf, hjust = -2, vjust = 2, label = neg, colour = 'red')+ 
  scale_size_manual(values = c(1,2))
p1

volc$EU_Urb_sig <- 'No'
volc$EU_Urb_sig[volc$EU_Urb_pval<0.05] <- 'Yes'
volc$EU_Urb_lab <- volc$ID
volc$EU_Urb_lab[!volc$EU_Urb_pval<0.05] <- NA
pos <- sum(volc$EU_Urb_sig == 'Yes' & log2(volc$EU_Urb_Diff)>0)
neg <- sum(volc$EU_Urb_sig == 'Yes' & log2(volc$EU_Urb_Diff)<0)

ggplot(volc, aes(x = log2(EU_Urb_Diff), y = -log10(EU_Urb_pval))) + 
    geom_point(aes(colour = EU_Urb_sig)) + theme_bw() + 
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed')+
  ggtitle('Urban vs European') + theme(aspect.ratio = 1)+
  scale_colour_manual(values = c('grey', 'blue')) + 
  theme(legend.position = 'none') +
  annotate('text', x = Inf, y = Inf, hjust = 2, vjust = 2, label = pos, colour = 'red') +
  annotate('text', x = -Inf, y = Inf, hjust = -2, vjust = 2, label = neg, colour = 'red')
ggsave('volcano_EU_Urb.pdf', width = 10, height = 10)



p2<-ggplot(volc, aes(x = log2(EU_Urb_Diff), y = -log10(EU_Urb_pval))) + 
  geom_point(aes(colour = EU_Urb_sig, size = EU_Urb_sig)) + theme_bw() + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  ggtitle('Urban vs European') + theme(aspect.ratio = 1)+
  scale_colour_manual(values = c('grey', 'blue')) + 
  theme(legend.position = 'none') +
  ylim(c(0,3)) +  xlim(c(-5,5))+
  annotate('text', x = Inf, y = Inf, hjust = 2, vjust = 2, label = pos, colour = 'red') +
  annotate('text', x = -Inf, y = Inf, hjust = -2, vjust = 2, label = neg, colour = 'red')+ 
  scale_size_manual(values = c(1,2))
p2


volc$Urb_Rur_sig <- 'No'
volc$Urb_Rur_sig[volc$Urb_Rur_pval<0.05] <- 'Yes'
volc$Urb_Rur_lab <- volc$ID
volc$Urb_Rur_lab[!volc$Urb_Rur_pval<0.05] <- NA
pos <- sum(volc$Urb_Rur_sig == 'Yes' & log2(volc$Urb_Rur_Diff)>0)
neg <- sum(volc$Urb_Rur_sig == 'Yes' & log2(volc$Urb_Rur_Diff)<0)

ggplot(volc, aes(x = log2(Urb_Rur_Diff), y = -log10(Urb_Rur_pval))) + 
    geom_point(aes(colour = Urb_Rur_sig)) + theme_bw() + 
  geom_text_repel(aes(label = Urb_Rur_lab), size = 3) + 
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed')+
  ggtitle('Rural vs Urban') + theme(aspect.ratio = 1)+
  scale_colour_manual(values = c('grey', 'blue')) + 
  theme(legend.position = 'none') +
  annotate('text', x = Inf, y = Inf, hjust = 2, vjust = 2, label = pos, colour = 'red') +
  annotate('text', x = -Inf, y = Inf, hjust = -2, vjust = 2, label = neg, colour = 'red')
ggsave('volcano_Urb_Rur.pdf', width = 10, height = 10)

p3<-ggplot(volc, aes(x = log2(Urb_Rur_Diff), y = -log10(Urb_Rur_pval))) + 
  geom_point(aes(colour = Urb_Rur_sig, size = Urb_Rur_sig)) + theme_bw() + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  ggtitle('Rural vs Urban') + theme(aspect.ratio = 1)+
  scale_colour_manual(values = c('grey', 'blue')) + 
  theme(legend.position = 'none') +
  ylim(c(0,3)) +  xlim(c(-5,5))+
  annotate('text', x = Inf, y = Inf, hjust = 2, vjust = 2, label = pos, colour = 'red') +
  annotate('text', x = -Inf, y = Inf, hjust = -2, vjust = 2, label = neg, colour = 'red')+ 
  scale_size_manual(values = c(1,2))
p3


plot_grid(p1,p2,p3, ncol=3)
ggsave('volcanoes.pdf', width = 9, height = 3)


tern <- freqs5
tern$Name <- tern$ID
tern$Sig <- 'No'
tern$Sig[tern$ID %in% freqs.sig$ID] <- 'Yes'
tern$Name[!tern$ID %in% freqs.sig$ID] <- NA


t1 <- ggtern(tern, aes(x=AverageEur, y=AverageUrban, z=AverageRural)) + 
  geom_point(aes(colour = Cytokine, alpha = Sig, size = Sig)) + 
  scale_alpha_manual(values = c(0.2, 1)) + 
  scale_size_manual(values = c(1.5,2)) + 
  theme_bw() +
  theme_nolabels() + 
  xlab('EU') + ylab('Urban') + zlab('Rural') + 
  ggtitle('Cytokine') + 
  theme(legend.position = 'bottom', 
        legend.text=element_text(size=6),
        legend.title = element_blank()) +
  guides(alpha='none', size='none')

t2 <- ggtern(tern, aes(x=AverageEur, y=AverageUrban, z=AverageRural)) + 
  geom_point(aes(colour = Population, alpha = Sig, size = Sig)) + 
  scale_alpha_manual(values = c(0.2, 1)) + 
  scale_size_manual(values = c(1.5,2)) + 
  theme_bw() +
  theme_nolabels() + 
  xlab('EU') + ylab('Urban') + zlab('Rural')+ 
  ggtitle('Population')+ 
  theme(legend.position = 'bottom', 
        legend.text=element_text(size=6),
        legend.title = element_blank()) +
  guides(alpha='none', size='none')

t3 <- ggtern(tern, aes(x=AverageEur, y=AverageUrban, z=AverageRural)) + 
  geom_point(aes(colour = Stimulation, alpha = Sig, size = Sig)) + 
  scale_alpha_manual(values = c(0.2, 1)) + 
  scale_size_manual(values = c(1.5,2)) + 
  theme_bw() +
  theme_nolabels() + 
  xlab('EU') + ylab('Urban') + zlab('Rural')+ 
  ggtitle('Stimulation')+ 
  theme(legend.position = 'bottom', 
        legend.text=element_text(size=6),
        legend.title = element_blank()) +
  guides(alpha='none', size='none')

pdf('ternary_all.pdf', width = 14, height = 5)
print(grid.arrange(t1,t2,t3, ncol = 3))
dev.off()


ggtern(tern, aes(x=AverageEur, y=AverageUrban, z=AverageRural)) + 
  geom_point(aes(colour = Cytokine, size = Sig, alpha = Sig)) + 
  facet_wrap(.~Population) + 
  scale_alpha_manual(values = c(0.5,1)) + 
  scale_size_manual(values = c(1,1.5)) + 
  theme_bw() +
  theme_nolabels() + 
  xlab('EU') + ylab('Urban') + zlab('Rural')
ggsave('ternary_per_pop.pdf', width = 10, height = 10)



#plot per population
freqs.melt2 <- melt(freqs, id = colnames(freqs)[1:4])
freqs.melt2 <- merge(meta[,1:7], freqs.melt2, by.x='ID', by.y='variable')


#add level above 1%
freqs.melt2$Suff <- 'No'
freqs.melt2$Suff[freqs.melt2$ID.y %in% freqs4$ID] <- 'Yes'

freqs.melt2$Sig <- 'No'
freqs.melt2$Sig[freqs.melt2$ID.y %in% freqs.sig$ID] <- 'Yes'


freqs.split <- split(freqs.melt2, f=freqs.melt2$Population)

pdf('Populations_freqs.pdf', width = 20, height = 10, onefile = T)

for(i in 1:length(freqs.split)){

pop <- freqs.split[[i]]
name <- names(freqs.split)[i]


p1 <- ggplot(pop, aes(y=Stimulation, x=value,  group = group, fill = group)) + 
  facet_grid(.~Cytokine) +  
  geom_bar(fun = "mean", stat = "summary", position = position_dodge(), 
           aes(alpha = Sig)) + 
  geom_point(position = position_dodge(width = 1), aes(size = Suff)) + 
  theme_bw() + 
  xlab('% of population') + 
  ylab('') + 
  ggtitle(name) + 
  scale_size_manual(values = c(0,1)) + 
  scale_alpha_manual(values = c(0.5, 1))
print(p1)


p1 <- ggplot(pop, aes(y=Stimulation, x=value,  group = group, fill = group)) + 
  facet_grid(.~Cytokine) +  
  geom_bar(fun = "mean", stat = "summary", position = position_dodge()) + 
  stat_summary(fun = mean,
               geom = "errorbar",
               width = 0.3,
               position = position_dodge(width = 1),
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x))) +
  theme_bw() + 
  xlab('% of population') + 
  ylab('') + 
  ggtitle(name) + 
  scale_size_manual(values = c(0,1)) + 
  scale_alpha_manual(values = c(0.5, 1))
print(p1)


}
dev.off()

