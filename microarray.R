BiocManager::install(c('beadarray', 'GEOquery', 'illuminaHumanv1.db',
                       'illuminaHumanv2.db', 'illuminaHumanv3.db', 'BeadArrayUseCases',
                       'GOstats','GenomicRanges','Biostrings'))
BiocManager::install('GEOquery')
BiocManager::install('enrichplot')
BiocManager::install('tidyverse')
install.packages('pheatmap')
install.packages('enrichplot')


library(GEOquery)
library(enrichplot)
library(DOSE)
library(c('beadarray', 'GEOquery', 'illuminaHumanv1.db',
          'illuminaHumanv2.db', 'illuminaHumanv3.db', 'BeadArrayUseCases',
          'GOstats','GenomicRanges','Biostrings'))
library(beadarray)
library(GEOquery)
library(illuminaHumanv1.db)
library(illuminaHumanv2.db)
library(illuminaHumanv3.db)
library(BeadArrayUseCases)
library(GOstats)
library(GenomicRanges)
library(Biostrings)

# import packages
library(limma)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(dplyr)


raw_data <- read.table('SCP2_Microarray.txt', header = T, row.names = 1, stringsAsFactors = F)
log2_data <- log2(raw_data)
data_norm <- as.matrix(normalizeQuantiles(log2_data))
boxplot(data_norm, outline=F)
design_mat <- cbind(1, c(0,0,1,1))
fit <- lmFit(data_norm, design_mat)
fit <- eBayes(fit)
tt <- topTable(fit, coef=2, adjust.method = 'BH', sort.by = 'p', number =50000)

col_break = seq(0.990, 1, by=0.001)
pheatmap(
  cor(data_norm), show_rownames=F, breaks = col_break, 
         color = colorRampPalette(c('blue', 'white', 'red'))(length(col_break)))
write.csv(tt, 'microarray.csv')
write.csv(data_norm, 'normalized_data.csv')

tt$log_trans <- -log10(tt$adj.P.Val)

tt <- tt %>%
  mutate(color=case_when(logFC<0 & adj.P.Val<0.05 ~ 'blue',
                         logFC>0 & adj.P.Val<0.05 ~ 'red',
                         TRUE ~ 'black'))



# volcano plot
volcano_plot <- 
ggplot(tt, aes(x=logFC, y=log_trans, color=color))+
  geom_point(size= 3, shape=16, position = 'jitter', alpha=0.6) +
  scale_color_manual(values = c('black', 'blue', 'red')) +
  geom_hline(yintercept = 1.30102999566398, size=1, linetype='dashed') + 
  xlim(-3.5,3.5) +
  theme_classic()

volcano_plot
png(filename = 'volcano_plot.png', res = 600, width = 6000, height = 3600)
jpeg(filename = 'volcano_plot.jpeg', res = 600, width = 6000, height = 3600)
volcano_plot
dev.off()




volcano_plot = 
  ggplot(tt, aes(x=logFC, y=log_trans, fill=color))+
  geom_point(size= 2, shape=16, position = 'jitter')+
  # geom_hline(yintercept= 2.65351010743292, color='gray', size=1)+
  # geom_vline(xintercept = c(-1.5,1.5), color='gray', size=1)+
  geom_point(data=tt[tt$logFC<0 & tt$adj.P.Val<0.05,], color='blue', shape=16, size=4) +
  geom_point(data=tt[tt$logFC>0 & tt$adj.P.Val<0.05,], color='red', shape=16, size=4) +
  xlab('Log2 Fold Change (LFC)')+
  ylab('Log10(FDR)')+
  geom_hline(yintercept = 1.30102999566398, size=1, linetype='dashed') + 
  theme_classic()

volcano_plot
png(filename = 'volcano_plot.png', res = 600, width = 4800, height = 3600)
volcano_plot
dev.off()

# GO biological process input and collapse the GO terms and pool the gene members. 
GO_BP <- read.table('GO_Biological_Process_2018_table.txt', header = F, sep = '\t', stringsAsFactors = F)
GO_BP_simple <- GO_BP[,c(1,9)]
# head(GO_BP_simple$Genes)
# typeof(GO_BP_simple$Genes[1])
split_fun <- function(x) {
  unlist(strsplit(x, split=';'))
}
GO_BP_simple$Comp <- lapply(GO_BP_simple$V9, split_fun)

# input normalized data
normalized_data <- read.csv('normalized_data.csv', header = T, row.names = 1)

# cell migration
GO_0030334 <- GO_BP_simple$Comp[[1]]
GO_0030335 <- GO_BP_simple$Comp[[4]]
GO_0030336 <- GO_BP_simple$Comp[[14]]
GO_2000147 <- GO_BP_simple$Comp[[3]]
GO_0014910 <- GO_BP_simple$Comp[[19]]
cell_migration <- unique(c(GO_0014910, GO_2000147, GO_0030336, GO_0030335, GO_0030334))

cell_migration_data <- normalized_data[match(cell_migration, rownames(normalized_data)),]
head(cell_migration_data)
pheatmap(cell_migration_data, scale = 'row', cluster_cols = F, border_color = NA)

# cell proliferation
GO_0048661 <- GO_BP_simple$Comp[[5]]
GO_0048660 <- GO_BP_simple$Comp[[11]]
GO_0043066 <- GO_BP_simple$Comp[[30]]
GO_0045931 <- GO_BP_simple$Comp[[17]]
GO_0008284 <- GO_BP_simple$Comp[[9]]
GO_0042127 <- GO_BP_simple$Comp[[6]]
GO_0048522 <- GO_BP_simple$Comp[[25]]
cell_proliferation <- unique(c(GO_0048661, GO_0048660, GO_0043066, GO_0045931, GO_0008284, GO_0042127, GO_0048522))

cell_proliferation_data <- normalized_data[match(cell_proliferation, rownames(normalized_data)),]
pheatmap(cell_proliferation_data, scale = 'row', cluster_cols = F, border_color = NA)

# cell differentiation 
GO_0045596 <- GO_BP_simple$Comp[[8]]
GO_0072182 <- GO_BP_simple$Comp[[26]]
cell_differentiation <- unique(c(GO_0045596, GO_0072182, 'BMP4'))
cell_differentiation_data <- normalized_data[match(cell_differentiation, rownames(normalized_data)),]
pheatmap(cell_differentiation_data, scale = 'row', cluster_cols = F, border_color = NA)

# cell adhesion
GO_0022409 <- GO_BP_simple$Comp[[18]]
GO_0007162 <- GO_BP_simple$Comp[[12]]
GO_0033629 <- GO_BP_simple$Comp[[16]]
cell_adhesion <- unique(c(GO_0022409, GO_0007162, GO_0033629))
cell_adhesion_data <- normalized_data[match(cell_adhesion, rownames(normalized_data)),]
pheatmap(cell_adhesion_data, scale = 'row', cluster_cols = F, border_color = NA)

# extracellular matrix organization 
GO_0030198 <- GO_BP_simple$Comp[[2]]
ECM <- GO_0030198
ECM_data <- normalized_data[match(ECM, rownames(normalized_data)),]
pheatmap(ECM_data, scale = 'row', cluster_cols = F, border_color = NA)

# blood coagulation
GO_0030195 <- GO_BP_simple$Comp[[7]]
GO_0030193 <- GO_BP_simple$Comp[[13]]
blood_coagulation <- unique(c(GO_0030195, GO_0030193))
blood_coagulation_data <- normalized_data[match(blood_coagulation, rownames(normalized_data)),]
pheatmap(blood_coagulation_data, scale = 'row', cluster_cols = F, border_color = NA)

# signal transduction
GO_0009968 <- GO_BP_simple$Comp[[29]]
GO_0090288 <- GO_BP_simple$Comp[[20]]
GO_0048017 <- GO_BP_simple$Comp[[27]]
GO_0040036 <- GO_BP_simple$Comp[[23]]
GO_0048015 <- GO_BP_simple$Comp[[21]]
GO_0043405 <- GO_BP_simple$Comp[[10]]
GO_0043408 <- GO_BP_simple$Comp[[24]]
signal_transduction <- unique(c(GO_0009968, GO_0090288, GO_0048017, GO_0040036, GO_0048015, GO_0043405, GO_0043408))
signal_transduction_data <- normalized_data[match(signal_transduction, rownames(normalized_data)),]
pheatmap(signal_transduction_data, scale = 'row', cluster_cols = F, border_color = NA)

# endoderm 
GO_0001655 <- GO_BP_simple$Comp[[28]]
GO_0001706 <- GO_BP_simple$Comp[[22]]
endoderm <- unique(c(GO_0001655, GO_0001706))
endoderm_data <- normalized_data[match(endoderm, rownames(normalized_data)),]
pheatmap(endoderm_data, scale = 'row', cluster_cols = F, border_color = NA)

# peptidyl-tyrosine modification
GO_0018212 <- GO_BP_simple$Comp[[15]]
modification <- GO_0018212
modification_data <- normalized_data[match(modification, rownames(normalized_data)),]
pheatmap(modification_data, scale = 'row', cluster_cols = F, border_color = NA)




# a <- unlist(strsplit(GO_BP_simple$Genes[1], split=';'))
# input bmp family
BMPs_data <- read.table('BMPs.txt', header = T, row.names = 1, stringsAsFactors = F)
pheatmap(BMPs_data, cluster_cols = T, cluster_rows = T , show_colnames = T, show_rownames = T, scale = 'column')

# 






design <- matrix(nrow=18,ncol=2)
design[,1] <- c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1)
design[,2] <- c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0)
colnames(design) <- c("Normal","Tumour")
