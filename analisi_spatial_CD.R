library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)
library(ggdendroplot)
library(rstatix)
library(patchwork)
library(RColorBrewer)
library(vegan)
library(ggrepel)
library(DESeq2)
library(Seurat)
library(pheatmap)
library(ggpubr)

### remove genes expressed in lots of cells

metadata <- read.csv('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/R5941_ColonTMA_metadata_file.csv')

metadata$Patients <- ifelse(
  metadata$fov %in% c(39,40,41), "Patient_20", ifelse(
    metadata$fov %in% c(37,38), "Patient_19", ifelse(
      metadata$fov %in% c(35,36), "Patient_18", ifelse(
        metadata$fov %in% c(33,34), "Patient_17", ifelse(
          metadata$fov %in% c(31,32), "Patient_16", ifelse(
            metadata$fov %in% c(29,30), "Patient_15", ifelse(
              metadata$fov %in% c(27,28), "Patient_14", ifelse(
                metadata$fov %in% c(25,26), "Patient_13", ifelse(
                  metadata$fov %in% c(23,24), "Patient_12", ifelse(
                    metadata$fov %in% c(21,22), "Patient_11", ifelse(
                      metadata$fov %in% c(18,19,20), "Patient_9", ifelse(
                        metadata$fov %in% c(14,15,16,17), "Patient_8", ifelse(
                          metadata$fov %in% c(12,13), "Patient_7", ifelse(
                            metadata$fov %in% c(9,10,11), "Patient_6", ifelse(
                              metadata$fov %in% c(6,7,8), "Patient_4", ifelse(
                                metadata$fov %in% c(4,5), "Patient_2",
                                "Patient_1"))))))))))))))))

metadata$group <- ifelse(
  metadata$Patients %in% c('Patient_11', 'Patient_12', 'Patient_13', 'Patient_14',
                           'Patient_15', 'Patient_16', 'Patient_17', 'Patient_18',
                           'Patient_19', 'Patient_20'), "CD", "K")


metadata <- metadata[metadata$group == 'CD',]

write.csv(metadata, row.names = FALSE, file = '/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/progetto_alberto_baeri/spatial/CD/R5941_ColonTMA_metadata_file_CD.csv')


counts <- read.csv('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/progetto_alberto_baeri/spatial/R5941_ColonTMA_exprMat_file.csv')

counts <- counts[counts$fov %in% unique(metadata$fov),]
counts[,which(colnames(counts) == 'MALAT1')] <- NULL
counts[,which(colnames(counts) == 'MHC.I')] <- NULL
counts[,which(colnames(counts) == 'NEAT1')] <- NULL
counts[,which(colnames(counts) == 'RPL37')] <- NULL
counts[,which(str_detect(colnames(counts), 'RPL.*'))] <- NULL

'RGS13' %in% colnames(counts)

write.csv(counts, row.names = FALSE,file='/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/progetto_alberto_baeri/spatial/CD/R5941_ColonTMA_exprMat_file_CD.csv')

positions <- read.csv('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/progetto_alberto_baeri/spatial/R5941_ColonTMA_fov_positions_file.csv')

positions <- positions[positions$fov %in% unique(metadata$fov),]
write.csv(positions, row.names = FALSE, file= '/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/progetto_alberto_baeri/spatial/CD/R5941_ColonTMA_fov_positions_file_CD.csv')

## correlation matrix all clusters not named
counts_per_cluster <- data.frame(readxl::read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/gene_counts_all_clusters_non_named_after_subclustering.xlsx'))
rownames(counts_per_cluster) <- counts_per_cluster[,1]
counts_per_cluster[,1] <- NULL

colnames(counts_per_cluster) <- gsub('X', '', colnames(counts_per_cluster))

counts_per_cluster_z <- counts_per_cluster

for (r in 1:nrow(counts_per_cluster_z)) {
  data <- as.numeric(counts_per_cluster_z[r,])
  counts_per_cluster_z[r,] <- (data-mean(data))/sd(data)
}

colnames(counts_per_cluster_z) <- c('0 - B reg/memory', '1 - Epithelial', '2 - T', '3 - Smooth muscle', '4 - Plasma IgG', '5 - Plasma IgG', '6 - Fibroblast', '7 - Plasma IgA', '8 - B',
                                    '9 - Macrophage','10 - Endothelial','11 - Plasma IgG','12 - Plasma IgA','13 - Plasma IgG','14 - Smooth muscle', '15 - Mast', '16 - Epithelial',
                                    '17 - Enteric glia', '18 - Epithelial', '19 - Epithelial', '20 - Plasma IgA', '21 - Epithelial')

corr_clusters <- cor(counts_per_cluster_z, method = 'spearman')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/correlation_matrix_all_genes_numbered_clusters.pdf', height = 10, width = 12)
pheatmap(corr_clusters, cluster_cols=TRUE, cluster_rows=TRUE, show_rownames = TRUE,show_colnames = TRUE)
dev.off()   

######### correlation matrix named clusters

counts_per_cluster <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/gene_counts_all_clusters_named.xlsx'))
rownames(counts_per_cluster) <- counts_per_cluster[,1]
counts_per_cluster[,1] <- NULL

markers <- c('CD3D',	'CD3E',	'CD3G', 'CD2', 'FYN','CD19', 'MS4A1',	'CD44',
             'IGKC','IGHA1','IGHG1', 
             'CD68', 'CD163', 'CLEC10A', 'CPA3', 'KIT',
             'TPSAB1.B2', 'MYH11', 'ACTG2', 'ACTA2', 'TAGLN',
             'PECAM1', 'CD93', 'TIE1','EPCAM','COL1A1', 'DCN', 'VIM', 'S100B')
              

counts_per_cluster_markers <- counts_per_cluster[markers,]
counts_per_cluster_markers_z <- counts_per_cluster_markers

for (r in 1:nrow(counts_per_cluster_markers_z)) {
  data <- as.numeric(counts_per_cluster_markers_z[r,])
  counts_per_cluster_markers_z[r,] <- (data-mean(data))/sd(data)
}

counts_per_cluster_z <- counts_per_cluster

for (r in 1:nrow(counts_per_cluster_z)) {
  data <- as.numeric(counts_per_cluster_z[r,])
  counts_per_cluster_z[r,] <- (data-mean(data))/sd(data)
}

colnames(counts_per_cluster_markers_z) <- c("B"   ,"B memory/reg" ,"Endothelial" ,"Enteric glia cells", "Epithelial","Fibroblast"  , "Macrophage" ,      
                                            "Mast cells", "Plasma IgA", "Plasma IgG"   ,"Smooth muscle cells" ,  "T" )

colnames(counts_per_cluster_z) <- c("B"   ,"B memory/reg" ,"Endothelial" ,"Enteric glia cells", "Epithelial","Fibroblast"  , "Macrophage" ,      
                                    "Mast cells", "Plasma IgA", "Plasma IgG"   ,"Smooth muscle cells" ,  "T" )

counts_per_cluster_markers_z <- counts_per_cluster_markers_z[,c("B"   ,"B memory/reg" ,"Endothelial" ,"Enteric glia cells", "Epithelial","Fibroblast"  , "Macrophage" ,      
                                                                "Mast cells", "Plasma IgA", "Plasma IgG"   ,"Smooth muscle cells" ,  "T" )]

corr_clusters_all_genes <- cor(counts_per_cluster_z, method = 'spearman')

ann <- data.frame(cluster = colnames(counts_per_cluster_z))
rownames(ann) <- colnames(counts_per_cluster_z)

ann_colors = list(cluster=c('B' = '#0F2080', 'B memory/reg' = '#228833', 'Endothelial' = '#AA3377', 'Enteric glia cells' = '#BBBBBB', 'Epithelial' = '#AEC7E8',
                            'Fibroblast' = '#CCBB44','Macrophage' = '#F7B6D2', 'Mast cells' = '#0077BB', 
                            'Plasma IgA' = '#EE7733', 'Plasma IgG'='#33BBEE','Smooth muscle cells' = '#CC3311','T' = '#7F7F7F'))

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/correlation_matrix_all_genes_named_clusters.pdf', height = 10, width = 13)
pheatmap(corr_clusters_all_genes, cluster_cols=TRUE, cluster_rows=TRUE, show_rownames = TRUE,show_colnames = TRUE, annotation_row = ann, annotation_col = ann,
         annotation_colors = ann_colors)
dev.off()


#### bubble chart

rowclus <- hclust(dist( counts_per_cluster_markers_z ))    #cluster the rows
colclus <- hclust(dist( t(counts_per_cluster_markers_z) )) #cluster the columns

hm <- hmReady(counts_per_cluster_markers_z, colclus=colclus,rowclus=rowclus)

hm$y <- hm$rowid
hm$y  <- factor(hm$y, levels = markers)
hm$x <- hm$variable
hm$x  <- factor(hm$x, levels = rev(c("T","B"   ,"B memory/reg" ,"Plasma IgA", "Plasma IgG","Macrophage" ,"Mast cells","Smooth muscle cells" ,"Endothelial" ,"Epithelial","Fibroblast"  , 
                                     "Enteric glia cells")))

my.lines<-data.frame(y=c(5.5, 7.5, 8.5, 11.5,14.5,17.5,21.5,24.5,25.5,27.5), 
                     x=rep(0,10), 
                     yend=c(5.5, 7.5, 8.5, 11.5,14.5,17.5,21.5,24.5,25.5,27.5), 
                     xend=rep(13.5,10))

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/heatmap.pdf', height = 5, width = 15)
ggplot() + 
  geom_tile(data=hm, aes(x=x, y=y, fill=value)) +
  scale_fill_gradient2(low = "blue",  mid = 'white',high = "red") +
  coord_flip() +
  theme_classic() +
  theme(axis.text.y = element_text(size = 15),axis.line=element_blank(), axis.text.x = element_text(angle=90, size = 15), axis.title = element_blank()) +
  geom_segment(data=my.lines, aes(x,y,xend=xend, yend=yend), linewidth=0.5, inherit.aes=F)

dev.off()

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/bubble.pdf', height = 6, width = 15)
ggplot() + 
  geom_point(data=hm, aes(x=x, y=y, fill=value, size=value), alpha=0.75, shape = 21) +
  scale_fill_gradient2(low = "blue",  mid = 'white',high = "red") +
  scale_size_continuous(range=c(1,6),limit=c(min(counts_per_cluster_markers_z), max(counts_per_cluster_markers_z)), breaks = seq(min(counts_per_cluster_markers_z),max(counts_per_cluster_markers_z), length.out=5)) +
  #geom_dendro(colclus,ylim=c(89,95)) +
  coord_flip() +
  theme_classic() +
  theme(axis.text.y = element_text(size = 15), axis.line=element_blank(), axis.text.x = element_text(angle=90, size = 15), axis.title = element_blank()) +
  geom_segment(data=my.lines, aes(x,y,xend=xend, yend=yend), linewidth=0.5, inherit.aes=F)

dev.off()

### stacked barplot
metadata <- read.csv('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/metadata_adata.csv', header = TRUE)

patients_low <- unique(metadata[metadata$level == 'LOW', 'Patients'])
patients_high <- unique(metadata[metadata$level == 'HIGH', 'Patients'])

df_high  <- data.frame(cell_type = rep(unique(metadata$new_leiden), each=length(patients_high)), patient = rep(patients_high, length(unique(metadata$new_leiden))), counts = 0)
lapply(as.list(patients_high), function(x) {
  counts <- table(metadata[metadata$Patients ==x, 'new_leiden'])
  for (n in names(counts)) {
    df_high[df_high$patient == x & df_high$cell_type == n, 3] <<- counts[n]
  }
  
})

df_low  <- data.frame(cell_type = rep(unique(metadata$new_leiden), each=length(patients_low)), patient = rep(patients_low, length(unique(metadata$new_leiden))), counts = 0)
lapply(as.list(patients_low), function(x) {
  counts <- table(metadata[metadata$Patients ==x, 'new_leiden'])
  for (n in names(counts)) {
    df_low[df_low$patient == x & df_low$cell_type == n, 3] <<- counts[n]
  }
  
})

df_tot <- rbind(cbind(df_high, level='HIGH'), cbind(df_low, level = 'LOW'))

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/stacked_barplot_high_low_patients.pdf', height = 8, width = 10)
ggplot(df_tot, aes(fill=cell_type, y=counts, x=patient)) + 
  geom_bar(position=position_fill(reverse = TRUE), stat="identity") +
  scale_fill_manual(values = c('B' = '#0F2080', 'B reg' = '#228833', 'Endothelial' = '#AA3377', 'Enteric glia cells' = '#BBBBBB', 'Epithelial' = '#AEC7E8',
                               'Fibroblast' = '#CCBB44','Macrophage' = '#F7B6D2', 'Mast cells' = '#0077BB', 
                               'Plasma IgA' = '#EE7733', 'Plasma IgG'='#33BBEE','Smooth muscle cells' = '#CC3311','T' = '#7F7F7F')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=14), legend.text=element_text(size=12), axis.title.x = element_blank(), 
        legend.title = element_blank(), legend.position = 'bottom') +
  facet_wrap(~level, nrow=1, ncol=2, scales = 'free', ) +
  theme(aspect.ratio=3) 
dev.off()

###### t-test cell type abundance high vs low

stats_abundance_cell_type <- df_tot %>% group_by(cell_type) %>% t_test(counts ~ level, var.equal = FALSE)
#####

df_low_mean  <- data.frame(cell_type = unique(metadata$new_leiden), counts = 0)
rownames(df_low_mean) <- df_low_mean[,1]
df_low_mean[,1] <- NULL

for (c in unique(df_low$cell_type)) {
  df_low_mean[c,] <- mean(df_low[df_low$cell_type == c,3])
}

df_low_mean$cell_type <- rownames(df_low_mean)
rownames(df_low_mean) <- 1:nrow(df_low_mean)

df_high_mean  <- data.frame(cell_type = unique(metadata$new_leiden), counts = 0)
rownames(df_high_mean) <- df_high_mean[,1]
df_high_mean[,1] <- NULL

for (c in unique(df_high$cell_type)) {
  df_high_mean[c,] <- mean(df_high[df_high$cell_type == c,3])
}

df_high_mean$cell_type <- rownames(df_high_meaon)
rownames(df_high_mean) <- 1:nrow(df_high_mean)

df_tot <- rbind(cbind(df_high_mean, level='HIGH'), cbind(df_low_mean, level = 'LOW'))

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/stacked_barplot_high_low.pdf', width = 10)
ggplot(df_tot, aes(fill=cell_type, y=counts, x=level)) + 
  geom_bar(position=position_fill(reverse = TRUE), stat="identity") +
  scale_fill_manual(values = c('B' = '#0F2080', 'B reg' = '#228833', 'Endothelial' = '#AA3377', 'Enteric glia cells' = '#BBBBBB', 'Epithelial' = '#AEC7E8',
                               'Fibroblast' = '#CCBB44','Macrophage' = '#F7B6D2', 'Mast cells' = '#0077BB', 
                               'Plasma IgA' = '#EE7733', 'Plasma IgG'='#33BBEE','Smooth muscle cells' = '#CC3311','T' = '#7F7F7F')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=14), legend.text=element_text(size=12), axis.title.x = element_blank(), 
        legend.title = element_blank(), legend.position = 'bottom') +
  theme(aspect.ratio=3) 
dev.off()

### stacked barplot by fov
metadata <- read.csv('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/metadata_adata.csv', header = TRUE)

df_fov  <- data.frame(cell_type = rep(unique(metadata$new_leiden), each=length(unique(metadata$fov))), fov = rep(unique(metadata$fov), length(unique(metadata$new_leiden))), counts = 0)
lapply(as.list(unique(metadata$fov)), function(x) {
  counts <- table(metadata[metadata$fov ==x, 'new_leiden'])
  for (n in names(counts)) {
    df_fov[df_fov$fov == x & df_fov$cell_type == n, 3] <<- counts[n]
  }
})

df_fov$level <- unlist(lapply(as.list(df_fov$fov), function(x) {unique(metadata[metadata$fov == x, 'level'])}))

df_fov$fov <- factor(df_fov$fov)
df_fov$cell_type <- factor(df_fov$cell_type)

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/stacked_barplot_fov.pdf', height = 6, width = 8)
ggplot(df_fov, aes(fill=cell_type, y=counts, x=fov)) + 
  geom_bar(position=position_fill(reverse = TRUE), stat="identity") +
  scale_fill_manual(values = c('B' = '#0F2080', 'B reg' = '#228833', 'Endothelial' = '#AA3377', 'Enteric glia cells' = '#BBBBBB', 'Epithelial' = '#AEC7E8',
                               'Fibroblast' = '#CCBB44','Macrophage' = '#F7B6D2', 'Mast cells' = '#0077BB', 
                               'Plasma IgA' = '#EE7733', 'Plasma IgG'='#33BBEE','Smooth muscle cells' = '#CC3311','T' = '#7F7F7F')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), legend.text=element_text(size=12), axis.title.x = element_blank(), 
        legend.title = element_blank(), legend.position = 'bottom') +
  facet_wrap(~level, nrow=1, ncol=2, scales = 'free', ) +
  theme(aspect.ratio=3) 
dev.off()



## stacked barplot + pca per patient

cells_per_patient <- metadata %>% group_by(Patients) %>%  dplyr::count(new_leiden)
# PCA
m <- matrix(nrow=length(unique(metadata$Patients)), ncol=length(unique(metadata$new_leiden)))
rownames(m) <- unique(metadata$Patients)
colnames(m) <- unique(metadata$new_leiden)

for (r in rownames(m)) {
  for (c in colnames(m)) {
    if (c %in% metadata[metadata$Patients == r,'new_leiden']) {
      m[r,c] <- cells_per_patient[(cells_per_patient$Patients == r) & (cells_per_patient$new_leiden == c), 'n']$n
    } else {
      m[r,c] <- 0
    }
    
  }
}

status <- unlist(lapply(as.list(rownames(m)), function(x) unique(metadata[metadata$Patients == x, 'level'])))

theme_set(theme_bw())
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             strip.background = element_blank())

dist<-vegdist(m, method="bray")
variance_exp <- ape::pcoa(dist)$values$Relative_eig

set.seed(1)
permANOVA<-adonis2(dist ~ status, permutations = 999)
# p = 0.549

fit <- cmdscale(dist,eig=TRUE, k=2)

fit_p <- hclust(dist, method="average")
pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/hierarchical_tree_patients_number_cells.pdf', width = 5, height = 5)
plot(fit_p, hang = -1, cex = 1, main = "Bray-Curtis Dissimilarity", 
     ylab = "Bray-Curtis Dissimilarity", xlab = "ID")
dev.off()

# prepare the data for plotting
mds <- data.frame(
  x    = fit$points[,1],
  y    = fit$points[,2],
  name = rownames(m),
  status = status)

# plot
pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial//spatial/CD/results_seed/PCA_patients_number_cells.pdf', width = 5, height = 5)
ggplot(mds, aes(x=x,y=y,label=name)) + 
  geom_point(aes(colour=status), alpha = 0.9, size=4) +
  xlab('PC1 [42.9%]') +
  ylab('PC2 [24.3%]') +
  scale_colour_manual(name = 'iNKTIL10 level', values = c( "#ffb703", "#f94144")) + theme(aspect.ratio=1) 
dev.off()

cells_per_patient$Patients <- paste(cells_per_patient$Patients, '-',
                                    unlist(lapply(as.list(cells_per_patient$Patients), function(x) unique(metadata[metadata$Patients == x, 'level']))))

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/stacked_barplot_num_cells_per_patient.pdf')
ggplot(cells_per_patient, aes(fill=new_leiden, y=n, x=Patients)) + 
  geom_bar(position=position_fill(reverse = TRUE), stat="identity") +
  theme_classic() +
  scale_fill_manual(values = c('B' = '#0F2080', 'B reg' = '#228833', 'Endothelial' = '#AA3377', 'Enteric glia cells' = '#BBBBBB', 'Epithelial' = '#AEC7E8',
                               'Fibroblast' = '#CCBB44','Macrophage' = '#F7B6D2', 'Mast cells' = '#0077BB', 
                               'Plasma IgA' = '#EE7733', 'Plasma IgG'='#33BBEE','Smooth muscle cells' = '#CC3311','T' = '#7F7F7F')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=14), legend.text=element_text(size=12), axis.title.x = element_blank(),axis.title.y = element_blank(), 
        legend.title = element_blank()) +
  theme(aspect.ratio = 2.2) +
  ggtitle('Number of cells per cell type in each patient')
dev.off()

### volcano plot for DEGs between epithelial/T cells within/further than 500 µm from iNKT10 cell in fov 38
df_epithelial_fov38 <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_epithelial_within_further_500_iNKT10_fov38.xlsx'))
df_epithelial_fov38[,1] <- NULL
df_epithelial_fov38$diff <- 'NO'
df_epithelial_fov38[df_epithelial_fov38$logfoldchanges>0 & df_epithelial_fov38$pvals<=0.05,]$diff <- 'within_500'
df_epithelial_fov38[df_epithelial_fov38$logfoldchanges<0 & df_epithelial_fov38$pvals<=0.05,]$diff <- 'farther_500'

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_epithelial_within_further_500_iNKT10_fov38.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_epithelial_fov38, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c(within_500="#28658f", farther_500="#ef857c")) +
  geom_point(data = subset(df_epithelial_fov38, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals<0.05,df_epithelial_fov38$names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs of epithelial cells within or farther than 500 µm from iNKT10 cell in fov 38')

dev.off()



df_T_fov38 <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_T_within_further_500_iNKT10_fov38.xlsx'))
df_T_fov38[,1] <- NULL
df_T_fov38$diff <- 'NO'
df_T_fov38[df_T_fov38$logfoldchanges>0 & df_T_fov38$pvals<=0.1,]$diff <- 'within_500'
df_T_fov38[df_T_fov38$logfoldchanges<0 & df_T_fov38$pvals<=0.1,]$diff <- 'farther_500'

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_T_within_further_500_iNKT10_fov38.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_T_fov38, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c(within_500="#28658f", farther_500="#ef857c")) +
  geom_point(data = subset(df_T_fov38, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals<0.1,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.1), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs of T cells within or farther than 500 µm from iNKT10 cell in fov 38')

dev.off()

#####
from_david_to_cyrcle <- function(file_BP, degs_high, cell_type) {
  david_BP <- read.table(paste0('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/', file_BP, '.txt'), header=TRUE, sep='\t')
  
  david_BP$up_in_high <- 0
  david_BP$up_in_low <- 0
  
  for (r in 1:nrow(david_BP)) {
    high <- c()
    low <- c()
    for (x in unlist(strsplit(david_BP[r,'Genes'], ', '))) {
      if (x %in% degs_high) {
        high <- c(high, x)
      } else {
        low <- c(low, x)
      }
    } 
    david_BP[r,]$up_in_high <- length(high)
    david_BP[r,]$up_in_low <- length(low)
  }
  
  david_BP2 <- data.frame(ID = david_BP$Term, classification = david_BP$Category, FDR = david_BP$FDR, 
                          all = david_BP$Pop.Hits, enriched_in_HIGH = david_BP$up_in_high, enriched_in_LOW = david_BP$up_in_low)
  david_BP2 <- david_BP2[order(david_BP2$FDR),]
  
  if (nrow(david_BP2) > 10) {
    david_BP2 <- david_BP2[1:10,]
  }
  write.table(david_BP2[,1], file=paste0('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/david_BP_', cell_type, '_modified_entire_GO.tsv'),col.names =TRUE, row.names =FALSE, quote = FALSE, sep='\t')
  
  david_BP2$ID <- gsub('~.*$', '', david_BP2$ID)
  
  write.table(david_BP2, file=paste0('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/david_BP_', cell_type, '_modified.tsv'),col.names =TRUE, row.names =FALSE, quote = FALSE, sep='\t')
  
  
### volcano plot for DEGs between cells of the same type, high vs low iNKT10 level
df_B <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_B_low_high.xlsx'))
df_B[,1] <- NULL
df_B$diff <- 'NO'
df_B[(df_B$logfoldchanges> 0 & df_B$pvals<=0.05),]$diff <- 'LOW'
df_B[(df_B$logfoldchanges< 0 & df_B$pvals<=0.05),]$diff <- 'HIGH'

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_B_low_high.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_B, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c(LOW="#28658f", HIGH="#ef857c")) +
  geom_point(data = subset(df_B, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals<0.0001,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs of B cells between high and low iNKT10 level')

dev.off()

write.table(df_B[df_B$diff != 'NO','names'], file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_B_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

# do' a david in input i geni con log2FC >/< 0 e pvalue < 0.05

df_B_fdr_0.05_high <- df_B[df_B$diff == "HIGH",'names']

df_B_fdr_0.05_high[df_B_fdr_0.05_high=='HLA.DRA'] <- 'HLA-DRA'
df_B_fdr_0.05_high[df_B_fdr_0.05_high=='HLA.DPB1'] <- 'HLA-DPB1'
df_B_fdr_0.05_high <- df_B_fdr_0.05_high[df_B_fdr_0.05_high != "HLA.DQB1.2"]
df_B_fdr_0.05_high <- c(df_B_fdr_0.05_high, 'HLA-DQB1','HLA-DQB2')
df_B_fdr_0.05_high[df_B_fdr_0.05_high=='HLA.DQA1'] <- 'HLA-DQA1'
df_B_fdr_0.05_high[df_B_fdr_0.05_high=='HLA.DPA1'] <- 'HLA-DPA1'

from_david_to_cyrcle('david_B_BP', df_B_fdr_0.05_high, 'B')
# prendi i geni arricchiti in ogni GO term tornato da David (quindi geni con log2FC >/< 0 e pvalue < 0.05). 
# per ogni gene se compare tra quelli arricchiti in High, lo consideri high, se no low


### volcano plot for DEGs between cells of the same type, high vs low iNKT10 level
df_Breg <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_B_reg_low_high.xlsx'))
df_Breg[,1] <- NULL
df_Breg$diff <- 'NO'
df_Breg[(df_Breg$logfoldchanges> 0 & df_Breg$pvals<=0.05),]$diff <- 'LOW'
df_Breg[(df_Breg$logfoldchanges< 0 & df_Breg$pvals<=0.05),]$diff <- 'HIGH'

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_Breg_low_high.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_Breg, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c(LOW="#28658f", HIGH="#ef857c")) +
  geom_point(data = subset(df_Breg, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals<0.0001,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs of B reg cells between high and low iNKT10 level')

dev.off()

write.table(df_Breg[df_Breg$diff != 'NO','names'], file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_Breg_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

df_Breg_fdr_0.05_high <- df_Breg[df_Breg$diff == "HIGH",'names']

df_Breg_fdr_0.05_high[df_Breg_fdr_0.05_high=='HLA.DPB1'] <- 'HLA-DPB1'

from_david_to_cyrcle('david_Breg_BP', df_Breg_fdr_0.05_high, 'Breg')



### volcano plot for DEGs between cells of the same type, high vs low iNKT10 level
df_endo <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_endo_low_high.xlsx'))
df_endo[,1] <- NULL
df_endo$diff <- 'NO'
df_endo[(df_endo$logfoldchanges> 0 & df_endo$pvals<=0.05),]$diff <- 'LOW'
df_endo[(df_endo$logfoldchanges< 0 & df_endo$pvals<=0.05),]$diff <- 'HIGH'

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_endo_low_high.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_endo, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c(LOW="#28658f", HIGH="#ef857c")) +
  geom_point(data = subset(df_endo, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals<0.0001,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs of endothelial cells between high and low iNKT10 level')

dev.off()

write.table(df_endo[df_endo$diff != 'NO','names'], file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_endo_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

df_endo_fdr_0.05_high <- df_endo[df_endo$diff == "HIGH",'names']

df_endo_fdr_0.05_high[df_endo_fdr_0.05_high=='MAP1LC3B.2'] <- 'MAP1LC3B'

from_david_to_cyrcle('david_endo_BP', df_endo_fdr_0.05_high, 'endo')


### volcano plot for DEGs between cells of the same type, high vs low iNKT10 level
df_glia <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_glia_low_high.xlsx'))
df_glia[,1] <- NULL
df_glia$diff <- 'NO'
df_glia[(df_glia$logfoldchanges> 0 & df_glia$pvals<=0.05),]$diff <- 'LOW'
df_glia[(df_glia$logfoldchanges< 0 & df_glia$pvals<=0.05),]$diff <- 'HIGH'

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_glia_low_high.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_glia, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c(LOW="#28658f", HIGH="#ef857c")) +
  geom_point(data = subset(df_glia, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals<0.01,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs of enteric glia cells between high and low iNKT10 level')

dev.off()

write.table(df_glia[df_glia$diff != 'NO','names'], file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_glia_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

df_glia_fdr_0.05_high <- df_glia[df_glia$diff == "HIGH",'names']

from_david_to_cyrcle('david_glia_BP', df_glia_fdr_0.05_high, 'glia')


### volcano plot for DEGs between cells of the same type, high vs low iNKT10 level
df_epi <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_epi_low_high.xlsx'))
df_epi[,1] <- NULL
df_epi$diff <- 'NO'
df_epi[(df_epi$logfoldchanges> 0.5 & df_epi$pvals<=0.05),]$diff <- 'LOW'
df_epi[(df_epi$logfoldchanges< -0.5 & df_epi$pvals<=0.05),]$diff <- 'HIGH'
write.table(df_epi, file = '/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/epitelial_DEGs_high_low.txt', quote = FALSE, sep='\t',row.names = F)

gene_list <- c('HSPA1A.B',	'CXCL1.2.3', 'S100P',	'PLAC8',	'HSPB1',	'LA.DQB1.2',	'S100A6',	'TPT1',	'TAGLN',  'LYZ',	'KRT20',	'COL1A1',	'CD164',
               'VIM',	'COL6A2', 'LCN2')

gene_list <- c('LCN2')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_epi_low_high.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_epi, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c(LOW="#28658f", HIGH="#ef857c")) +
  geom_point(data = subset(df_epi, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(names %in% gene_list,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) +
  geom_vline(xintercept=c(-0.5,0.5), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs of epithelial cells between high and low iNKT10 level')

dev.off()

write.table(df_epi[df_epi$diff != 'NO','names'], file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_epi_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

df_epi_fdr_0.05_high <- df_epi[df_epi$diff == "HIGH",'names']

from_david_to_cyrcle('david_epi_BP', df_epi_fdr_0.05_high, 'epi')

### volcano plot for DEGs between cells of the same type, high vs low iNKT10 level
df_fibro <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_fibro_low_high.xlsx'))
df_fibro[,1] <- NULL
df_fibro$diff <- 'NO'
df_fibro[(df_fibro$logfoldchanges> 0 & df_fibro$pvals<=0.05),]$diff <- 'LOW'
df_fibro[(df_fibro$logfoldchanges< 0 & df_fibro$pvals<=0.05),]$diff <- 'HIGH'
write.table(df_fibro, file = '/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/fibroblast_DEGs_high_low.txt', quote = FALSE, sep='\t',row.names = F)

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_fibro_low_high.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_fibro, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c(LOW="#28658f", HIGH="#ef857c")) +
  geom_point(data = subset(df_fibro, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals<0.0001,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs of fibroblasts between high and low iNKT10 level')

dev.off()

write.table(df_fibro[df_fibro$diff != 'NO','names'], file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_fibro_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

df_fibro_fdr_0.05_high <- df_fibro[df_fibro$diff == "HIGH",'names']

from_david_to_cyrcle('david_fibro_BP', df_fibro_fdr_0.05_high, 'fibro')


### volcano plot for DEGs between cells of the same type, high vs low iNKT10 level
df_macrophage <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_macro_low_high.xlsx'))
df_macrophage[,1] <- NULL
df_macrophage$diff <- 'NO'
df_macrophage[(df_macrophage$logfoldchanges> 0 & df_macrophage$pvals<=0.05),]$diff <- 'LOW'
df_macrophage[(df_macrophage$logfoldchanges< 0 & df_macrophage$pvals<=0.05),]$diff <- 'HIGH'

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_macro_low_high.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_macrophage, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c(LOW="#28658f", HIGH="#ef857c")) +
  geom_point(data = subset(df_macrophage, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals<0.0001,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs of macrophages between high and low iNKT10 level')

dev.off()

write.table(df_macrophage[df_macrophage$diff != 'NO','names'], file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_macro_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

df_macro_fdr_0.05_high <- df_macrophage[df_macrophage$diff == "HIGH",'names']
df_macro_fdr_0.05_high[df_macro_fdr_0.05_high=='HLA.DPB1'] <- 'HLA-DPB1'
df_macro_fdr_0.05_high <- df_macro_fdr_0.05_high[df_macro_fdr_0.05_high !='HLA.DRB'] 

from_david_to_cyrcle('david_macro_BP', df_macro_fdr_0.05_high, 'macro')


### volcano plot for DEGs between cells of the same type, high vs low iNKT10 level
df_mast <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_mast_low_high.xlsx'))
df_mast[,1] <- NULL
df_mast$diff <- 'NO'
df_mast[(df_mast$logfoldchanges> 0 & df_mast$pvals<=0.05),]$diff <- 'LOW'
df_mast[(df_mast$logfoldchanges< 0 & df_mast$pvals<=0.05),]$diff <- 'HIGH'

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_mast_low_high.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_mast, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c(LOW="#28658f", HIGH="#ef857c")) +
  geom_point(data = subset(df_mast, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals<0.0001,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs of mast cells between high and low iNKT10 level')

dev.off()

write.table(df_mast[df_mast$diff != 'NO','names'], file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_mast_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

df_mast_fdr_0.05_high <- df_mast[df_mast$diff == "HIGH",'names']

from_david_to_cyrcle('david_mast_BP', df_mast_fdr_0.05_high, 'mast')


### volcano plot for DEGs between cells of the same type, high vs low iNKT10 level
df_plasma_A <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_plasma_A_low_high.xlsx'))
df_plasma_A[,1] <- NULL
df_plasma_A$diff <- 'NO'
df_plasma_A[(df_plasma_A$logfoldchanges> 0 & df_plasma_A$pvals<=0.05),]$diff <- 'LOW'
df_plasma_A[(df_plasma_A$logfoldchanges< 0 & df_plasma_A$pvals<=0.05),]$diff <- 'HIGH'

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_plasma_A_low_high.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_plasma_A, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c(LOW="#28658f", HIGH="#ef857c")) +
  geom_point(data = subset(df_plasma_A, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals<0.00001,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs of plasma IgA cells between high and low iNKT10 level')

dev.off()

write.table(df_plasma_A[df_plasma_A$diff != 'NO','names'], file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_plasma_A_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

df_plasma_A_fdr_0.05_high <- df_plasma_A[df_plasma_A$diff == "HIGH",'names']

from_david_to_cyrcle('david_plasma_IgA_BP', df_plasma_A_fdr_0.05_high, 'Plasma IgA')


### volcano plot for DEGs between cells of the same type, high vs low iNKT10 level
df_plasma_G <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_plasma_G_low_high.xlsx'))
df_plasma_G[,1] <- NULL
df_plasma_G$diff <- 'NO'
df_plasma_G[(df_plasma_G$logfoldchanges> 0 & df_plasma_G$pvals<=0.05),]$diff <- 'LOW'
df_plasma_G[(df_plasma_G$logfoldchanges< 0 & df_plasma_G$pvals<=0.05),]$diff <- 'HIGH'

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_plasma_G_low_high.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_plasma_G, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c(LOW="#28658f", HIGH="#ef857c")) +
  geom_point(data = subset(df_plasma_G, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals<0.0001,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs of plasma IgG cells between high and low iNKT10 level')

dev.off()

write.table(df_plasma_G[df_plasma_G$diff != 'NO','names'], file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_plasma_G_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

df_plasma_G_fdr_0.05_high <- df_plasma_G[df_plasma_G$diff == "HIGH",'names']

from_david_to_cyrcle('david_plasma_IgG_BP', df_plasma_G_fdr_0.05_high, 'Plasma IgG')


### volcano plot for DEGs between cells of the same type, high vs low iNKT10 level
df_muscle <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_muscle_low_high.xlsx'))
df_muscle[,1] <- NULL
df_muscle$diff <- 'NO'
df_muscle[(df_muscle$logfoldchanges> 0 & df_muscle$pvals<=0.05),]$diff <- 'LOW'
df_muscle[(df_muscle$logfoldchanges< 0 & df_muscle$pvals<=0.05),]$diff <- 'HIGH'

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_muscle_low_high.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_muscle, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c(LOW="#28658f", HIGH="#ef857c")) +
  geom_point(data = subset(df_muscle, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals<0.0001,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs of smooth muscle cells between high and low iNKT10 level')

dev.off()

write.table(df_muscle[df_muscle$diff != 'NO','names'], file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_muscle_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

df_muscle_fdr_0.05_high <- df_muscle[df_muscle$diff == "HIGH",'names']

from_david_to_cyrcle('david_muscle_BP', df_muscle_fdr_0.05_high, 'muscle')


### volcano plot for DEGs between cells of the same type, high vs low iNKT10 level
df_T <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_T_low_high.xlsx'))
df_T[,1] <- NULL
df_T$diff <- 'NO'
df_T[(df_T$logfoldchanges> 0.5 & df_T$pvals<=0.05),]$diff <- 'LOW'
df_T[(df_T$logfoldchanges< -0.5 & df_T$pvals<=0.05),]$diff <- 'HIGH'
write.table(df_T, file = '/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/T_DEGs_high_low.txt', quote = FALSE, sep='\t',row.names = F)

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_T_low_high.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_T, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c(LOW="#28658f", HIGH="#ef857c")) +
  geom_point(data = subset(df_T, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals<0.0001 & abs(df_T$logfoldchanges)> 0.5,names,'')), size =4, max.overlaps = 1000)  +
  geom_vline(xintercept=c(-0.5,0.5), col="red", linewidth=0.3) + 
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs of T cells between high and low iNKT10 level')

dev.off()

write.table(df_T[df_T$diff != 'NO','names'], file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_T_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

df_T_fdr_0.05_high <- df_T[df_T$diff == "HIGH",'names']

from_david_to_cyrcle('david_T_BP', df_T_fdr_0.05_high, 'T')

### distances between epithelial and each other cell type considering the means computed for each fov
all_distances_epi_mean <- matrix(ncol=3)
colnames(all_distances_epi_mean) <- c('iNKT10_level', 'Distance','Cell_type')
for (f in list.files('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_epi_mean', full.names = TRUE)) {
  dist <- data.frame(readxl::read_excel(f))
  colnames(dist) <- c('iNKT10_level', 'Distance')
  dist$Cell_type <- gsub('.xlsx', '',strsplit(f, split='/')[[1]][12]) 
  all_distances_epi_mean <- rbind(all_distances_epi_mean, dist)
}
all_distances_epi_mean <-  all_distances_epi_mean[-1,] 

all_distances_epi_mean2 <- all_distances_epi_mean[-which(all_distances_epi_mean$Cell_type %in% c('B')),]

stats_epi <- all_distances_epi_mean2 %>%
  group_by(Cell_type) %>%
  t_test(Distance~iNKT10_level,) %>%
  adjust_pvalue(method = 'fdr')

stats_epi$y.position <- unlist(lapply(seq_along(1:length(stats_epi$Cell_type)), function(x) {
  max(all_distances_epi_mean2[all_distances_epi_mean2$Cell_type == stats_epi$Cell_type[x],'Distance'][! is.na(all_distances_epi_mean2[all_distances_epi_mean2$Cell_type == stats_epi$Cell_type[x],'Distance'])]) +max(all_distances_epi_mean2[all_distances_epi_mean2$Cell_type == stats_epi$Cell_type[x],'Distance'][! is.na(all_distances_epi_mean2[all_distances_epi_mean2$Cell_type == stats_epi$Cell_type[x],'Distance'])])/14
}))

stats_epi$p <- round(stats_epi$p, 3)

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_epi_distances_by_fov_points.pdf', width = 30, height = 20)
ggplot(all_distances_epi_mean, aes(x= iNKT10_level, y=Distance,color=iNKT10_level)) +
  geom_boxplot(aes(colour=iNKT10_level),width=0.4, lwd=1, outliers = FALSE) +
  geom_point(aes(colour=iNKT10_level), position=position_jitter(seed=1,width = .3),alpha = 1, size=3) +
  scale_colour_manual(values = c("#ef857c","#28658f")) +
  facet_wrap(~Cell_type,nrow = 3, ncol = 5,scale = "free_y") +
  theme_classic() +
  scale_shape(solid=T)  +
  ylab("") +
  xlab("") +
  theme(aspect.ratio=1.5, strip.text=element_text(size = 25, face = 'bold'), axis.text.x = element_text(size=25, angle = 45,hjust = 1), 
        axis.text.y = element_text(size=25),
        legend.text = element_text(size = 25),legend.title =  element_text(size = 25), strip.text.x = element_text(size=20)) +
  stat_pvalue_manual(stats_epi, label='p', size=7)
dev.off()


## bubble plot with mean distance between the mean distances between epi and each other cell type in the different fovs
distances_low <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_low_mean.txt', header=TRUE, sep=',')
distances_high <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_high_mean.txt', header=TRUE, sep=',')

distances_high_low <- cbind(distances_low, mean_high = NA)

for (r in 1:nrow(distances_high)) {
  if (length(distances_high_low[(distances_high_low$cell_type1 == distances_high[r,1]) & (distances_high_low$cell_type2 == distances_high[r,2]),'mean']) != 0) {
    distances_high_low[(distances_high_low$cell_type1 == distances_high[r,1]) & (distances_high_low$cell_type2 == distances_high[r,2]), 'mean_high'] <- distances_high[r, 'mean']
  } else if (length(distances_high_low[(distances_high_low$cell_type2 == distances_high[r,1]) & (distances_high_low$cell_type1 == distances_high[r,2]), 'mean']) != 0) {
    distances_high_low[(distances_high_low$cell_type2 == distances_high[r,1]) & (distances_high_low$cell_type1 == distances_high[r,2]), 'mean_high'] <- distances_high[r, 'mean']
  } else if ((length(distances_high_low[(distances_high_low$cell_type1 == distances_high[r,1]) & (distances_high_low$cell_type2 == distances_high[r,2]),'mean']) == 0) & (length(distances_high_low[(distances_high_low$cell_type2 == distances_high[r,1]) & (distances_high_low$cell_type1 == distances_high[r,2]), 'mean']) == 0)) {
    distances_high_low <- rbind(distances_high_low, c(distances_high[r,1],distances_high[r,2], 0, distances_high[r,'mean']))
  }
}

distances_high_low_epi <- distances_high_low[(distances_high_low$cell_type1 == 'Epithelial') |(distances_high_low$cell_type2 == 'Epithelial') ,]

lapply(as.list(1:nrow(distances_high_low_epi)), function(x) {if (distances_high_low_epi[x,2] != 'Epithelial') {distances_high_low_epi[x,1] <<- distances_high_low_epi[x,2]
distances_high_low_epi[x,2] <<- 'Epithelial'}})

distances_df <- data.frame(cell = rep(distances_high_low_epi$cell_type1, 2), condition = rep(c('High iNKT10', 'Low iNKT10'), each = nrow(distances_high_low_epi)),
                           value = c(distances_high_low_epi$mean_high, distances_high_low_epi$mean))

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/epithelial_high_low_mean.pdf', width = 5, height = 10)
ggplot(distances_df, aes(x=condition , y=cell, size = -value)) +
  geom_point(alpha=0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6, size = 15), axis.title = element_blank(), axis.text.y = element_text(size=15))
dev.off()

colnames(distances_high_low_epi) <- c('Cell type1', 'Cell type2', 'Mean distance in low iNKT10', 'Mean distance in high iNKT10')


### compute the % of ZBTB16+ T cells and of ZBTB16+IL10+ T cells over the total amount of cells per fov
num_cells_each_fov <- c(1037,2161,3193,829,1184,1466,1734,2719, 1662,1826, 3211,1395, 1419, 1483, 2152, 2020, 3844, 3010, 3487)
names(num_cells_each_fov) <- c('21','23','24','25','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41')

fov_low_ZBTB16 = c(1, 8, 6, 162, 12, 5, 197, 6, 120)
names(fov_low_ZBTB16) <- c('21','25','29','30','37','38', '39','40','41')
fov_low_ZBTB16_IL10 = c(0, 1, 1, 6, 1, 0, 6, 0, 3)
names(fov_low_ZBTB16_IL10) <- c('21','25','29','30','37','38', '39','40','41')

fov_high_ZBTB16 = c(9, 42, 9, 5, 10, 9, 272, 14, 4, 16)
names(fov_high_ZBTB16) <- c('23','24','27','28','31','32','33','34','35','36')
fov_high_ZBTB16_IL10 = c(1, 0, 1,0,0,0,4,0,0,1)
names(fov_high_ZBTB16_IL10) <- c('23','24','27','28','31','32','33','34','35','36')

fov_low_ZBTB16_perc <- c()
for (i in 1:length(fov_low_ZBTB16)) {
  fov_low_ZBTB16_perc <- c(fov_low_ZBTB16_perc, fov_low_ZBTB16[i]*100/num_cells_each_fov[names(fov_low_ZBTB16)[i]])
}

fov_low_ZBTB16_IL10_perc <- c()
for (i in 1:length(fov_low_ZBTB16_IL10)) {
  fov_low_ZBTB16_IL10_perc <- c(fov_low_ZBTB16_IL10_perc, fov_low_ZBTB16_IL10[i]*100/num_cells_each_fov[names(fov_low_ZBTB16_IL10)[i]])
}

fov_high_ZBTB16_perc <- c()
for (i in 1:length(fov_high_ZBTB16)) {
  fov_high_ZBTB16_perc <- c(fov_high_ZBTB16_perc, fov_high_ZBTB16[i]*100/num_cells_each_fov[names(fov_high_ZBTB16)[i]])
}

fov_high_ZBTB16_IL10_perc <- c()
for (i in 1:length(fov_high_ZBTB16_IL10)) {
  fov_high_ZBTB16_IL10_perc <- c(fov_high_ZBTB16_IL10_perc, fov_high_ZBTB16_IL10[i]*100/num_cells_each_fov[names(fov_high_ZBTB16_IL10)[i]])
}

df = data.frame(iNKT10_level = rep(rep(c('low', 'high'), c(length(fov_low_ZBTB16_perc), length(fov_high_ZBTB16_perc))),2), 
                value = c(fov_low_ZBTB16_perc,fov_high_ZBTB16_perc,fov_low_ZBTB16_IL10_perc,fov_high_ZBTB16_IL10_perc),
                type2 = rep(c('T expressing ZBTB16', 'T expressing ZBTB16 and IL10'), c(length(c(fov_low_ZBTB16_perc,fov_high_ZBTB16_perc)), 
                                                                                        length(c(fov_low_ZBTB16_IL10_perc,fov_high_ZBTB16_IL10_perc)))), 
                type = rep(c('Low iNKT10 - T expressing ZBTB16', 'High iNKT10 - T expressing ZBTB16',
                             'Low iNKT10 - T expressing ZBTB16 and IL10', 'High iNKT10 - T expressing ZBTB16 and IL10'), 
                           c(length(fov_low_ZBTB16_perc), length(fov_high_ZBTB16_perc), length(fov_low_ZBTB16_IL10_perc), length(fov_high_ZBTB16_IL10_perc))))

stats <- df %>%
  group_by(type2) %>%
  t_test(value~iNKT10_level) %>%
  adjust_pvalue(method = 'fdr')

stats$y.position <- c(10.2,1)

df$iNKT10_level <- factor(df$iNKT10_level)
df$type2 <- factor(df$type2)

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_T_ZBTB16_IL10_per_fov.pdf')
ggplot(df, aes(x= iNKT10_level, y=value)) +
  geom_boxplot(aes(colour=iNKT10_level),width=0.4, lwd=1, outliers = FALSE) +
  geom_point(aes(colour=iNKT10_level), position=position_jitter(seed=1,width = .3),alpha = 1, size=3) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  ylab("% of cells per fov") +
  xlab("") +
  stat_pvalue_manual(stats,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_wrap(~type2) 

dev.off()

### number of epithelial/fibroblasts within x micro from T ZBTB16+IL10+, IL10+, ZTBT16+

percentage_cells_within <- function(distance) {
  if (length(readBin(file(paste0('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/number_epithelial_less_', as.character(distance), '_ZBTB16_IL10'),open="rb"),"raw", 65536)) > 0) {
    epi_less_ZBTB16_IL10 <- read.table(paste0('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/number_epithelial_less_', as.character(distance), '_ZBTB16_IL10'))
  }
  if (length(readBin(file(paste0('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/number_fibroblast_less_', as.character(distance), '_ZBTB16_IL10'),open="rb"),"raw", 65536)) > 0) {
    fibro_less_ZBTB16_IL10 <- read.table(paste0('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/number_fibroblast_less_', as.character(distance), '_ZBTB16_IL10'))
  }
  epi_less_ZBTB16 <- read.table(paste0('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/number_epithelial_less_', as.character(distance),'_ZBTB16'))
  fibro_less_ZBTB16 <- read.table(paste0('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/number_fibroblast_less_', as.character(distance),'_ZBTB16'))
  epi_less_IL10 <- read.table(paste0('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/number_epithelial_less_', as.character(distance),'_IL10'))
  fibro_less_IL10 <- read.table(paste0('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/number_fibroblast_less_', as.character(distance),'_IL10'))
  
  if (exists(x='epi_less_ZBTB16_IL10')) {
    epi <- data.frame(percentage = c(epi_less_ZBTB16_IL10$V1,epi_less_ZBTB16$V1,epi_less_IL10$V1), 
                      type = c(rep('within ZBTB16+IL10+ T cells', length(epi_less_ZBTB16_IL10$V1)), rep('within ZBTB16+ T cells', length(epi_less_ZBTB16$V1)),
                               rep('within IL10+ T cells', length(epi_less_IL10$V1))))
  } else {
    epi <- data.frame(percentage = c(epi_less_ZBTB16$V1,epi_less_IL10$V1), 
                      type = c(rep('within ZBTB16+ T cells', length(epi_less_ZBTB16$V1)),
                               rep('within IL10+ T cells', length(epi_less_IL10$V1))))
  }
  
  p <- ggplot(epi, aes(x= type, y=percentage)) +
    geom_boxplot(aes(colour=type),width=0.4, lwd=1, outliers = FALSE) +
    geom_point(size=3) +
    theme_classic() +
    scale_shape(solid=T)  +
    ylab("% of cells over the total amount of epithelial cells per fov") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5), legend.position = 'none') +
    ggtitle(paste('Epithelial cells within', as.character(distance), 'µm\nfrom iNKT cells'))
  
  ggsave(paste0('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/number_epi_cells_',as.character(distance), '_T.pdf'),p, height = 5, width =3.5)
  
  if (exists('fibro_less_ZBTB16_IL10')) {
    fibro <- data.frame(percentage = c(fibro_less_ZBTB16_IL10$V1,fibro_less_ZBTB16$V1,fibro_less_IL10$V1), 
                        type = c(rep('within ZBTB16+IL10+ T cells', length(fibro_less_ZBTB16_IL10$V1)), rep('within ZBTB16+ T cells', length(fibro_less_ZBTB16$V1)),
                                 rep('within IL10+ T cells', length(fibro_less_IL10$V1))))
  } else {
    fibro <- data.frame(percentage = c(fibro_less_ZBTB16$V1,fibro_less_IL10$V1), 
                        type = c(rep('within ZBTB16+ T cells', length(fibro_less_ZBTB16$V1)),
                                 rep('within IL10+ T cells', length(fibro_less_IL10$V1))))
  }
  
  
  p <- ggplot(fibro, aes(x= type, y=percentage)) +
    geom_boxplot(aes(colour=type),width=0.4, lwd=1, outliers = FALSE) +
    geom_point(size=3) +
    theme_classic() +
    scale_shape(solid=T)  +
    ylab("% of cells over the total amount of fibroblasts per fov") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5), legend.position = 'none') +
    ggtitle(paste('Fibroblast cells within', as.character(distance), 'µm\nfrom iNKT cells'))
  
  ggsave(paste0('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/number_fibro_cells_', as.character(distance), '_T.pdf'),p, height = 5, width =3.5)
  
}

percentage_cells_within(100)
percentage_cells_within(500)

### DEGs between epithelial cells within 500 micro from T ZBTB16+IL10 and epithelial cells more far than 500 micro from T ZBTB16+IL10 
df_epi <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_epithelial_within_farther_500_ZBTB16_IL10.xlsx'))
df_epi[,1] <- NULL
df_epi$diff <- 'NO'
df_epi[(df_epi$logfoldchanges> 0 & df_epi$pvals<=0.05),]$diff <- 'within'
df_epi[(df_epi$logfoldchanges< 0 & df_epi$pvals<=0.05),]$diff <- 'farther'

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_epithelial_less_more_500_T_ZBTB16_IL10.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_epi, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('within'="#28658f", 'farther'="#ef857c")) +
  geom_point(data = subset(df_epi, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals_adj<0.05,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs between epithelial cells closer and farther than 500 µm from T ZBTB16+IL10+ cells')

dev.off()


### DEGs between fibroblast cells within 500 micro from T ZBTB16+IL10 and fibro more far than 500 micro from T ZBTB16+IL10 
df_fibro <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_fibroblast_within_farther_500_ZBTB16_IL10.xlsx'))
df_fibro[,1] <- NULL
df_fibro$diff <- 'NO'
df_fibro[(df_fibro$logfoldchanges> 0 & df_fibro$pvals<=0.05),]$diff <- 'within'
df_fibro[(df_fibro$logfoldchanges< 0 & df_fibro$pvals<=0.05),]$diff <- 'farther'

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_fibroblast_less_more_500_T_ZBTB16_IL10.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_fibro, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('within'="#28658f", 'farther'="#ef857c")) +
  geom_point(data = subset(df_fibro, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals_adj<0.05,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs between fibroblasts closer and farther than 500 µm from T ZBTB16+IL10+ cells')

dev.off()


### DEGs between epithelial cells within 500 micro from T ZBTB16+IL10 and epithelial cells within 500 micro from T ZBTB16
df_epi <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_epithelial_within_500_ZBTB16_IL10_ZBTB16.xlsx'))
df_epi[,1] <- NULL
df_epi$diff <- 'NO'
df_epi[(df_epi$logfoldchanges> 0 & df_epi$pvals<=0.05),]$diff <- 'within 500µm from T ZBTB16+IL10+ cells'
df_epi[(df_epi$logfoldchanges< 0 & df_epi$pvals<=0.05),]$diff <- 'within 500µm from T ZBTB16+ cells'

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_epithelial_less_500_T_ZBTB16_IL10_less_500_T_ZBTB10.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_epi, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('within 500µm from T ZBTB16+IL10+ cells'="#28658f", 'within 500µm from T ZBTB16+ cells'="#ef857c")) +
  geom_point(data = subset(df_epi, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals_adj<0.05,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs between epithelial cells within 500 µm from T ZBTB16+IL10+ cells and\nepithelial cells within 500 µm from T ZBTB16+ cells')

dev.off()


### DEGs between fibroblast cells within 500 micro from T ZBTB16+IL10 and fibroblast cells within 500 micro from T ZBTB16
df_fibro <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_fibroblast_within_500_ZBTB16_IL10_ZBTB16.xlsx'))
df_fibro[,1] <- NULL
df_fibro$diff <- 'NO'
df_fibro[(df_fibro$logfoldchanges> 0 & df_fibro$pvals<=0.05),]$diff <- 'within 500µm from T ZBTB16+IL10+ cells'
df_fibro[(df_fibro$logfoldchanges< 0 & df_fibro$pvals<=0.05),]$diff <- 'within 500µm from T ZBTB16+ cells'

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_fibroblast_less_500_T_ZBTB16_IL10_less_500_T_ZBTB10.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_fibro, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('within 500µm from T ZBTB16+IL10+ cells'="#28658f", 'within 500µm from T ZBTB16+ cells'="#ef857c")) +
  geom_point(data = subset(df_fibro, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals_adj<0.05,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs between fibroblasts within 500 µm from T ZBTB16+IL10+ cells and\nfibroblasts within 500 µm from T ZBTB16+ cells')

dev.off()


### DEGs between epithelial cells within 500 micro from T ZBTB16+IL10 and epithelial cells within 500 micro from T IL10
df_epi <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_epithelial_within_500_ZBTB16_IL10_IL10.xlsx'))
df_epi[,1] <- NULL
df_epi$diff <- 'NO'
df_epi[(df_epi$logfoldchanges> 0 & df_epi$pvals<=0.05),]$diff <- 'within 500µm from T ZBTB16+IL10+ cells'
df_epi[(df_epi$logfoldchanges< 0 & df_epi$pvals<=0.05),]$diff <- 'within 500µm from T IL10+ cells'

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_epithelial_less_500_T_ZBTB16_IL10_less_500_T_IL10.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_epi, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('within 500µm from T ZBTB16+IL10+ cells'="#28658f", 'within 500µm from T IL10+ cells'="#ef857c")) +
  geom_point(data = subset(df_epi, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals_adj<0.05,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs between epithelial cells within 500 µm from T ZBTB16+IL10+ cells and\nepithelial cells within 500 µm from T IL10+ cells')

dev.off()

### DEGs between fibroblasts within 500 nm from T ZBTB16+IL10 and fibroblasts within 500 nm from T IL10
df_fibro <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_fibroblast_within_500_ZBTB16_IL10_IL10.xlsx'))
df_fibro[,1] <- NULL
df_fibro$diff <- 'NO'
df_fibro[(df_fibro$logfoldchanges> 0 & df_fibro$pvals<=0.05),]$diff <- 'within 500µm from T ZBTB16+IL10+ cells'
df_fibro[(df_fibro$logfoldchanges< 0 & df_fibro$pvals<=0.05),]$diff <- 'within 500µm from T IL10+ cells'

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_fibroblast_less_500_T_ZBTB16_IL10_less_500_T_IL10.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_fibro, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('within 500µm from T ZBTB16+IL10+ cells'="#28658f", 'within 500µm from T IL10+ cells'="#ef857c")) +
  geom_point(data = subset(df_fibro, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals_adj<0.05,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs between fibroblasts within 500 µm from T ZBTB16+IL10+ cells and\nfibroblasts within 500 µm from T IL10+ cells')

dev.off()

### ripley
ripley_low <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/ripley_per_fov_low.xlsx'))
ripley_high <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/ripley_per_fov_high.xlsx'))

rownames(ripley_low) <- ripley_low$...1
ripley_low$...1 <- NULL

rownames(ripley_high) <- ripley_high$...1
ripley_high$...1 <- NULL

ripley_low <- ripley_low[,colnames(ripley_high)]
colnames(ripley_high) <- colnames(ripley_low) <- c("Epithelial","Endothelial","Macrophage","Mast cells" ,"Fibroblast","Plasma IgA","T","Plasma IgG","B reg","B","Smooth muscle cells", "Enteric glia cells" )

ripley_stat <- data.frame(stats = c(as.vector(as.matrix(ripley_low)),as.vector(as.matrix(ripley_high))),
                          cell_type = c(rep(colnames(ripley_low), each = nrow(ripley_low)), rep(colnames(ripley_high), each = nrow(ripley_high))),
                          level = c(rep('low', length(rep(colnames(ripley_low), each = nrow(ripley_low)))), rep('high', length(rep(colnames(ripley_high), each = nrow(ripley_high))))),
                          fov = c(rep(rownames(ripley_low), ncol(ripley_low)), rep(rownames(ripley_high), ncol(ripley_high))))

ripley_stat <- ripley_stat[-which(ripley_stat$stats == '[]'),]
ripley_stat$stats <- as.numeric(ripley_stat$stats)

stats <- ripley_stat %>%
  group_by(cell_type) %>%
  wilcox_test(stats~level) %>%
  adjust_pvalue(method = 'fdr')

stats$y.position <- unlist(lapply(seq_along(1:length(stats$cell_type)), function(x) {
  max(ripley_stat[ripley_stat$cell_type == stats$cell_type[x],'stats'][! is.na(ripley_stat[ripley_stat$cell_type == stats$cell_type[x],'stats'])]) +max(ripley_stat[ripley_stat$cell_type == stats$cell_type[x],'stats'][! is.na(ripley_stat[ripley_stat$cell_type == stats$cell_type[x],'stats'])])/12
}))

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_ripley_low_high_wilcoxon.pdf', height = 7, width = 10)
ggplot(ripley_stat, aes(x= level, y=stats)) +
  geom_boxplot(aes(colour=level),width=0.4, lwd=0.5, outliers = FALSE) +
  geom_point(aes(colour=level), position=position_jitter(seed=1,width = .3),alpha = 1, size=1.5) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  xlab("") +
  stat_pvalue_manual(stats,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_wrap(~cell_type, nrow = 2, ncol = 6) +
  ylab("Ripley’s L function") 

dev.off()

#colors <- c('#332288', '#CC6677', '#88CCEE', '#999933', '#882255', '#44AA99', '#DDCC77', '#AA4499', '#648FFF', '#FE6100', '#785EF0', '#FFB000','#DC267F', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7')

ripley_stat$fov <- factor(ripley_stat$fov, levels= c(21,25,29,30,37,38,39,40,41,23,24,27,28,31,32,33,34,35,36))
pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/grouped_barplot_ripley_low_high.pdf', height = 7, width = 15)
ggplot(ripley_stat, aes(fill=fov, y=stats, x=cell_type)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c(rep('#AA3377',9), rep('#0F2080',10))) +
  theme_classic() +
  xlab('Cell type')+
  ylab('Ripley\'s statistic') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

## mean ripley score considering all the cell types
ripley_mean_by_fov <- matrix(nrow = 19, ncol = 1)
colnames(ripley_mean_by_fov) <- 'mean_Ripley'
rownames(ripley_mean_by_fov) <- c(rownames(ripley_high), rownames(ripley_low))

lapply(as.list(levels(ripley_stat$fov)), function (x) {
  y <- mean(ripley_stat[ripley_stat$fov == as.numeric(x),]$stats)
  ripley_mean_by_fov[x,1] <<- y})

ripley_mean_by_fov <- data.frame(ripley_mean_by_fov)
ripley_mean_by_fov$level <- ifelse(rownames(ripley_mean_by_fov) %in% c(21,25,29,30,37,38,39,40,41), 'low', 'high')

stats <- ripley_mean_by_fov %>%
  wilcox_test(mean_Ripley~level) 

stats$y.position <- max(ripley_mean_by_fov$mean_Ripley) +5

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_mean_ripley_low_high_wilcoxon.pdf', height = 5, width = 3)
ggplot(ripley_mean_by_fov, aes(x= level, y=mean_Ripley)) +
  geom_boxplot(aes(colour=level),width=0.4, lwd=0.5, outliers = FALSE) +
  geom_point(aes(colour=level), position=position_jitter(seed=1,width = .3),alpha = 1, size=1.5) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  xlab("") +
  stat_pvalue_manual(stats,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  ylab("Mean Ripley’s L function") 

dev.off()

## mean ripley score without considering B, Breg, mast, macrophage, T, enteric glia
ripley_mean_by_fov <- matrix(nrow = 19, ncol = 1)
colnames(ripley_mean_by_fov) <- 'mean_Ripley'
rownames(ripley_mean_by_fov) <- c(rownames(ripley_high), rownames(ripley_low))

lapply(as.list(levels(ripley_stat$fov)), function (x) {
  y <- mean(ripley_stat[(ripley_stat$fov == as.numeric(x)) & (!ripley_stat$cell_type %in% c("Enteric glia cells", 'B','B reg', 'T', 'Macrophage', 
                                                                                          'Mast cells')),]$stats)
  ripley_mean_by_fov[x,1] <<- y})

ripley_mean_by_fov <- data.frame(ripley_mean_by_fov)
ripley_mean_by_fov$level <- ifelse(rownames(ripley_mean_by_fov) %in% c(21,25,29,30,37,38,39,40,41), 'low', 'high')

stats <- ripley_mean_by_fov %>%
  wilcox_test(mean_Ripley~level) 

stats$y.position <- max(ripley_mean_by_fov$mean_Ripley) +3

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_mean_ripley_low_high_wilcoxon_without_immune.pdf', height = 5, width = 3)
ggplot(ripley_mean_by_fov, aes(x= level, y=mean_Ripley)) +
  geom_boxplot(aes(colour=level),width=0.4, lwd=0.5, outliers = FALSE) +
  geom_point(aes(colour=level), position=position_jitter(seed=1,width = .3),alpha = 1, size=1.5) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  xlab("") +
  stat_pvalue_manual(stats,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  ylab("Mean Ripley’s L function") 

dev.off()

## mean ripley score only considering fibro, muscle, epi and endo
ripley_mean_by_fov <- matrix(nrow = 19, ncol = 1)
colnames(ripley_mean_by_fov) <- 'mean_Ripley'
rownames(ripley_mean_by_fov) <- c(rownames(ripley_high), rownames(ripley_low))

lapply(as.list(levels(ripley_stat$fov)), function (x) {
  y <- mean(ripley_stat[(ripley_stat$fov == as.numeric(x)) & (ripley_stat$cell_type %in% c('Fibroblast', 'Epithelial', 'Endothelial', "Smooth muscle cells")),]$stats)
  ripley_mean_by_fov[x,1] <<- y})

ripley_mean_by_fov <- data.frame(ripley_mean_by_fov)
ripley_mean_by_fov$level <- ifelse(rownames(ripley_mean_by_fov) %in% c(21,25,29,30,37,38,39,40,41), 'low', 'high')

stats <- ripley_mean_by_fov %>%
  wilcox_test(mean_Ripley~level) 

stats$y.position <- max(ripley_mean_by_fov$mean_Ripley) +3

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_mean_ripley_low_high_wilcoxon_fibro_endo_epi_muscle.pdf', height = 5, width = 3)
ggplot(ripley_mean_by_fov, aes(x= level, y=mean_Ripley)) +
  geom_boxplot(aes(colour=level),width=0.4, lwd=0.5, outliers = FALSE) +
  geom_point(aes(colour=level), position=position_jitter(seed=1,width = .3),alpha = 1, size=1.5) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  xlab("") +
  stat_pvalue_manual(stats,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  ylab("Mean Ripley’s L function") 

dev.off()

### centrality measures
degree_low <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/degree_per_fov_low.xlsx'))
degree_high <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/degree_per_fov_high.xlsx'))

rownames(degree_low) <- degree_low$...1
rownames(degree_high) <- degree_high$...1
degree_low[,1] <- NULL
degree_high[,1] <- NULL
degree_high <- degree_high[,colnames(degree_low)]

degree_stat <- data.frame(stats = c(as.vector(as.matrix(degree_low)),as.vector(as.matrix(degree_high))),
                          cell_type = c(rep(colnames(degree_low), each = nrow(degree_low)), rep(colnames(degree_high), each = nrow(degree_high))),
                          level = c(rep('low', length(rep(colnames(degree_low), each = nrow(degree_low)))), rep('high', length(rep(colnames(degree_high), each = nrow(degree_high))))),
                          fov = c(rep(rownames(degree_low), ncol(degree_low)), rep(rownames(degree_high), ncol(degree_high))))

degree_stat$stats <- as.numeric(degree_stat$stats)

stats <- degree_stat %>%
  group_by(cell_type) %>%
  wilcox_test(stats~level) %>%
  adjust_pvalue(method = 'fdr')

stats$y.position <- unlist(lapply(seq_along(1:length(stats$cell_type)), function(x) {
  max(degree_stat[degree_stat$cell_type == stats$cell_type[x],'stats'][! is.na(degree_stat[degree_stat$cell_type == stats$cell_type[x],'stats'])]) +max(degree_stat[degree_stat$cell_type == stats$cell_type[x],'stats'][! is.na(degree_stat[degree_stat$cell_type == stats$cell_type[x],'stats'])])/12
}))

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_degree_low_high.pdf', height = 7, width = 10)
ggplot(degree_stat, aes(x= level, y=stats)) +
  geom_boxplot(aes(colour=level),width=0.4, lwd=0.5, outliers = FALSE) +
  geom_point(aes(colour=level), position=position_jitter(seed=1,width = .3),alpha = 1, size=1.5) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  xlab("") +
  stat_pvalue_manual(stats,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_wrap(~cell_type, nrow = 2, ncol = 6) +
  ylab("Degree centrality") 

dev.off()

degree_stat$fov <- factor(degree_stat$fov, levels= c(21,25,29,30,37,38,39,40,41,23,24,27,28,31,32,33,34,35,36))
pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/grouped_barplot_degree_low_high.pdf', height = 7, width = 15)
ggplot(degree_stat, aes(fill=fov, y=stats, x=cell_type)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c(rep('#AA3377',9), rep('#0F2080',10))) +
  theme_classic() +
  xlab('Cell type')+
  ylab('Degree centrality') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

## mean degree score considering all the cell types
degree_mean_by_fov <- matrix(nrow = 19, ncol = 1)
colnames(degree_mean_by_fov) <- 'mean_degree'
rownames(degree_mean_by_fov) <- c(rownames(degree_high), rownames(degree_low))

lapply(as.list(levels(degree_stat$fov)), function (x) {
  y <- mean(degree_stat[degree_stat$fov == as.numeric(x),]$stats[!is.na(degree_stat[degree_stat$fov == as.numeric(x),]$stats)])
  degree_mean_by_fov[x,1] <<- y})

degree_mean_by_fov <- data.frame(degree_mean_by_fov)
degree_mean_by_fov$level <- ifelse(rownames(degree_mean_by_fov) %in% c(21,25,29,30,37,38,39,40,41), 'low', 'high')

stats <- degree_mean_by_fov %>%
  wilcox_test(mean_degree~level) 

stats$y.position <- max(degree_mean_by_fov$mean_degree) +0.05

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_mean_degree_low_high_wilcoxon.pdf', height = 5, width = 3)
ggplot(degree_mean_by_fov, aes(x= level, y=mean_degree)) +
  geom_boxplot(aes(colour=level),width=0.4, lwd=0.5, outliers = FALSE) +
  geom_point(aes(colour=level), position=position_jitter(seed=1,width = .3),alpha = 1, size=1.5) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  xlab("") +
  stat_pvalue_manual(stats,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  ylab("Mean degree centrality") 

dev.off()

## mean degree score without considering B, Breg, mast, macrophage, T, enteric glia
degree_mean_by_fov <- matrix(nrow = 19, ncol = 1)
colnames(degree_mean_by_fov) <- 'mean_degree'
rownames(degree_mean_by_fov) <- c(rownames(degree_high), rownames(degree_low))

lapply(as.list(levels(degree_stat$fov)), function (x) {
  y <- mean(degree_stat[(degree_stat$fov == as.numeric(x)) & (!degree_stat$cell_type %in% c("Enteric.glia.cells", 'B','B.reg', 'T', 'Macrophage', 
                                                                                             'Mast.cells')),]$stats[!is.na(degree_stat[(degree_stat$fov == as.numeric(x)) & (!degree_stat$cell_type %in% c("Enteric.glia.cells", 'B','B.reg', 'T', 'Macrophage', 
                                                                                                                                                                                                           'Mast.cells')),]$stats)])
  degree_mean_by_fov[x,1] <<- y})

degree_mean_by_fov <- data.frame(degree_mean_by_fov)
degree_mean_by_fov$level <- ifelse(rownames(degree_mean_by_fov) %in% c(21,25,29,30,37,38,39,40,41), 'low', 'high')

stats <- degree_mean_by_fov %>%
  wilcox_test(mean_degree~level) 

stats$y.position <- max(degree_mean_by_fov$mean_degree) +0.05

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_mean_degree_low_high_wilcoxon_without_immune.pdf', height = 5, width = 3)
ggplot(degree_mean_by_fov, aes(x= level, y=mean_degree)) +
  geom_boxplot(aes(colour=level),width=0.4, lwd=0.5, outliers = FALSE) +
  geom_point(aes(colour=level), position=position_jitter(seed=1,width = .3),alpha = 1, size=1.5) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  xlab("") +
  stat_pvalue_manual(stats,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  ylab("Mean degree centrality") 

dev.off()

## mean degree score only considering fibro, muscle, epi and endo
degree_mean_by_fov <- matrix(nrow = 19, ncol = 1)
colnames(degree_mean_by_fov) <- 'mean_degree'
rownames(degree_mean_by_fov) <- c(rownames(degree_high), rownames(degree_low))

lapply(as.list(levels(degree_stat$fov)), function (x) {
  y <- mean(degree_stat[(degree_stat$fov == as.numeric(x)) & (degree_stat$cell_type %in% c('Fibroblast', 'Epithelial', 'Endothelial', "Smooth.muscle.cells")),]$stats[!is.na(degree_stat[(degree_stat$fov == as.numeric(x)) & (degree_stat$cell_type %in% c('Fibroblast', 'Epithelial', 'Endothelial', "Smooth.muscle.cells")),]$stats)])
  degree_mean_by_fov[x,1] <<- y})

degree_mean_by_fov <- data.frame(degree_mean_by_fov)
degree_mean_by_fov$level <- ifelse(rownames(degree_mean_by_fov) %in% c(21,25,29,30,37,38,39,40,41), 'low', 'high')

stats <- degree_mean_by_fov %>%
  wilcox_test(mean_degree~level) 

stats$y.position <- max(degree_mean_by_fov$mean_degree) +0.05

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_mean_degree_low_high_wilcoxon_fibro_endo_epi_muscle.pdf', height = 5, width = 3)
ggplot(degree_mean_by_fov, aes(x= level, y=mean_degree)) +
  geom_boxplot(aes(colour=level),width=0.4, lwd=0.5, outliers = FALSE) +
  geom_point(aes(colour=level), position=position_jitter(seed=1,width = .3),alpha = 1, size=1.5) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  xlab("") +
  stat_pvalue_manual(stats,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  ylab("Mean degree centrality") 

dev.off()


##

average_low <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/average_per_fov_low.xlsx'))
average_high <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/average_per_fov_high.xlsx'))

rownames(average_low) <- average_low$...1
rownames(average_high) <- average_high$...1
average_low[,1] <- NULL
average_high[,1] <- NULL
average_high <- average_high[,colnames(average_low)]
average_stat <- data.frame(stats = c(as.vector(as.matrix(average_low)),as.vector(as.matrix(average_high))),
                          cell_type = c(rep(colnames(average_low), each = nrow(average_low)), rep(colnames(average_high), each = nrow(average_high))),
                          level = c(rep('low', length(rep(colnames(average_low), each = nrow(average_low)))), rep('high', length(rep(colnames(average_high), each = nrow(average_high))))),
                          fov = c(rep(rownames(average_low), ncol(average_low)), rep(rownames(average_high), ncol(average_high))))

average_stat$stats <- as.numeric(average_stat$stats)

stats <- average_stat %>%
  group_by(cell_type) %>%
  wilcox_test(stats~level) %>%
  adjust_pvalue(method = 'fdr')

stats$y.position <- unlist(lapply(seq_along(1:length(stats$cell_type)), function(x) {
  max(average_stat[average_stat$cell_type == stats$cell_type[x],'stats'][! is.na(average_stat[average_stat$cell_type == stats$cell_type[x],'stats'])]) +max(average_stat[average_stat$cell_type == stats$cell_type[x],'stats'][! is.na(average_stat[average_stat$cell_type == stats$cell_type[x],'stats'])])/12
}))

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_average_low_high.pdf', height = 7, width = 10)
ggplot(average_stat, aes(x= level, y=stats)) +
  geom_boxplot(aes(colour=level),width=0.4, lwd=0.5, outliers = FALSE) +
  geom_point(aes(colour=level), position=position_jitter(seed=1,width = .3),alpha = 1, size=1.5) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  xlab("") +
  stat_pvalue_manual(stats,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_wrap(~cell_type, nrow = 2, ncol = 6) +
  ylab("Average clustering") 

dev.off()

average_stat$fov <- factor(average_stat$fov, levels= c(21,25,29,30,37,38,39,40,41,23,24,27,28,31,32,33,34,35,36))
pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/grouped_barplot_clustering_coef_low_high.pdf', height = 7, width = 15)
ggplot(average_stat, aes(fill=fov, y=stats, x=cell_type)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c(rep('#AA3377',9), rep('#0F2080',10))) +
  theme_classic() +
  xlab('Cell type')+
  ylab('Clustering coefficient') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

## mean clustering coef considering all the cell types
clustering_mean_by_fov <- matrix(nrow = 19, ncol = 1)
colnames(clustering_mean_by_fov) <- 'mean_clustering_coefficient'
rownames(clustering_mean_by_fov) <- c(rownames(average_high), rownames(average_low))

lapply(as.list(levels(average_stat$fov)), function (x) {
  y <- mean(average_stat[average_stat$fov == as.numeric(x),]$stats[!is.na(average_stat[average_stat$fov == as.numeric(x),]$stats)])
  clustering_mean_by_fov[x,1] <<- y})

clustering_mean_by_fov <- data.frame(clustering_mean_by_fov)
clustering_mean_by_fov$level <- ifelse(rownames(clustering_mean_by_fov) %in% c(21,25,29,30,37,38,39,40,41), 'low', 'high')

stats <- clustering_mean_by_fov %>%
  wilcox_test(mean_clustering_coefficient~level) 

stats$y.position <- max(clustering_mean_by_fov$mean_clustering_coefficient) +0.01

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_mean_clustering_low_high_wilcoxon.pdf', height = 5, width = 3)
ggplot(clustering_mean_by_fov, aes(x= level, y=mean_clustering_coefficient)) +
  geom_boxplot(aes(colour=level),width=0.4, lwd=0.5, outliers = FALSE) +
  geom_point(aes(colour=level), position=position_jitter(seed=1,width = .3),alpha = 1, size=1.5) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  xlab("") +
  stat_pvalue_manual(stats,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  ylab("Mean clustering coefficient") 

dev.off()

## mean clustering coef without considering B, Breg, mast, macrophage, T, enteric glia
clustering_mean_by_fov <- matrix(nrow = 19, ncol = 1)
colnames(clustering_mean_by_fov) <- 'mean_clustering_coefficient'
rownames(clustering_mean_by_fov) <- c(rownames(average_high), rownames(average_low))

lapply(as.list(levels(average_stat$fov)), function (x) {
  y <- mean(average_stat[(average_stat$fov == as.numeric(x)) & (!average_stat$cell_type %in% c("Enteric.glia.cells", 'B','B.reg', 'T', 'Macrophage', 
                                                                                            'Mast.cells')),]$stats[!is.na(average_stat[(average_stat$fov == as.numeric(x)) & (!average_stat$cell_type %in% c("Enteric.glia.cells", 'B','B.reg', 'T', 'Macrophage', 
                                                                                                                                                                                                          'Mast.cells')),]$stats)])
  clustering_mean_by_fov[x,1] <<- y})

clustering_mean_by_fov <- data.frame(clustering_mean_by_fov)
clustering_mean_by_fov$level <- ifelse(rownames(clustering_mean_by_fov) %in% c(21,25,29,30,37,38,39,40,41), 'low', 'high')

stats <- clustering_mean_by_fov %>%
  wilcox_test(mean_clustering_coefficient~level) 

stats$y.position <- max(clustering_mean_by_fov$mean_clustering_coefficient) +0.01

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_mean_clustering_low_high_wilcoxon_without_immune.pdf', height = 5, width = 3)
ggplot(clustering_mean_by_fov, aes(x= level, y=mean_clustering_coefficient)) +
  geom_boxplot(aes(colour=level),width=0.4, lwd=0.5, outliers = FALSE) +
  geom_point(aes(colour=level), position=position_jitter(seed=1,width = .3),alpha = 1, size=1.5) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  xlab("") +
  stat_pvalue_manual(stats,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  ylab("Mean clustering coefficient") 

dev.off()


## mean clustering coef only considering fibro, muscle, epi and endo
clustering_mean_by_fov <- matrix(nrow = 19, ncol = 1)
colnames(clustering_mean_by_fov) <- 'mean_clustering_coefficient'
rownames(clustering_mean_by_fov) <- c(rownames(average_high), rownames(average_low))

lapply(as.list(levels(average_stat$fov)), function (x) {
  y <- mean(average_stat[(average_stat$fov == as.numeric(x)) & (!average_stat$cell_type %in% c('Fibroblast', 'Epithelial', 'Endothelial', 
                                                                                               "Smooth.muscle.cells")),]$stats[!is.na(average_stat[(average_stat$fov == as.numeric(x)) & (!average_stat$cell_type %in% c('Fibroblast', 'Epithelial', 'Endothelial', "Smooth.muscle.cells")),]$stats)])
  clustering_mean_by_fov[x,1] <<- y})

clustering_mean_by_fov <- data.frame(clustering_mean_by_fov)
clustering_mean_by_fov$level <- ifelse(rownames(clustering_mean_by_fov) %in% c(21,25,29,30,37,38,39,40,41), 'low', 'high')

stats <- clustering_mean_by_fov %>%
  wilcox_test(mean_clustering_coefficient~level) 

stats$y.position <- max(clustering_mean_by_fov$mean_clustering_coefficient) +0.01

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_mean_clustering_low_high_wilcoxon_fibro_endo_epi_muscle.pdf', height = 5, width = 3)
ggplot(clustering_mean_by_fov, aes(x= level, y=mean_clustering_coefficient)) +
  geom_boxplot(aes(colour=level),width=0.4, lwd=0.5, outliers = FALSE) +
  geom_point(aes(colour=level), position=position_jitter(seed=1,width = .3),alpha = 1, size=1.5) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  xlab("") +
  stat_pvalue_manual(stats,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  ylab("Mean clustering coefficient") 

dev.off()


##

closeness_low <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/closeness_per_fov_low.xlsx'))
closeness_high <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/closeness_per_fov_high.xlsx'))

rownames(closeness_low) <- closeness_low$...1
rownames(closeness_high) <- closeness_high$...1
closeness_low[,1] <- NULL
closeness_high[,1] <- NULL
closeness_high <- closeness_high[,colnames(closeness_low)]

closeness_stat <- data.frame(stats = c(as.vector(as.matrix(closeness_low)),as.vector(as.matrix(closeness_high))),
                           cell_type = c(rep(colnames(closeness_low), each = nrow(closeness_low)), rep(colnames(closeness_high), each = nrow(closeness_high))),
                           level = c(rep('low', length(rep(colnames(closeness_low), each = nrow(closeness_low)))), rep('high', length(rep(colnames(closeness_high), each = nrow(closeness_high))))),
                           fov = c(rep(rownames(closeness_low), ncol(closeness_low)), rep(rownames(closeness_high), ncol(closeness_high))))

closeness_stat$stats <- as.numeric(closeness_stat$stats)

stats <- closeness_stat %>%
  group_by(cell_type) %>%
  wilcox_test(stats~level) %>%
  adjust_pvalue(method = 'fdr')

stats$y.position <- unlist(lapply(seq_along(1:length(stats$cell_type)), function(x) {
  max(closeness_stat[closeness_stat$cell_type == stats$cell_type[x],'stats'][! is.na(closeness_stat[closeness_stat$cell_type == stats$cell_type[x],'stats'])]) +max(closeness_stat[closeness_stat$cell_type == stats$cell_type[x],'stats'][! is.na(closeness_stat[closeness_stat$cell_type == stats$cell_type[x],'stats'])])/12
}))

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_closeness_low_high.pdf', height = 7, width = 10)
ggplot(closeness_stat, aes(x= level, y=stats)) +
  geom_boxplot(aes(colour=level),width=0.4, lwd=0.5, outliers = FALSE) +
  geom_point(aes(colour=level), position=position_jitter(seed=1,width = .3),alpha = 1, size=1.5) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  xlab("") +
  stat_pvalue_manual(stats,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_wrap(~cell_type, nrow = 2, ncol = 6) +
  ylab("Closeness centrality") 

dev.off()

closeness_stat$fov <- factor(closeness_stat$fov, levels= c(21,25,29,30,37,38,39,40,41,23,24,27,28,31,32,33,34,35,36))
pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/grouped_barplot_closeness_low_high.pdf', height = 7, width = 15)
ggplot(closeness_stat, aes(fill=fov, y=stats, x=cell_type)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = c(rep('#AA3377',9), rep('#0F2080',10))) +
  theme_classic() +
  xlab('Cell type')+
  ylab('Closeness centrality') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

## mean closeness considering all the cell types
closeness_mean_by_fov <- matrix(nrow = 19, ncol = 1)
colnames(closeness_mean_by_fov) <- 'mean_closeness_centrality'
rownames(closeness_mean_by_fov) <- c(rownames(closeness_high), rownames(closeness_low))

lapply(as.list(levels(closeness_stat$fov)), function (x) {
  y <- mean(closeness_stat[closeness_stat$fov == as.numeric(x),]$stats[!is.na(closeness_stat[closeness_stat$fov == as.numeric(x),]$stats)])
  closeness_mean_by_fov[x,1] <<- y})

closeness_mean_by_fov <- data.frame(closeness_mean_by_fov)
closeness_mean_by_fov$level <- ifelse(rownames(closeness_mean_by_fov) %in% c(21,25,29,30,37,38,39,40,41), 'low', 'high')

stats <- closeness_mean_by_fov %>%
  wilcox_test(mean_closeness_centrality~level) 

stats$y.position <- max(closeness_mean_by_fov$mean_closeness_centrality) +0.01

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_mean_closeness_low_high_wilcoxon.pdf', height = 5, width = 3)
ggplot(closeness_mean_by_fov, aes(x= level, y=mean_closeness_centrality)) +
  geom_boxplot(aes(colour=level),width=0.4, lwd=0.5, outliers = FALSE) +
  geom_point(aes(colour=level), position=position_jitter(seed=1,width = .3),alpha = 1, size=1.5) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  xlab("") +
  stat_pvalue_manual(stats,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  ylab("Mean closeness centrality") 

dev.off()

## mean closeness without considering B, Breg, mast, macrophage, T, enteric glia
closeness_mean_by_fov <- matrix(nrow = 19, ncol = 1)
colnames(closeness_mean_by_fov) <- 'mean_closeness_centrality'
rownames(closeness_mean_by_fov) <- c(rownames(closeness_high), rownames(closeness_low))

lapply(as.list(levels(closeness_stat$fov)), function (x) {
  y <- mean(closeness_stat[(closeness_stat$fov == as.numeric(x)) & (!closeness_stat$cell_type %in% c("Enteric.glia.cells", 'B','B.reg', 'T', 'Macrophage', 
                                                                                               'Mast.cells')),]$stats[!is.na(closeness_stat[(closeness_stat$fov == as.numeric(x)) & (!closeness_stat$cell_type %in% c("Enteric.glia.cells", 'B','B.reg', 'T', 'Macrophage', 
                                                                                                                                                                                                                'Mast.cells')),]$stats)])
  closeness_mean_by_fov[x,1] <<- y})

closeness_mean_by_fov <- data.frame(closeness_mean_by_fov)
closeness_mean_by_fov$level <- ifelse(rownames(closeness_mean_by_fov) %in% c(21,25,29,30,37,38,39,40,41), 'low', 'high')

stats <- closeness_mean_by_fov %>%
  wilcox_test(mean_closeness_centrality~level) 

stats$y.position <- max(closeness_mean_by_fov$mean_closeness_centrality) +0.01

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_mean_closeness_low_high_wilcoxon_without_immune.pdf', height = 5, width = 3)
ggplot(closeness_mean_by_fov, aes(x= level, y=mean_closeness_centrality)) +
  geom_boxplot(aes(colour=level),width=0.4, lwd=0.5, outliers = FALSE) +
  geom_point(aes(colour=level), position=position_jitter(seed=1,width = .3),alpha = 1, size=1.5) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  xlab("") +
  stat_pvalue_manual(stats,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  ylab("Mean closeness centrality") 

dev.off()


## mean closeness only considering fibro, muscle, epi and endo
closeness_mean_by_fov <- matrix(nrow = 19, ncol = 1)
colnames(closeness_mean_by_fov) <- 'mean_closeness_centrality'
rownames(closeness_mean_by_fov) <- c(rownames(closeness_high), rownames(closeness_low))

lapply(as.list(levels(closeness_stat$fov)), function (x) {
  y <- mean(closeness_stat[(closeness_stat$fov == as.numeric(x)) & (!closeness_stat$cell_type %in% c('Fibroblast', 'Epithelial', 'Endothelial', 
                                                                                               "Smooth.muscle.cells")),]$stats[!is.na(closeness_stat[(closeness_stat$fov == as.numeric(x)) & (!closeness_stat$cell_type %in% c('Fibroblast', 'Epithelial', 'Endothelial', "Smooth.muscle.cells")),]$stats)])
  closeness_mean_by_fov[x,1] <<- y})

closeness_mean_by_fov <- data.frame(closeness_mean_by_fov)
closeness_mean_by_fov$level <- ifelse(rownames(closeness_mean_by_fov) %in% c(21,25,29,30,37,38,39,40,41), 'low', 'high')

stats <- closeness_mean_by_fov %>%
  wilcox_test(mean_closeness_centrality~level) 

stats$y.position <- max(closeness_mean_by_fov$mean_closeness_centrality) +0.01

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_mean_closeness_low_high_wilcoxon_fibro_endo_epi_muscle.pdf', height = 5, width = 3)
ggplot(closeness_mean_by_fov, aes(x= level, y=mean_closeness_centrality)) +
  geom_boxplot(aes(colour=level),width=0.4, lwd=0.5, outliers = FALSE) +
  geom_point(aes(colour=level), position=position_jitter(seed=1,width = .3),alpha = 1, size=1.5) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  xlab("") +
  stat_pvalue_manual(stats,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  ylab("Mean closeness centrality") 

dev.off()


## neighboorhod enrichment analysis

neighboorhod_results <- data.frame(matrix(ncol=4))
colnames(neighboorhod_results) <- c('cell_type_pair', 'score', 'level','fov')

for (f in list.files('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/neighboorhod/', full.names = TRUE)) {
  scores <- data.frame(read_excel(f),row.names = 1)
  
  fov <- as.numeric(gsub('.xlsx','',gsub('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/neighboorhod//neighboorhod_','',f)))
  
  level <- ifelse(fov %in% c(21,25,29,30,37,38,39,40,41), 'low', 'high')
  
  for (r in rownames(scores)) {
    for (c in colnames(scores)) {
      if (! is.na(scores[r,c])) {
        new_row <- c(paste(r,'-',c),as.character(scores[r,c]), level, as.character(fov))
        neighboorhod_results <- rbind(neighboorhod_results, new_row)
      }
    }    
  }
}

neighboorhod_results <- neighboorhod_results[-1,]
neighboorhod_results$score <- as.numeric(neighboorhod_results$score)

stats <- neighboorhod_results %>%
  group_by(cell_type_pair) %>%
  wilcox_test(score~level) %>%
  adjust_pvalue(method = 'fdr')

stats$y.position <- unlist(lapply(seq_along(1:length(stats$cell_type_pair)), function(x) {
  max(neighboorhod_results[neighboorhod_results$cell_type_pair == stats$cell_type_pair[x],'score'][! is.na(neighboorhod_results[neighboorhod_results$cell_type_pair == stats$cell_type_pair[x],'score'])]) +max(neighboorhod_results[neighboorhod_results$cell_type_pair == stats$cell_type_pair[x],'score'][! is.na(neighboorhod_results[neighboorhod_results$cell_type_pair == stats$cell_type_pair[x],'score'])])/10
}))

stats_significant <- stats[stats$p < 0.1,]
stats_significant2 <- stats_significant[!stats_significant$cell_type_pair %in% c('Epithelial - Epithelial', 'Smooth muscle cells - Smooth.muscle.cells'),]
stats_significant3 <- stats_significant[stats_significant$cell_type_pair %in% c('Epithelial - Epithelial', 'Smooth muscle cells - Smooth.muscle.cells'),]

neighboorhod_results_sign <-neighboorhod_results[neighboorhod_results$cell_type_pair %in% stats_significant$cell_type_pair,]
neighboorhod_results_sign2 <-neighboorhod_results[neighboorhod_results$cell_type_pair %in% stats_significant2$cell_type_pair,]
neighboorhod_results_sign3 <-neighboorhod_results[neighboorhod_results$cell_type_pair %in% stats_significant3$cell_type_pair,]

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_neigh_analysis_low_high.pdf', height = 5, width = 10)
ggplot(neighboorhod_results_sign2, aes(x= level, y=score)) +
  geom_boxplot(aes(colour=level),width=0.4, lwd=0.5, outliers = FALSE) +
  geom_point(aes(colour=level), position=position_jitter(seed=1,width = .3),alpha = 1, size=1.5) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  xlab("") +
  stat_pvalue_manual(stats_significant2,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_wrap(~cell_type_pair, nrow = 2, ncol = 6) +
  ylab("Neighboorhod enrichment analysis") +
  theme(strip.text=element_text(size = 8))

dev.off()

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/boxplot_neigh_analysis_low_high2.pdf', height =3, width = 5)
ggplot(neighboorhod_results_sign3, aes(x= level, y=score)) +
  geom_boxplot(aes(colour=level),width=0.4, lwd=0.5, outliers = FALSE) +
  geom_point(aes(colour=level), position=position_jitter(seed=1,width = .3),alpha = 1, size=1.5) +
  scale_colour_manual(values = c("#ffb703", "#f94144")) +
  theme_classic() +
  scale_shape(solid=T)  +
  xlab("") +
  stat_pvalue_manual(stats_significant3,label = "p", tip.length=0.01 ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_wrap(~cell_type_pair, nrow = 1, ncol = 2) +
  ylab("Neighboorhod enrichment analysis") +
  theme(strip.text=element_text(size = 8))

dev.off()



## novae 10 domains

# check how many niches are inside each fov and if there are differences between high and low
metadata_niche <- read.csv('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/metadata_adata_niche_novae_10domains.csv')

fov_niche_cell_type <- metadata_niche[,c('fov', 'novae_domains_10','new_leiden')]

# plot the proportion of cells belonging to each niche in each fov 
cells_each_niche <- data.frame(fov = rep(unique(fov_niche_cell_type$fov), length(unique(fov_niche_cell_type$novae_domains_10))),
                                             niche = rep(unique(fov_niche_cell_type$novae_domains_10), each = length(unique(fov_niche_cell_type$fov))),
                                             count = NA)

for (f in unique(cells_each_niche$fov)) {
  for (n in unique(cells_each_niche$niche)) {
    cells_each_niche[cells_each_niche$fov ==f & cells_each_niche$niche ==n, 3] <- nrow(fov_niche_cell_type[fov_niche_cell_type$fov == f & fov_niche_cell_type$novae_domains_10 ==n,])
  }
}

cells_each_niche$fov <- factor(cells_each_niche$fov)

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/amount_of_niche_each_fov_novae_10domains.pdf')
ggplot(cells_each_niche, aes(fill=niche, y=count, x=fov)) + 
  geom_bar(position=position_fill(reverse = TRUE), stat="identity") +
  labs(y = "% of niches in each fov") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.text=element_text(size=8), axis.title.x = element_blank(), 
        legend.title = element_blank()) +
  theme(aspect.ratio=2) 
dev.off()

# plot the composition of each niche
cell_types_each_niche <- data.frame(cell_type = rep(unique(fov_niche_cell_type$new_leiden), length(unique(fov_niche_cell_type$novae_domains_10))),
                               niche = rep(unique(fov_niche_cell_type$novae_domains_10), each = length(unique(fov_niche_cell_type$new_leiden))),
                               count = NA)

for (c in unique(cell_types_each_niche$cell_type)) {
  for (n in unique(cell_types_each_niche$niche)) {
    cell_types_each_niche[cell_types_each_niche$cell_type ==c & cell_types_each_niche$niche ==n, 3] <- nrow(fov_niche_cell_type[fov_niche_cell_type$new_leiden == c & fov_niche_cell_type$novae_domains_10 ==n,])
  }
}

cell_types_each_niche$cell_type <- factor(cell_types_each_niche$cell_type)

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/cellular_composition_each_niche_novae_10domains.pdf')
ggplot(cell_types_each_niche, aes(fill=cell_type, y=count, x=niche)) + 
  geom_bar(position=position_fill(reverse = TRUE), stat="identity") +
  scale_fill_manual(values = c('B' = '#0F2080', 'B reg' = '#228833', 'Endothelial' = '#AA3377', 'Enteric glia cells' = '#BBBBBB', 'Epithelial' = '#AEC7E8',
                      'Fibroblast' = '#CCBB44','Macrophage' = '#F7B6D2', 'Mast cells' = '#0077BB', 
                      'Plasma IgA' = '#EE7733', 'Plasma IgG'='#33BBEE','Smooth muscle cells' = '#CC3311','T' = '#7F7F7F')) +
  labs(y = "% of cell types in each niche") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.text=element_text(size=8), axis.title.x = element_blank(), 
        legend.title = element_blank()) +
  theme(aspect.ratio=2) 
dev.off()

# check the proportion of cells belonging to each niche in each fov and if there are differences between high and low
proportion_of_cells_each_niche <- data.frame(fov = rep(unique(fov_niche_cell_type$fov), length(unique(fov_niche_cell_type$novae_domains_10))),
                                             niche = rep(unique(fov_niche_cell_type$novae_domains_10), each = length(unique(fov_niche_cell_type$fov))),
                                             prop = NA)

for (f in unique(proportion_of_cells_each_niche$fov)) {
  for (n in unique(proportion_of_cells_each_niche$niche)) {
    proportion_of_cells_each_niche[proportion_of_cells_each_niche$fov ==f & proportion_of_cells_each_niche$niche ==n, 3] <- nrow(fov_niche_cell_type[fov_niche_cell_type$fov == f & fov_niche_cell_type$novae_domains_10 ==n,])/nrow(fov_niche_cell_type[fov_niche_cell_type$fov ==f,])
  }
}

proportion_of_cells_each_niche$level <- unlist(lapply(as.list(proportion_of_cells_each_niche$fov), function(x) {ifelse(x %in% c(21,25,29,30,37,38,39,40,41), 'low', 'high')}))

stats <- proportion_of_cells_each_niche %>%
  group_by(niche) %>%
  wilcox_test(prop~level)



## novae 10 domains

# check how many niches are inside each fov and if there are differences between high and low
metadata_niche <- read.csv('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/metadata_adata_niche_novae_10domains.csv')

fov_niche_cell_type <- metadata_niche[,c('fov', 'novae_domains_10','new_leiden')]

# plot the proportion of cells belonging to each niche in each fov 
cells_each_niche <- data.frame(fov = rep(unique(fov_niche_cell_type$fov), length(unique(fov_niche_cell_type$novae_domains_10))),
                               niche = rep(unique(fov_niche_cell_type$novae_domains_10), each = length(unique(fov_niche_cell_type$fov))),
                               count = NA)

for (f in unique(cells_each_niche$fov)) {
  for (n in unique(cells_each_niche$niche)) {
    cells_each_niche[cells_each_niche$fov ==f & cells_each_niche$niche ==n, 3] <- nrow(fov_niche_cell_type[fov_niche_cell_type$fov == f & fov_niche_cell_type$novae_domains_10 ==n,])
  }
}

cells_each_niche$fov <- factor(cells_each_niche$fov)

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/amount_of_niche_each_fov_novae_10domains.pdf')
ggplot(cells_each_niche, aes(fill=niche, y=count, x=fov)) + 
  geom_bar(position=position_fill(reverse = TRUE), stat="identity") +
  labs(y = "% of niches in each fov") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.text=element_text(size=8), axis.title.x = element_blank(), 
        legend.title = element_blank()) +
  theme(aspect.ratio=2) 
dev.off()

# plot the composition of each niche
cell_types_each_niche <- data.frame(cell_type = rep(unique(fov_niche_cell_type$new_leiden), length(unique(fov_niche_cell_type$novae_domains_10))),
                                    niche = rep(unique(fov_niche_cell_type$novae_domains_10), each = length(unique(fov_niche_cell_type$new_leiden))),
                                    count = NA)

for (c in unique(cell_types_each_niche$cell_type)) {
  for (n in unique(cell_types_each_niche$niche)) {
    cell_types_each_niche[cell_types_each_niche$cell_type ==c & cell_types_each_niche$niche ==n, 3] <- nrow(fov_niche_cell_type[fov_niche_cell_type$new_leiden == c & fov_niche_cell_type$novae_domains_10 ==n,])
  }
}

cell_types_each_niche$cell_type <- factor(cell_types_each_niche$cell_type)

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/cellular_composition_each_niche_novae_10domains.pdf')
ggplot(cell_types_each_niche, aes(fill=cell_type, y=count, x=niche)) + 
  geom_bar(position=position_fill(reverse = TRUE), stat="identity") +
  scale_fill_manual(values = c('B' = '#0F2080', 'B reg' = '#228833', 'Endothelial' = '#AA3377', 'Enteric glia cells' = '#BBBBBB', 'Epithelial' = '#AEC7E8',
                               'Fibroblast' = '#CCBB44','Macrophage' = '#F7B6D2', 'Mast cells' = '#0077BB', 
                               'Plasma IgA' = '#EE7733', 'Plasma IgG'='#33BBEE','Smooth muscle cells' = '#CC3311','T' = '#7F7F7F')) +
  labs(y = "% of cell types in each niche") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.text=element_text(size=8), axis.title.x = element_blank(), 
        legend.title = element_blank()) +
  theme(aspect.ratio=2) 
dev.off()

# check the proportion of cells belonging to each niche in each fov and if there are differences between high and low
proportion_of_cells_each_niche <- data.frame(fov = rep(unique(fov_niche_cell_type$fov), length(unique(fov_niche_cell_type$novae_domains_10))),
                                             niche = rep(unique(fov_niche_cell_type$novae_domains_10), each = length(unique(fov_niche_cell_type$fov))),
                                             prop = NA)

for (f in unique(proportion_of_cells_each_niche$fov)) {
  for (n in unique(proportion_of_cells_each_niche$niche)) {
    proportion_of_cells_each_niche[proportion_of_cells_each_niche$fov ==f & proportion_of_cells_each_niche$niche ==n, 3] <- nrow(fov_niche_cell_type[fov_niche_cell_type$fov == f & fov_niche_cell_type$novae_domains_10 ==n,])/nrow(fov_niche_cell_type[fov_niche_cell_type$fov ==f,])
  }
}

proportion_of_cells_each_niche$level <- unlist(lapply(as.list(proportion_of_cells_each_niche$fov), function(x) {ifelse(x %in% c(21,25,29,30,37,38,39,40,41), 'low', 'high')}))

stats <- proportion_of_cells_each_niche %>%
  group_by(niche) %>%
  wilcox_test(prop~level)


## novae 20 domains

# check how many niches are inside each fov and if there are differences between high and low
metadata_niche <- read.csv('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/metadata_adata_niche_novae_20domains.csv')

fov_niche_cell_type <- metadata_niche[,c('fov', 'novae_domains_20','new_leiden')]

# plot the proportion of cells belonging to each niche in each fov 
cells_each_niche <- data.frame(fov = rep(unique(fov_niche_cell_type$fov), length(unique(fov_niche_cell_type$novae_domains_20))),
                               niche = rep(unique(fov_niche_cell_type$novae_domains_20), each = length(unique(fov_niche_cell_type$fov))),
                               count = NA)

for (f in unique(cells_each_niche$fov)) {
  for (n in unique(cells_each_niche$niche)) {
    cells_each_niche[cells_each_niche$fov ==f & cells_each_niche$niche ==n, 3] <- nrow(fov_niche_cell_type[fov_niche_cell_type$fov == f & fov_niche_cell_type$novae_domains_20 ==n,])
  }
}

cells_each_niche$fov <- factor(cells_each_niche$fov)

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/amount_of_niche_each_fov_novae_15domains.pdf')
ggplot(cells_each_niche, aes(fill=niche, y=count, x=fov)) + 
  geom_bar(position=position_fill(reverse = TRUE), stat="identity") +
  labs(y = "% of cells belonging to each niche in each fov") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.text=element_text(size=8), axis.title.x = element_blank(), 
        legend.title = element_blank()) +
  theme(aspect.ratio=2) 
dev.off()

# plot the composition of each niche
cell_types_each_niche <- data.frame(cell_type = rep(unique(fov_niche_cell_type$new_leiden), length(unique(fov_niche_cell_type$novae_domains_20))),
                                    niche = rep(unique(fov_niche_cell_type$novae_domains_20), each = length(unique(fov_niche_cell_type$new_leiden))),
                                    count = NA)

for (c in unique(cell_types_each_niche$cell_type)) {
  for (n in unique(cell_types_each_niche$niche)) {
    cell_types_each_niche[cell_types_each_niche$cell_type ==c & cell_types_each_niche$niche ==n, 3] <- nrow(fov_niche_cell_type[fov_niche_cell_type$new_leiden == c & fov_niche_cell_type$novae_domains_20 ==n,])
  }
}

cell_types_each_niche$cell_type <- factor(cell_types_each_niche$cell_type)

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/cellular_composition_each_niche_novae_20domains.pdf')
ggplot(cell_types_each_niche, aes(fill=cell_type, y=count, x=niche)) + 
  geom_bar(position=position_fill(reverse = TRUE), stat="identity") +
  scale_fill_manual(values = c('B' = '#0F2080', 'B reg' = '#228833', 'Endothelial' = '#AA3377', 'Enteric glia cells' = '#BBBBBB', 'Epithelial' = '#AEC7E8',
                               'Fibroblast' = '#CCBB44','Macrophage' = '#F7B6D2', 'Mast cells' = '#0077BB', 
                               'Plasma IgA' = '#EE7733', 'Plasma IgG'='#33BBEE','Smooth muscle cells' = '#CC3311','T' = '#7F7F7F')) +
  labs(y = "% of cells belonging to each\ncell type in each niche") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.text=element_text(size=8), axis.title.x = element_blank(), 
        legend.title = element_blank()) +
  theme(aspect.ratio=2) 
dev.off()

# check the proportion of cells belonging to each niche in each fov and if there are differences between high and low
proportion_of_cells_each_niche <- data.frame(fov = rep(unique(fov_niche_cell_type$fov), length(unique(fov_niche_cell_type$novae_domains_20))),
                                             niche = rep(unique(fov_niche_cell_type$novae_domains_20), each = length(unique(fov_niche_cell_type$fov))),
                                             prop = NA)

for (f in unique(proportion_of_cells_each_niche$fov)) {
  for (n in unique(proportion_of_cells_each_niche$niche)) {
    proportion_of_cells_each_niche[proportion_of_cells_each_niche$fov ==f & proportion_of_cells_each_niche$niche ==n, 3] <- nrow(fov_niche_cell_type[fov_niche_cell_type$fov == f & fov_niche_cell_type$novae_domains_20 ==n,])/nrow(fov_niche_cell_type[fov_niche_cell_type$fov ==f,])
  }
}

proportion_of_cells_each_niche$level <- unlist(lapply(as.list(proportion_of_cells_each_niche$fov), function(x) {ifelse(x %in% c(21,25,29,30,37,38,39,40,41), 'low', 'high')}))

stats_greater_high <- proportion_of_cells_each_niche %>%
  group_by(niche) %>%
  wilcox_test(prop~level, alternative = 'greater',ref.group = 'high')

stats_greater_low <- proportion_of_cells_each_niche %>%
  group_by(niche) %>%
  wilcox_test(prop~level, alternative = 'greater',ref.group = 'low')

a <- aggregate(proportion_of_cells_each_niche$prop, list(proportion_of_cells_each_niche$niche,proportion_of_cells_each_niche$level), FUN=mean) 

colnames(a) <- c('Niche','Level', 'prop')
write.table(a, file='/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/niche_abundance_high_low.txt', quote = F, sep = ',', row.names = F)

a_sub <- a[a$Niche %in% c('D1002','D980','D993','D997'),]
write.table(a_sub, file='/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/niche_abundance_high_low_only_sig.txt', quote = F, sep = ',', row.names = F)

proportion_of_cells_each_niche_only_sig <- proportion_of_cells_each_niche[proportion_of_cells_each_niche$niche %in% c('D1002','D980','D993','D997'),]
write.table(proportion_of_cells_each_niche_only_sig, file='/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/proportion_of_cells_each_niche_only_sig.txt', quote = F, sep = ',', row.names = F)

stats_sign <- rbind(stats_greater_high[stats_greater_high$niche %in% c('D1002','D993','D997'),],
                    stats_greater_low[stats_greater_low$niche == 'D980',])

stats_sign$y.position <- unlist(lapply(seq_along(1:length(stats_sign$niche)), function(x) {
  max(proportion_of_cells_each_niche_only_sig[proportion_of_cells_each_niche_only_sig$niche == stats_sign$niche[x],'prop'][! is.na(proportion_of_cells_each_niche_only_sig[proportion_of_cells_each_niche_only_sig$niche == stats_sign$niche[x],'prop'])])+max(proportion_of_cells_each_niche_only_sig[proportion_of_cells_each_niche_only_sig$niche == stats_sign$niche[x],'prop'][! is.na(proportion_of_cells_each_niche_only_sig[proportion_of_cells_each_niche_only_sig$niche == stats_sign$niche[x],'prop'])])/12
}))

proportion_of_cells_each_niche_only_sig$niche <- factor(proportion_of_cells_each_niche_only_sig$niche)

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/mean_amount_of_niche_each_fov_novae_20domains.pdf')
ggplot(a, aes(fill=Niche, y=prop, x=Level)) + 
  geom_bar(position=position_fill(reverse = TRUE), stat="identity") +
  labs(y = "mean % of cells belonging to each niche\nin high and low iNKT10 FOVs") +
  theme_classic() +
  scale_fill_manual(values = c('#0F2080', '#006600', '#AA3377', '#BBBBBB', '#AEC7E8', '#CCBB44', '#F7B6D2','#000000',
                               '#0077BB', '#EE7733', '#33BBEE', '#CC3311', '#7F7F7F', '#FFB000','#601A4A','#785EF0','#66CC00', '#CCFF99', 
                               '#660000', '#E5CCFF','#FFFFCC','#CCFFE5')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.text=element_text(size=8), axis.title.x = element_blank(), 
        legend.title = element_blank()) +
  theme(aspect.ratio=2) 
dev.off()

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/amount_of_niche_each_fov_sign_novae_20domains.pdf', height = 7,width = 15)
ggplot(proportion_of_cells_each_niche_only_sig, aes(color=level, y=prop, x=level)) + 
  geom_boxplot(aes(color=level),width=0.4, lwd=1, outliers = F) +
  geom_jitter(aes(color=level), size=2,width = 0.25) +
  scale_color_manual(values = c("#28658f", "#ef857c")) +
  facet_wrap(~niche,nrow = 1, ncol = 4,scale = "free_y") +
  theme_classic() +
  scale_shape(solid=T)  +
  ylab("% of cells belonging to each niche in each fov") +
  xlab("") +
  theme(aspect.ratio = 1.5,strip.text=element_text(size = 25, face = 'bold'), axis.text.x = element_text(size=25, angle = 45,hjust = 1), 
        axis.text.y = element_text(size=25),axis.title.y = element_text(size=20),
        legend.text = element_text(size = 25),legend.title =  element_text(size = 25), strip.text.x = element_text(size=20)) +
  stat_pvalue_manual(stats_sign, label='p', size=7)
dev.off()


### DEGs between enriched niches high-low
df_D997 <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_niche_D997.xlsx'))
df_D997[,1] <- NULL
df_D997$diff <- 'NO'
df_D997[(df_D997$logfoldchanges> 0 & df_D997$pvals<=0.05),]$diff <- 'LOW'
df_D997[(df_D997$logfoldchanges< 0 & df_D997$pvals<=0.05),]$diff <- 'HIGH'

write.table(df_D997, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_D997.txt',col.names =TRUE, row.names =FALSE, quote = FALSE, sep = '\t')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_D997_low_high.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_D997, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c(LOW="#28658f", HIGH="#ef857c")) +
  geom_point(data = subset(df_D997, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals<0.01,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs of niche D997 between high and low iNKT10 level')

dev.off()

write.table(df_D997[df_D997$diff != 'NO','names'], file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_D997_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep = '\n')

df_D997_fdr_0.05_high <- df_D997[df_D997$diff == "HIGH",'names']

from_david_to_cyrcle('david_D997_BP', df_D997_fdr_0.05_high, 'D997')


##
df_D993 <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_niche_D993.xlsx'))
df_D993[,1] <- NULL
df_D993$diff <- 'NO'
df_D993[(df_D993$logfoldchanges> 0 & df_D993$pvals<=0.05),]$diff <- 'LOW'
df_D993[(df_D993$logfoldchanges< 0 & df_D993$pvals<=0.05),]$diff <- 'HIGH'

write.table(df_D993, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_D993.txt',col.names =TRUE, row.names =FALSE, quote = FALSE, sep = '\t')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_D993_low_high.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_D993, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c(LOW="#28658f", HIGH="#ef857c")) +
  geom_point(data = subset(df_D993, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals<0.001,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs of niche D993 between high and low iNKT10 level')

dev.off()

write.table(df_D993[df_D993$diff != 'NO','names'], file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_D993_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep = '\n')

df_D993_fdr_0.05_high <- df_D993[df_D993$diff == "HIGH",'names']

df_D993_fdr_0.05_high <- df_D993_fdr_0.05_high[df_D993_fdr_0.05_high != 'HLA.DRB']
df_D993_fdr_0.05_high[df_D993_fdr_0.05_high == "MZT2A.B"] <- 'MZT2A'
df_D993_fdr_0.05_high <- c(df_D993_fdr_0.05_high, 'MZT2B')
df_D993_fdr_0.05_high[df_D993_fdr_0.05_high == "EIF5A.L1"] <- 'EIF5A'

from_david_to_cyrcle('david_D993_BP', df_D993_fdr_0.05_high, 'D993')


##
df_D1002 <- data.frame(read_excel('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_niche_D1002.xlsx'))
df_D1002[,1] <- NULL
df_D1002$diff <- 'NO'
df_D1002[(df_D1002$logfoldchanges> 0 & df_D1002$pvals<=0.05),]$diff <- 'LOW'
df_D1002[(df_D1002$logfoldchanges< 0 & df_D1002$pvals<=0.05),]$diff <- 'HIGH'

write.table(df_D1002, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_D1002.txt',col.names =TRUE, row.names =FALSE, quote = FALSE, sep = '\t')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/volcano_differential_D1002_low_high.pdf', height = 7, width =10)

ggplot2::ggplot(data=df_D1002, aes(x=logfoldchanges, y=-log10(pvals), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c(LOW="#28658f", HIGH="#ef857c")) +
  geom_point(data = subset(df_D1002, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(pvals<0.01,names,'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs of niche D1002 between high and low iNKT10 level')

dev.off()

write.table(df_D1002[df_D1002$diff != 'NO','names'], file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_D1002_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

df_D1002_fdr_0.05_high <- df_D1002[df_D1002$diff == "HIGH",'names']

from_david_to_cyrcle('david_D1002_BP', df_D1002_fdr_0.05_high, 'D1002')


##
david_BP <- read.table('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/david_D980_BP.txt', header=TRUE, sep='\t')

david_BP$up_in_high <- 0
david_BP$up_in_low <- 0

for (r in 1:nrow(david_BP)) {
  low <- c()
  for (x in unlist(strsplit(david_BP[r,'Genes'], ', '))) {
    low <- c(low, x)
  } 
  david_BP[r,]$up_in_low <- length(low)
}

david_BP2 <- data.frame(ID = david_BP$Term, classification = david_BP$Category, FDR = david_BP$FDR, 
                        all = david_BP$Pop.Hits, enriched_in_HIGH = david_BP$up_in_high, enriched_in_LOW = david_BP$up_in_low)
david_BP2 <- david_BP2[order(david_BP2$FDR),]

if (nrow(david_BP2) > 10) {
  david_BP2 <- david_BP2[1:10,]
}
write.table(david_BP2[,1], file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/david_BP_D980_modified_entire_GO.tsv',col.names =TRUE, row.names =FALSE, quote = FALSE, sep='\t')

david_BP2$ID <- gsub('~.*$', '', david_BP2$ID)

write.table(david_BP2, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/david_BP_D980_modified.tsv',col.names =TRUE, row.names =FALSE, quote = FALSE, sep='\t')

