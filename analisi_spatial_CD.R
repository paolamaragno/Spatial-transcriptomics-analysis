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

# # heatmap
# NanoOBJ <- readRDS('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/data/RNA_seurat_object.Rds')
# 
# new_column_values <- numeric(69526)
# NanoOBJ <- AddMetaData(object = NanoOBJ, metadata = new_column_values, 
#                        col.name = "Patients")
# NanoOBJ$Patients <- ifelse(
#   NanoOBJ$fov %in% c(39,40,41), "Patient_20", ifelse(
#     NanoOBJ$fov %in% c(37,38), "Patient_19", ifelse(
#       NanoOBJ$fov %in% c(35,36), "Patient_18", ifelse(
#         NanoOBJ$fov %in% c(33,34), "Patient_17", ifelse(
#           NanoOBJ$fov %in% c(31,32), "Patient_16", ifelse(
#             NanoOBJ$fov %in% c(29,30), "Patient_15", ifelse(
#               NanoOBJ$fov %in% c(27,28), "Patient_14", ifelse(
#                 NanoOBJ$fov %in% c(25,26), "Patient_13", ifelse(
#                   NanoOBJ$fov %in% c(23,24), "Patient_12", ifelse(
#                     NanoOBJ$fov %in% c(21,22), "Patient_11", ifelse(
#                       NanoOBJ$fov %in% c(18,19,20), "Patient_9", ifelse(
#                         NanoOBJ$fov %in% c(14,15,16,17), "Patient_8", ifelse(
#                           NanoOBJ$fov %in% c(12,13), "Patient_7", ifelse(
#                             NanoOBJ$fov %in% c(9,10,11), "Patient_6", ifelse(
#                               NanoOBJ$fov %in% c(6,7,8), "Patient_4", ifelse(
#                                 NanoOBJ$fov %in% c(4,5), "Patient_2",
#                                 "Patient_1"))))))))))))))))
# NanoOBJ$new_column <- ifelse(
#   NanoOBJ$Patients %in% c('Patient_11', 'Patient_12', 'Patient_13', 'Patient_14',
#                           'Patient_15', 'Patient_16', 'Patient_17', 'Patient_18',
#                           'Patient_19', 'Patient_20'), "CD", "K")
# 
# CD <- subset(NanoOBJ, subset= new_column == 'CD')
# 
# new_column_values <- numeric(69526)
# CD <- AddMetaData(object = CD, metadata = new_column_values, 
#                   col.name = "Level")
# CD$Level <- ifelse(
#   CD$Patients %in% c('Patient_11', 'Patient_13', 'Patient_15', 'Patient_19', 'Patient_20'), "LOW", "HIGH")
# 
# metadata <- read.csv('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/metadata_adata.csv', header = TRUE)
# rownames(metadata) <- metadata[,1]
# metadata[,1] <- NULL
# 
# CD <-  SCTransform(CD, assay = "Nanostring")
# # this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
# 
# CD2 <- CD
# new_column_values <- numeric(69526)
# CD2 <- AddMetaData(object = CD2, metadata = new_column_values, 
#                    col.name = "cell")
# CD2$cell <- gsub('R5941.ColonTMA_', '', rownames(CD2@meta.data))
# CD2 <- subset(CD2, subset= cell %in% rownames(metadata))
# 
# CD2@meta.data$nb_clus <- unlist(lapply(as.list(CD2@meta.data$cell), function(x) {metadata[x,'new_leiden']}))
# 
# a <- names(CD2@active.ident)
# CD2@active.ident <- factor(unlist(lapply(as.list(a), function(x) {CD2@meta.data[x,'nb_clus']})))
# names(CD2@active.ident) <- a
# 
# CD2.markers <- FindAllMarkers(CD2, only.pos = TRUE)
# 
# CD2.markers %>%
#   group_by(cluster) %>%
#   dplyr::filter(avg_log2FC > 1) %>%
#   slice_head(n = 10) %>%
#   ungroup() -> top10
# 
# genes <- c(top10[top10$cluster =='B','gene']$gene,'CD19', 
#           top10[top10$cluster =='B reg','gene']$gene, 'CD44','CD19',
#            top10[top10$cluster =='Endothelial','gene']$gene, 'TIE1',
#            top10[top10$cluster =='Enteric glia cell','gene']$gene,'VIM', 
#            top10[top10$cluster =='Epithelial','gene']$gene,
#            top10[top10$cluster =='Fibroblast','gene']$gene,
#            top10[top10$cluster =='Macrophage','gene']$gene, 'CD68', 'CD163','CLEC10A', 
#            top10[top10$cluster =='Mast','gene']$gene,
#            top10[top10$cluster =='Plasma IgA','gene']$gene,
#           top10[top10$cluster =='Plasma IgG','gene']$gene,
#            top10[top10$cluster =='Smooth muscle cell','gene']$gene,
#            top10[top10$cluster =='T','gene']$gene,'CD3D', 'CD3G',
#           top10[top10$cluster =='Undifferentiated','gene']$gene, 'TUBB', 'MTOR'
# )
# 
# pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/expression_heatmap.pdf', height =16, width=16)
# DoHeatmap(CD2, features = genes, size = 8) + NoLegend()
# dev.off() 

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
  
  # file_KEGG <- read.table(paste0('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/', file_KEGG, '.txt'), header=TRUE, sep='\t')
  # 
  # file_KEGG$up_in_high <- 0
  # file_KEGG$up_in_low <- 0
  # 
  # for (r in 1:nrow(file_KEGG)) {
  #   high <- c()
  #   low <- c()
  #   for (x in unlist(strsplit(file_KEGG[r,'Genes'], ', '))) {
  #     if (x %in% degs_high) {
  #       high <- c(high, x)
  #     } else {
  #       low <- c(low, x)
  #     }
  #   } 
  #   file_KEGG[r,]$up_in_high <- length(high)
  #   file_KEGG[r,]$up_in_low <- length(low)
  # }
  # 
  # file_KEGG2 <- data.frame(ID = file_KEGG$Term, classification = file_KEGG$Category, FDR = file_KEGG$FDR, 
  #                          all = file_KEGG$Pop.Hits, enriched_in_HIGH = file_KEGG$up_in_high, enriched_in_LOW = file_KEGG$up_in_low)
  # file_KEGG2 <- file_KEGG2[order(file_KEGG2$FDR),]
  # 
  # if (nrow(file_KEGG2) > 10) {
  #   file_KEGG2 <- file_KEGG2[1:10,]
  # }
  # 
  # write.table(file_KEGG2[,1], file=paste0('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/david_KEGG_', cell_type,'_modified_GO.tsv'),col.names =TRUE, row.names =FALSE, quote = FALSE, sep='\t')
  # 
  # file_KEGG2$ID <- gsub(':.*$', '', file_KEGG2$ID)
  # 
  # write.table(file_KEGG2, file=paste0('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/david_KEGG_', cell_type,'_modified.tsv'),col.names =TRUE, row.names =FALSE, quote = FALSE, sep='\t')
}

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

## gsea on python
# df_epi <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/GSA_epi_DEGs.txt', sep=',', header = T)
# gsea_results_epi <- read.csv('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/GSA_epi_selected_table.csv')
# gsea_results_epi[,1] <- NULL
# df_epi[,1] <- NULL
# 
# colnames(gsea_results_epi)[11] <- '-log10(Adjusted p-value)'
# 
# gsea_results_epi$up_high <- NA
# gsea_results_epi$up_low <- NA
# 
# for (p in gsea_results_epi$Term) {
#   low <- 0
#   up <- 0
#   for (g in unlist(strsplit(gsea_results_epi[gsea_results_epi$Term == p, 'Genes'], split=';'))) {
#     if (df_epi[df_epi$names == g, "logfoldchanges"] > 0) {
#       low <- low+1}
#     else {
#       up <- up +1
#     }
#   }
#   gsea_results_epi[gsea_results_epi$Term == p, 'up_high'] <- up
#   gsea_results_epi[gsea_results_epi$Term == p, 'up_low'] <- low
#   }
#  
# a <- gsea_results_epi[,c('Term', '-log10(Adjusted p-value)', 'Genes','up_high', 'up_low')]
# 
# 
# df_final <- data.frame(term = rep(gsea_results_epi$Term,2),
#                        '-log10(Adjusted p-value)' = rep(gsea_results_epi$`-log10(Adjusted p-value)`,2),
#                        condition = c(rep('enriched in high', length(gsea_results_epi$up_high)),rep('enriched in low', length(gsea_results_epi$up_low))),
#                        count = c(gsea_results_epi$up_high,gsea_results_epi$up_low))
# 
# df_final$term <- factor(order(df_final$X.log10.Adjusted.p.value., decreasing = TRUE))
# 
# ggplot(df_final, aes(fill=condition, y=X.log10.Adjusted.p.value., x=term)) + 
#   geom_bar(position="stack", stat="identity")

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

# m <- data.frame(level = rep(unique(proportion_of_cells_each_niche$level), length(unique(proportion_of_cells_each_niche$niche))),
#                                              niche = rep(unique(proportion_of_cells_each_niche$niche), each=2),
#                                              m = NA)
# for (l in unique(proportion_of_cells_each_niche$level)) {
#   for (n in unique(proportion_of_cells_each_niche$niche)) {
#     m[m$level == l & m$niche ==n,3] <- mean(proportion_of_cells_each_niche[proportion_of_cells_each_niche$level ==l & proportion_of_cells_each_niche$niche ==n, 3])
#   }
# }

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

# check how many different niches are inside each fov and if there are differences between high and low

# num_niche_per_fov <- data.frame(fov = unique(metadata_niche$fov),
#                                 num_niche = NA)
# 
# for (f in num_niche_per_fov$fov) {
#   num_niche_per_fov[num_niche_per_fov$fov ==f, 2] <- length(unique(metadata_niche[metadata_niche$fov ==f,"novae_domains_20"]))
# }
# 
# num_niche_per_fov$level <- unlist(lapply(as.list(num_niche_per_fov$fov), function(x) {ifelse(x %in% c(21,25,29,30,37,38,39,40,41), 'low', 'high')}))
# 
# wilcox.test(num_niche_per_fov[num_niche_per_fov$level =='high',2],num_niche_per_fov[num_niche_per_fov$level =='low',2])
# 
# mean(num_niche_per_fov[num_niche_per_fov$level =='high',2])
# mean(num_niche_per_fov[num_niche_per_fov$level =='low',2])
# 
# # pathway score 
# pathways <- read.csv('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/pathways_novae_20domains_score_by_cells.csv')
# rownames(pathways) <- pathways[,1]
# pathways[,1] <- NULL
# pathways$fov <- unlist(lapply(as.list(rownames(pathways)), function(x) {metadata_niche[metadata_niche$X == x, 'fov']}))
# 
# niche_fov_pathway <- data.frame(pathway =rep(colnames(pathways)[1:34], length(unique(pathways$fov))*length(unique(pathways$novae_domains_20))),
#                                 niche = rep(unique(pathways$novae_domains_20), each = length(colnames(pathways)[1:34])),
#                                 fov = rep(unique(pathways$fov),each=length(colnames(pathways)[1:34])),
#                                 mean_score = NA)
# 
# 
# for (r in 1:nrow(niche_fov_pathway)) {
#   f <- niche_fov_pathway[r,'fov']
#   p <- niche_fov_pathway[r,'pathway']
#   n <- niche_fov_pathway[r,'niche']
#   niche_fov_pathway[r,4] <- mean(pathways[pathways$novae_domains_20 == n & pathways$fov == f, p])
# }
# write.table(niche_fov_pathway, file='/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/mean_pathway_score_by_niche_and_fov.txt', quote = F, sep = ',', row.names = F)
# 
# pathway_score_by_domain <- read.csv('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/pathways_novae_20domains.csv')
# rownames(pathway_score_by_domain) <- pathway_score_by_domain[,1]
# pathway_score_by_domain[,1] <- NULL
# 
# pathway_score_by_domain2 <- datNULLpathway_score_by_domain2 <- data.frame(niche = rep(pathway_score_by_domain$novae_domains_20, ncol(pathway_score_by_domain)-1),
#                                        pathway = rep(colnames(pathway_score_by_domain)[2:ncol(pathway_score_by_domain)], each = nrow(pathway_score_by_domain)), 
#                                        score = as.vector(as.matrix(pathway_score_by_domain[,2:ncol(pathway_score_by_domain)])))
# 
# pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/epithelial_high_low_mean.pdf', width = 5, height = 10)
# ggplot(pathway_score_by_domain2, aes(x=pathway , y=niche, size = score, fill = score)) +
#   geom_point(alpha=0.5) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.6), axis.title = element_blank())
# dev.off()
# 
# rowclus <- hclust(dist( pathway_score_by_domain ))    #cluster the rows
# colclus <- hclust(dist( t(pathway_score_by_domain) )) #cluster the columns
# 
# hm <- hmReady(pathway_score_by_domain, colclus=colclus,rowclus=rowclus)
# 
# hm$variable <- factor(hm$variable, levels = c('HALLMARK_MYOGENESIS','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','HALLMARK_MYC_TARGETS_V1',
#                                               'HALLMARK_PROTEIN_SECRETION','HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_XENOBIOTIC_METABOLISM',
#                                               'HALLMARK_UNFOLDED_PROTEIN_RESPONSE','HALLMARK_ANGIOGENESIS','HALLMARK_COAGULATION','HALLMARK_TGF_BETA_SIGNALING',
#                                               'HALLMARK_HEDGEHOG_SIGNALING','HALLMARK_KRAS_SIGNALING_UP','HALLMARK_KRAS_SIGNALING_DN','HALLMARK_DNA_REPAIR',
#                                               'HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','HALLMARK_IL6_JAK_STAT3_SIGNALING','HALLMARK_HEME_METABOLISM',
#                                               'HALLMARK_APICAL_JUNCTION','HALLMARK_APICAL_SURFACE','HALLMARK_NOTCH_SIGNALING','HALLMARK_INFLAMMATORY_RESPONSE',
#                                               'HALLMARK_WNT_BETA_CATENIN_SIGNALING','HALLMARK_IL2_STAT5_SIGNALING','HALLMARK_INTERFERON_ALPHA_RESPONSE',
#                                               'HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_COMPLEMENT','HALLMARK_PI3K_AKT_MTOR_SIGNALING','HALLMARK_GLYCOLYSIS',
#                                               'HALLMARK_MTORC1_SIGNALING','HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_P53_PATHWAY','HALLMARK_APOPTOSIS','HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY'))
# 
# hm$rowid <- factor(hm$rowid, levels =c("D890",  "D994", "D997" ,"D985","D999" ,"D935", "D989" , "D984" , "D987","D1003","D968", "D980" ,
#                                        "D998", "D1002","D967" ,"D959","D993","D1000" ,"D1001","D990"))
# 
# pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_seed/bubble_pathways_novae.pdf', height = 8, width = 20)
# ggplot() + 
#   geom_point(data=hm, aes(x=variable, y=rowid,fill=value), alpha=0.75, shape = 21, size = 5) +
#   scale_fill_gradient2(midpoint=0.05, low="black", mid="yellow",high="white") +
#   theme_classic() +
#   scale_size_continuous(range=c(1,6),limit=c(min(pathway_score_by_domain), max(pathway_score_by_domain)), breaks = seq(min(pathway_score_by_domain),max(pathway_score_by_domain), length.out=5)) +
#   theme(aspect.ratio = 0.45,panel.border = element_rect(linetype = "solid", fill = NA), axis.line=element_blank(), axis.text.x = element_text(angle=90, hjust = 1), axis.title = element_blank()) 
# dev.off()

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




################################### old


# DEG
new_column_values <- numeric(69526)
CD2 <- AddMetaData(object = CD2, metadata = new_column_values, 
                   col.name = "cluster_level")

CD2$cluster_level <- paste0(CD2$nb_clus,'_',CD2$Level)
a <- names(CD2@active.ident)
CD2@active.ident <- factor(unlist(lapply(as.list(a), function(x) {CD2$cluster_level})))
names(CD2@active.ident) <- a

######
from_david_to_cyrcle <- function(file_BP, degs_high, cell_type) {
  david_BP <- read.table(paste0('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/', file_BP, '.txt'), header=TRUE, sep='\t')
  
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
  write.table(david_BP2[,1], file=paste0('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/david_BP_', cell_type, '_modified_entire_GO.tsv'),col.names =TRUE, row.names =FALSE, quote = FALSE, sep='\t')
  
  david_BP2$ID <- gsub('~.*$', '', david_BP2$ID)
  
  write.table(david_BP2, file=paste0('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/david_BP_', cell_type, '_modified.tsv'),col.names =TRUE, row.names =FALSE, quote = FALSE, sep='\t')
  
  # file_KEGG <- read.table(paste0('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/', file_KEGG, '.txt'), header=TRUE, sep='\t')
  # 
  # file_KEGG$up_in_high <- 0
  # file_KEGG$up_in_low <- 0
  # 
  # for (r in 1:nrow(file_KEGG)) {
  #   high <- c()
  #   low <- c()
  #   for (x in unlist(strsplit(file_KEGG[r,'Genes'], ', '))) {
  #     if (x %in% degs_high) {
  #       high <- c(high, x)
  #     } else {
  #       low <- c(low, x)
  #     }
  #   } 
  #   file_KEGG[r,]$up_in_high <- length(high)
  #   file_KEGG[r,]$up_in_low <- length(low)
  # }
  # 
  # file_KEGG2 <- data.frame(ID = file_KEGG$Term, classification = file_KEGG$Category, FDR = file_KEGG$FDR, 
  #                          all = file_KEGG$Pop.Hits, enriched_in_HIGH = file_KEGG$up_in_high, enriched_in_LOW = file_KEGG$up_in_low)
  # file_KEGG2 <- file_KEGG2[order(file_KEGG2$FDR),]
  # 
  # if (nrow(file_KEGG2) > 10) {
  #   file_KEGG2 <- file_KEGG2[1:10,]
  # }
  # 
  # write.table(file_KEGG2[,1], file=paste0('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/david_KEGG_', cell_type,'_modified_GO.tsv'),col.names =TRUE, row.names =FALSE, quote = FALSE, sep='\t')
  # 
  # file_KEGG2$ID <- gsub(':.*$', '', file_KEGG2$ID)
  # 
  # write.table(file_KEGG2, file=paste0('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/david_KEGG_', cell_type,'_modified.tsv'),col.names =TRUE, row.names =FALSE, quote = FALSE, sep='\t')
}

###### B
deg_wilcoxon_B <- FindMarkers(CD2,test.use ='wilcox',ident.2 = 'B_LOW', ident.1 = 'B_HIGH')
deg_wilcoxon_B$diff <- "NO"
deg_wilcoxon_B$diff[deg_wilcoxon_B$avg_log2FC > 0 & deg_wilcoxon_B$p_val_adj<0.05] <- "B_HIGH"
deg_wilcoxon_B$diff[deg_wilcoxon_B$avg_log2FC < (-0) & deg_wilcoxon_B$p_val_adj<0.05] <- "B_LOW"
deg_wilcoxon_B$diff <- factor(deg_wilcoxon_B$diff)
write.table(deg_wilcoxon_B, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_B.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')
write.table(rownames(deg_wilcoxon_B[deg_wilcoxon_B$diff != 'NO',]), file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_B_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

deg_wilcoxon_B_fdr_0.05_high <- rownames(deg_wilcoxon_B[deg_wilcoxon_B$diff == "B_HIGH",])

from_david_to_cyrcle('david_B_BP', 'david_B_KEGG', deg_wilcoxon_B_fdr_0.05_high, 'B')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/volcano_differential_B_LOW_B_HIGH.pdf', height = 7, width =10)

ggplot2::ggplot(data=deg_wilcoxon_B, aes(x=avg_log2FC, y=-log10(p_val), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c(B_LOW="#28658f", NO="grey",B_HIGH="#ef857c")) +
  geom_point(data = subset(deg_wilcoxon_B, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.0001,rownames(deg_wilcoxon_B),'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('B')

dev.off()


####### B reg
deg_wilcoxon_Bmemory <- FindMarkers(CD2,test.use ='wilcox',ident.2 = 'B reg_LOW', ident.1 = 'B reg_HIGH')
deg_wilcoxon_Bmemory$diff <- "NO"
deg_wilcoxon_Bmemory$diff[deg_wilcoxon_Bmemory$avg_log2FC > 0 & deg_wilcoxon_Bmemory$p_val_adj<0.05] <- "B reg_HIGH"
deg_wilcoxon_Bmemory$diff[deg_wilcoxon_Bmemory$avg_log2FC < 0 & deg_wilcoxon_Bmemory$p_val_adj<0.05] <- "B reg_LOW"
deg_wilcoxon_Bmemory$diff <- factor(deg_wilcoxon_Bmemory$diff)
write.table(deg_wilcoxon_Bmemory, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_B_reg.txt',col.names =TRUE, row.names =TRUE, quote = FALSE, sep='\t')
write.table(rownames(deg_wilcoxon_Bmemory[deg_wilcoxon_Bmemory$diff != 'NO',]), file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_B_reg_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

deg_wilcoxon_Bmemory_fdr_0.05_high <- rownames(deg_wilcoxon_Bmemory[deg_wilcoxon_Bmemory$diff == "B memory/reg_HIGH",])

from_david_to_cyrcle('david_Bmem_BP', 'david_Bmem_KEGG', deg_wilcoxon_Bmemory_fdr_0.05_high, 'Bmem')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/volcano_differential_Bmemory_LOW_HIGH.pdf', height = 7, width =10)

ggplot2::ggplot(data=deg_wilcoxon_Bmemory, aes(x=avg_log2FC, y=-log10(p_val), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('B memory/reg_LOW'="#28658f", NO="grey",'B memory/reg_HIGH'="#ef857c")) +
  geom_point(data = subset(deg_wilcoxon_Bmemory, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.01,rownames(deg_wilcoxon_Bmemory),'')), size =4, max.overlaps = 100)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('B memory/reg')

dev.off()


###### endothelial
deg_wilcoxon_endo <- FindMarkers(CD2,test.use ='wilcox',ident.2 = 'Endothelial_LOW', ident.1 = 'Endothelial_HIGH')
deg_wilcoxon_endo$diff <- "NO"
deg_wilcoxon_endo$diff[deg_wilcoxon_endo$avg_log2FC > 0 & deg_wilcoxon_endo$p_val_adj<0.05] <- "Endothelial_HIGH"
deg_wilcoxon_endo$diff[deg_wilcoxon_endo$avg_log2FC < 0 & deg_wilcoxon_endo$p_val_adj<0.05] <- "Endothelial_LOW"
deg_wilcoxon_endo$diff <- factor(deg_wilcoxon_endo$diff)
write.table(deg_wilcoxon_endo, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_endothelial.txt',col.names =TRUE, row.names =TRUE, quote = FALSE, sep='\t')
write.table(rownames(deg_wilcoxon_endo[deg_wilcoxon_endo$diff != 'NO',]), file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_endothelial_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

deg_wilcoxon_endo_fdr_0.05_high <- rownames(deg_wilcoxon_endo[deg_wilcoxon_endo$diff == "Endothelial_HIGH" ,])

from_david_to_cyrcle('david_endo_BP', 'david_endo_KEGG', deg_wilcoxon_endo_fdr_0.05_high, 'endo')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/volcano_differential_endo_LOW_HIGH.pdf', height = 7, width =10)

ggplot2::ggplot(data=deg_wilcoxon_endo, aes(x=avg_log2FC, y=-log10(p_val), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('Endothelial_LOW'="#28658f", NO="grey",'Endothelial_HIGH'="#ef857c")) +
  geom_point(data = subset(deg_wilcoxon_endo, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.01,rownames(deg_wilcoxon_endo),'')), size =4,  max.overlaps = 100)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('Endothelial')

dev.off()

###### glia
deg_wilcoxon_glia <- FindMarkers(CD2,test.use ='wilcox',ident.2 = 'Enteric glia cell_LOW', ident.1 = 'Enteric glia cell_HIGH')
deg_wilcoxon_glia$diff <- "NO"
deg_wilcoxon_glia$diff[deg_wilcoxon_glia$avg_log2FC > 0 & deg_wilcoxon_glia$p_val_adj<0.05] <- "Enteric glia cell_HIGH"
deg_wilcoxon_glia$diff[deg_wilcoxon_glia$avg_log2FC < 0 & deg_wilcoxon_glia$p_val_adj<0.05] <- "Enteric glia cell_LOW"
deg_wilcoxon_glia$diff <- factor(deg_wilcoxon_glia$diff)
write.table(deg_wilcoxon_glia, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_glia.txt',col.names =TRUE, row.names =TRUE, quote = FALSE, sep='\t')

deg_wilcoxon_glia_fdr_0.05 <- rownames(deg_wilcoxon_glia[deg_wilcoxon_glia$p_val_adj<0.05,])
# gsea not possible becouse there is only one gene with fdr < 0.05

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/volcano_differential_glia_LOW_HIGH.pdf', height = 7, width =10)

ggplot2::ggplot(data=deg_wilcoxon_glia, aes(x=avg_log2FC, y=-log10(p_val), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('Enteric glia cell_LOW'="#28658f", NO="grey",'Enteric glia cell_HIGH'="#ef857c")) +
  geom_point(data = subset(deg_wilcoxon_glia, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.01 ,rownames(deg_wilcoxon_glia),'')), size =4)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('Enteric glia cells')

dev.off()


####### epithelial
deg_wilcoxon_epithelial <- FindMarkers(CD2,test.use ='wilcox',ident.2 = 'Epithelial_LOW', ident.1 = 'Epithelial_HIGH')
deg_wilcoxon_epithelial$diff <- "NO"
deg_wilcoxon_epithelial$diff[deg_wilcoxon_epithelial$avg_log2FC > 0 & deg_wilcoxon_epithelial$p_val_adj<0.05] <- "Epithelial_HIGH"
deg_wilcoxon_epithelial$diff[deg_wilcoxon_epithelial$avg_log2FC < 0 & deg_wilcoxon_epithelial$p_val_adj<0.05] <- "Epithelial_LOW"
deg_wilcoxon_epithelial$diff <- factor(deg_wilcoxon_epithelial$diff)
write.table(deg_wilcoxon_epithelial, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_epithelial.txt',col.names =TRUE, row.names =TRUE, quote = FALSE, sep='\t')
write.table(rownames(deg_wilcoxon_epithelial[deg_wilcoxon_epithelial$diff != 'NO',]), file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_epithelial_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

deg_wilcoxon_epithelial_fdr_0.05_high <- rownames(deg_wilcoxon_epithelial[deg_wilcoxon_epithelial$diff == "Epithelial_HIGH",])

from_david_to_cyrcle('david_epithelial_BP', deg_wilcoxon_epithelial_fdr_0.05_high, 'epithelial')

sign_genes <- read_excel(path = '/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs epiteliali volcano plot Elena.xlsx')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/volcano_differential_epithelial_LOW_HIGH.pdf', height =7, width =10)

ggplot2::ggplot(data=deg_wilcoxon_epithelial, aes(x=avg_log2FC, y=-log10(p_val), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('Epithelial_LOW'="#28658f", NO="grey",'Epithelial_HIGH'="#ef857c")) +
  geom_point(data = subset(deg_wilcoxon_epithelial, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(rownames(deg_wilcoxon_epithelial) %in% sign_genes$...1, rownames(deg_wilcoxon_epithelial),'')), size =4, max.overlaps =1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('Epithelial')

dev.off()

### check
degs <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_epithelial.txt', sep='\t')

david_BP <- read.table(paste0('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/david_epithelial_BP.txt'), header=TRUE, sep='\t')

enriched_go <- unlist(strsplit(david_BP[david_BP$Term =='GO:0006955~immune response','Genes'], ', '))
length(enriched_go[which(enriched_go %in% rownames(degs[degs$diff == 'Epithelial_LOW',]))])

length(enriched_go[which(enriched_go %in% rownames(degs[degs$diff == 'Epithelial_HIGH',]))])



#### fibroblast
deg_wilcoxon_fibro <- FindMarkers(CD2,test.use ='wilcox',ident.2 = 'Fibroblast_LOW', ident.1 = 'Fibroblast_HIGH')
deg_wilcoxon_fibro$diff <- "NO"
deg_wilcoxon_fibro$diff[deg_wilcoxon_fibro$avg_log2FC > 0 & deg_wilcoxon_fibro$p_val_adj<0.05] <- "Fibroblast_HIGH"
deg_wilcoxon_fibro$diff[deg_wilcoxon_fibro$avg_log2FC < 0 & deg_wilcoxon_fibro$p_val_adj<0.05] <- "Fibroblast_LOW"
deg_wilcoxon_fibro$diff <- factor(deg_wilcoxon_fibro$diff)
write.table(deg_wilcoxon_fibro, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_fibroblast.txt',col.names =TRUE, row.names =TRUE, quote = FALSE, sep='\t')
write.table(rownames(deg_wilcoxon_fibro[deg_wilcoxon_fibro$diff != 'NO',]), file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_fibroblast_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

deg_wilcoxon_fibro_fdr_0.05_high <- rownames(deg_wilcoxon_fibro[deg_wilcoxon_fibro$diff == "Fibroblast_HIGH",])

from_david_to_cyrcle('david_fibro_BP', deg_wilcoxon_fibro_fdr_0.05_high, 'fibroblast')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/volcano_differential_fibro_LOW_HIGH.pdf', height =7, width =10)

ggplot2::ggplot(data=deg_wilcoxon_fibro, aes(x=avg_log2FC, y=-log10(p_val), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('Fibroblast_LOW'="#28658f", NO="grey",'Fibroblast_HIGH'="#ef857c")) +
  geom_point(data = subset(deg_wilcoxon_fibro, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.01,rownames(deg_wilcoxon_fibro),'')), size =4, max.overlaps =40)  +
 geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('Fibroblast')

dev.off()


### macrophage
deg_wilcoxon_macrophage <- FindMarkers(CD2,test.use ='wilcox',ident.2 = 'Macrophage_LOW', ident.1 = 'Macrophage_HIGH')
deg_wilcoxon_macrophage$diff <- "NO"
deg_wilcoxon_macrophage$diff[deg_wilcoxon_macrophage$avg_log2FC > 0 & deg_wilcoxon_macrophage$p_val_adj<0.05] <- "Macrophage_HIGH"
deg_wilcoxon_macrophage$diff[deg_wilcoxon_macrophage$avg_log2FC < 0 & deg_wilcoxon_macrophage$p_val_adj<0.05] <- "Macrophage_LOW"
deg_wilcoxon_macrophage$diff <- factor(deg_wilcoxon_macrophage$diff)
write.table(deg_wilcoxon_macrophage, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_macrophage.txt',col.names =TRUE, row.names =TRUE, quote = FALSE, sep='\t')
write.table(rownames(deg_wilcoxon_macrophage[deg_wilcoxon_macrophage$diff != 'NO',]), file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_macrophage_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

deg_wilcoxon_macrophage_fdr_0.05_high <- rownames(deg_wilcoxon_macrophage[deg_wilcoxon_macrophage$diff == "Macrophage_HIGH",])

from_david_to_cyrcle('david_macrophage_BP', 'david_macrophage_KEGG', deg_wilcoxon_macrophage_fdr_0.05_high, 'macrophage')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/volcano_differential_Macrophage_LOW_HIGH.pdf', height =7, width =10)

ggplot2::ggplot(data=deg_wilcoxon_macrophage, aes(x=avg_log2FC, y=-log10(p_val), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('Macrophage_LOW'="#28658f", NO="grey",'Macrophage_HIGH'="#ef857c")) +
  geom_point(data = subset(deg_wilcoxon_macrophage, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.01 ,rownames(deg_wilcoxon_macrophage),'')), size =4, max.overlaps =40)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('Macrophage')

dev.off()


## mast cells
deg_wilcoxon_mast <- FindMarkers(CD2,test.use ='wilcox',ident.2 = 'Mast_LOW', ident.1 = 'Mast_HIGH')
deg_wilcoxon_mast$diff <- "NO"
deg_wilcoxon_mast$diff[deg_wilcoxon_mast$avg_log2FC > 0 & deg_wilcoxon_mast$p_val_adj<0.05] <- "Mast_HIGH"
deg_wilcoxon_mast$diff[deg_wilcoxon_mast$avg_log2FC < 0 & deg_wilcoxon_mast$p_val_adj<0.05] <- "Mast_LOW"
deg_wilcoxon_mast$diff <- factor(deg_wilcoxon_mast$diff)
write.table(deg_wilcoxon_mast, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_mast.txt',col.names =TRUE, row.names =TRUE, quote = FALSE, sep='\t')
write.table(rownames(deg_wilcoxon_mast[deg_wilcoxon_mast$diff != 'NO',]), file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_mast_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

deg_wilcoxon_mast_fdr_0.05_high <- rownames(deg_wilcoxon_mast[deg_wilcoxon_mast$diff == "Mast cell_HIGH",])

from_david_to_cyrcle('david_mast_BP', 'david_mast_KEGG', deg_wilcoxon_mast_fdr_0.05_high, 'mast')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/volcano_differential_Mast_LOW_HIGH.pdf', height =7, width =10)

ggplot2::ggplot(data=deg_wilcoxon_mast, aes(x=avg_log2FC, y=-log10(p_val), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('Mast cell_LOW'="#28658f", NO="grey",'Mast cell_HIGH'="#ef857c")) +
  geom_point(data = subset(deg_wilcoxon_mast, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.01,rownames(deg_wilcoxon_mast),'')), size =4, max.overlaps =40)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('Mast cells')

dev.off()


## plasma IgA cells
deg_wilcoxon_plasma <- FindMarkers(CD2,test.use ='wilcox',ident.2 = 'Plasma IgA_LOW', ident.1 = 'Plasma IgA_HIGH')
deg_wilcoxon_plasma$diff <- "NO"
deg_wilcoxon_plasma$diff[deg_wilcoxon_plasma$avg_log2FC > 0 & deg_wilcoxon_plasma$p_val_adj<0.05] <- "Plasma IgA_HIGH"
deg_wilcoxon_plasma$diff[deg_wilcoxon_plasma$avg_log2FC < 0 & deg_wilcoxon_plasma$p_val_adj<0.05] <- "Plasma IgA_LOW"
deg_wilcoxon_plasma$diff <- factor(deg_wilcoxon_plasma$diff)
write.table(deg_wilcoxon_plasma, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_plasma_IgA.txt',col.names =TRUE, row.names =TRUE, quote = FALSE, sep='\t')
write.table(rownames(deg_wilcoxon_plasma[deg_wilcoxon_plasma$diff != 'NO',]), file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_plasma_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

deg_wilcoxon_plasma_fdr_0.05_high <- rownames(deg_wilcoxon_plasma[deg_wilcoxon_plasma$diff == "Plasma_HIGH",])

from_david_to_cyrcle('david_plasma_BP', 'david_plasma_KEGG', deg_wilcoxon_plasma_fdr_0.05_high, 'plasma')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/volcano_differential_Plasma_LOW_HIGH.pdf', height =7, width =10)

ggplot2::ggplot(data=deg_wilcoxon_plasma, aes(x=avg_log2FC, y=-log10(p_val), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('Plasma_LOW'="#28658f", NO="grey",'Plasma_HIGH'="#ef857c")) +
  geom_point(data = subset(deg_wilcoxon_plasma, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.01,rownames(deg_wilcoxon_plasma),'')), size =4, max.overlaps =100)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('Plasma')

dev.off()

## plasma IgG cells
deg_wilcoxon_plasma <- FindMarkers(CD2,test.use ='wilcox',ident.2 = 'Plasma IgG_LOW', ident.1 = 'Plasma IgG_HIGH')
deg_wilcoxon_plasma$diff <- "NO"
deg_wilcoxon_plasma$diff[deg_wilcoxon_plasma$avg_log2FC > 0 & deg_wilcoxon_plasma$p_val_adj<0.05] <- "Plasma IgG_HIGH"
deg_wilcoxon_plasma$diff[deg_wilcoxon_plasma$avg_log2FC < 0 & deg_wilcoxon_plasma$p_val_adj<0.05] <- "Plasma IgG_LOW"
deg_wilcoxon_plasma$diff <- factor(deg_wilcoxon_plasma$diff)
write.table(deg_wilcoxon_plasma, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_plasma_IgG.txt',col.names =TRUE, row.names =TRUE, quote = FALSE, sep='\t')
write.table(rownames(deg_wilcoxon_plasma[deg_wilcoxon_plasma$diff != 'NO',]), file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_plasma_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

deg_wilcoxon_plasma_fdr_0.05_high <- rownames(deg_wilcoxon_plasma[deg_wilcoxon_plasma$diff == "Plasma_HIGH",])

from_david_to_cyrcle('david_plasma_BP', 'david_plasma_KEGG', deg_wilcoxon_plasma_fdr_0.05_high, 'plasma')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/volcano_differential_Plasma_LOW_HIGH.pdf', height =7, width =10)

ggplot2::ggplot(data=deg_wilcoxon_plasma, aes(x=avg_log2FC, y=-log10(p_val), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('Plasma_LOW'="#28658f", NO="grey",'Plasma_HIGH'="#ef857c")) +
  geom_point(data = subset(deg_wilcoxon_plasma, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.01,rownames(deg_wilcoxon_plasma),'')), size =4, max.overlaps =100)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('Plasma')

dev.off()


## smooth muscle cells
deg_wilcoxon_smooth <- FindMarkers(CD2,test.use ='wilcox',ident.2 = 'Smooth muscle cell_LOW', ident.1 = 'Smooth muscle cell_HIGH')
deg_wilcoxon_smooth$diff <- "NO"
deg_wilcoxon_smooth$diff[deg_wilcoxon_smooth$avg_log2FC > 0 & deg_wilcoxon_smooth$p_val_adj<0.05] <- "Smooth muscle cell_HIGH"
deg_wilcoxon_smooth$diff[deg_wilcoxon_smooth$avg_log2FC < 0 & deg_wilcoxon_smooth$p_val_adj<0.05] <- "Smooth muscle cell_LOW"
deg_wilcoxon_smooth$diff <- factor(deg_wilcoxon_smooth$diff)
write.table(deg_wilcoxon_smooth, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_smooth_muscle.txt',col.names =TRUE, row.names =TRUE, quote = FALSE, sep='\t')
write.table(rownames(deg_wilcoxon_smooth[deg_wilcoxon_smooth$diff != 'NO',]), file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_smooth_muscle_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

deg_wilcoxon_smooth_fdr_0.05_high <- rownames(deg_wilcoxon_smooth[deg_wilcoxon_smooth$diff == "Smooth muscle_HIGH",])

from_david_to_cyrcle('david_smooth_BP', 'david_smooth_KEGG', deg_wilcoxon_smooth_fdr_0.05_high, 'smooth')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/volcano_differential_Smooth_muscle_LOW_HIGH.pdf', height =7, width =10)

ggplot2::ggplot(data=deg_wilcoxon_smooth, aes(x=avg_log2FC, y=-log10(p_val), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('Smooth muscle_LOW'="#28658f", NO="grey",'Smooth muscle_HIGH'="#ef857c")) +
  geom_point(data = subset(deg_wilcoxon_smooth, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.01,rownames(deg_wilcoxon_smooth),'')), size =4, max.overlaps =40)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('Smooth muscle')

dev.off()


### T cells
deg_wilcoxon_T <- FindMarkers(CD2,test.use ='wilcox',ident.2 = 'T_LOW', ident.1 = 'T_HIGH')
deg_wilcoxon_T$diff <- "NO"
deg_wilcoxon_T$diff[deg_wilcoxon_T$avg_log2FC > 0 & deg_wilcoxon_T$p_val_adj<0.05] <- "T_HIGH"
deg_wilcoxon_T$diff[deg_wilcoxon_T$avg_log2FC < 0 & deg_wilcoxon_T$p_val_adj<0.05] <- "T_LOW"
deg_wilcoxon_T$diff <- factor(deg_wilcoxon_T$diff)
write.table(deg_wilcoxon_T, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_T.txt',col.names =TRUE, row.names =TRUE, quote = FALSE, sep='\t')
write.table(rownames(deg_wilcoxon_T[deg_wilcoxon_T$diff != 'NO',]), file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_T_names.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

deg_wilcoxon_T_fdr_0.05_high <- rownames(deg_wilcoxon_T[deg_wilcoxon_T$diff == "T_HIGH",])

from_david_to_cyrcle('david_T_BP', deg_wilcoxon_T_fdr_0.05_high, 'T')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/volcano_differential_T_LOW_HIGH.pdf', height =7, width =10)

ggplot2::ggplot(data=deg_wilcoxon_T, aes(x=avg_log2FC, y=-log10(p_val), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('T_LOW'="#28658f", NO="grey",'T_HIGH'="#ef857c")) +
  geom_point(data = subset(deg_wilcoxon_T, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.01,rownames(deg_wilcoxon_T),'')), size =4, max.overlaps =50)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('T cells')

dev.off()


### Undifferentiated cells
deg_wilcoxon_undif <- FindMarkers(CD2,test.use ='wilcox',ident.2 = 'Undifferentiated_LOW', ident.1 = 'Undifferentiated_HIGH')
deg_wilcoxon_undif$diff <- "NO"
deg_wilcoxon_undif$diff[deg_wilcoxon_undif$avg_log2FC > 0 & deg_wilcoxon_undif$p_val_adj<0.05] <- "Undifferentiated_HIGH"
deg_wilcoxon_undif$diff[deg_wilcoxon_undif$avg_log2FC < 0 & deg_wilcoxon_undif$p_val_adj<0.05] <- "Undifferentiated_LOW"
deg_wilcoxon_undif$diff <- factor(deg_wilcoxon_undif$diff)
write.table(deg_wilcoxon_undif, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_Undifferentiated.txt',col.names =TRUE, row.names =TRUE, quote = FALSE, sep='\t')


### distances

### distances between epithelial and each other cell type without considering the fovs separately
all_distances_epi <- matrix(ncol=3)
colnames(all_distances_epi) <- c('iNKT10 level', 'Distance','Cell_type')
for (f in list.files('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_epi', full.names = TRUE)) {
  dist <- data.frame(readxl::read_excel(f))
  colnames(dist) <- c('iNKT10 level', 'Distance')
  dist$Cell_type <- gsub('.xlsx', '',strsplit(f, split='/')[[1]][12]) 
  all_distances_epi <- rbind(all_distances_epi, dist)
}
all_distances_epi <-  all_distances_epi[-1,] 
colnames(all_distances_epi) <- c('iNKT10_level', 'Distance','Cell_type')

stats_epi <- all_distances_epi %>%
  group_by(Cell_type) %>%
  t_test(Distance~iNKT10_level) %>%
  adjust_pvalue(method = 'bonferroni')

stats_epi$y.position <- unlist(lapply(seq_along(1:length(stats_epi$Cell_type)), function(x) {
  max(all_distances_epi[all_distances_epi$Cell_type == stats_epi$Cell_type[x],'Distance'][! is.na(all_distances_epi[all_distances_epi$Cell_type == stats_epi$Cell_type[x],'Distance'])]) +max(all_distances_epi[all_distances_epi$Cell_type == stats_epi$Cell_type[x],'Distance'][! is.na(all_distances_epi[all_distances_epi$Cell_type == stats_epi$Cell_type[x],'Distance'])])/14
}))

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/boxplot_epi_distances.pdf', width = 30, height = 20)
ggplot(all_distances_epi, aes(x= iNKT10_level, y=Distance,color=iNKT10_level)) +
  geom_boxplot(aes(color=iNKT10_level),width=0.4, lwd=1, outliers = FALSE) +
  scale_colour_manual(values = c("#ef857c","#28658f")) +
  facet_wrap(~Cell_type,nrow = 3, ncol = 5,scale = "free_y") +
  theme_classic() +
  scale_shape(solid=T)  +
  ylab("") +
  xlab("") +
  theme(aspect.ratio=1.5, strip.text=element_text(size = 25, face = 'bold'), axis.text.x = element_text(size=25, angle = 45,hjust = 1), 
        axis.text.y = element_text(size=25),
        legend.text = element_text(size = 25),legend.title =  element_text(size = 25), strip.text.x = element_text(size=20)) +
  stat_pvalue_manual(stats_epi, label='p.adj', size=7)
dev.off()

### distances between epithelial and each other cell type considering the means computed for each fov
all_distances_epi_mean <- matrix(ncol=3)
colnames(all_distances_epi_mean) <- c('iNKT10_level', 'Distance','Cell_type')
for (f in list.files('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_epi_mean', full.names = TRUE)) {
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

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/boxplot_epi_distances_by_fov_points.pdf', width = 30, height = 20)
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


## bubble plot with mean distance between epi and each other cell type without considering different fovs
distances_low <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_low.txt', header=TRUE, sep=',')
distances_high <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_high.txt', header=TRUE, sep=',')

distances_high_low <- cbind(distances_low, mean_high = NA)

for (r in 1:nrow(distances_low)) {
  if (length(distances_high[(distances_high$cell_type1 == distances_low[r,1]) & (distances_high$cell_type2 == distances_low[r,2]),'mean']) != 0) {
      distances_high_low[(distances_high_low$cell_type1 == distances_low[r,1]) & (distances_high_low$cell_type2 == distances_low[r,2]), 'mean_high'] <- distances_high[(distances_high$cell_type1 == distances_low[r,1]) & (distances_high$cell_type2 == distances_low[r,2]), 'mean']
    } else if (length(distances_high[(distances_high$cell_type2 == distances_low[r,1]) & (distances_high$cell_type1 == distances_low[r,2]), 'mean']) != 0) {
      distances_high_low[(distances_high_low$cell_type1 == distances_low[r,1]) & (distances_high_low$cell_type2 == distances_low[r,2]), 'mean_high'] <- distances_high[(distances_high$cell_type2 == distances_low[r,1]) & (distances_high$cell_type1 == distances_low[r,2]), 'mean']
    } 
}

distances_high_low_epi <- distances_high_low[(distances_high_low$cell_type1 == 'Epithelial') |(distances_high_low$cell_type2 == 'Epithelial') ,]

lapply(as.list(1:nrow(distances_high_low_epi)), function(x) {if (distances_high_low_epi[x,2] != 'Epithelial') {
  distances_high_low_epi[x,1] <<- distances_high_low_epi[x,2]
distances_high_low_epi[x,2] <<- 'Epithelial'}})

distances_df <- data.frame(cell = rep(distances_high_low_epi$cell_type1, 2), condition = rep(c('High iNKT10', 'Low iNKT10'), each = nrow(distances_high_low_epi)),
                           value = c(distances_high_low_epi$mean_high, distances_high_low_epi$mean))

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/epithelial_high_low.pdf', width = 5, height = 10)
ggplot(distances_df, aes(x=condition , y=cell, size = -value)) +
  geom_point(alpha=0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6, size = 15), axis.title = element_blank(), axis.text.y = element_text(size=15))
dev.off()

colnames(distances_high_low_epi) <- c('Cell type1', 'Cell type2', 'Mean distance in low iNKT10', 'Mean distance in high iNKT10')

## bubble plot with mean distance between the mean distances between epi and each other cell type in the different fovs
distances_low <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_low_mean.txt', header=TRUE, sep=',')
distances_high <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_high_mean.txt', header=TRUE, sep=',')

distances_high_low <- cbind(distances_low, mean_high = NA)

for (r in 1:nrow(distances_low)) {
  if (length(distances_high[(distances_high$cell_type1 == distances_low[r,1]) & (distances_high$cell_type2 == distances_low[r,2]),'mean']) != 0) {
    distances_high_low[(distances_high_low$cell_type1 == distances_low[r,1]) & (distances_high_low$cell_type2 == distances_low[r,2]), 'mean_high'] <- distances_high[(distances_high$cell_type1 == distances_low[r,1]) & (distances_high$cell_type2 == distances_low[r,2]), 'mean']
  } else if (length(distances_high[(distances_high$cell_type2 == distances_low[r,1]) & (distances_high$cell_type1 == distances_low[r,2]), 'mean']) != 0) {
    distances_high_low[(distances_high_low$cell_type1 == distances_low[r,1]) & (distances_high_low$cell_type2 == distances_low[r,2]), 'mean_high'] <- distances_high[(distances_high$cell_type2 == distances_low[r,1]) & (distances_high$cell_type1 == distances_low[r,2]), 'mean']
  } 
}

distances_high_low_epi <- distances_high_low[(distances_high_low$cell_type1 == 'Epithelial') |(distances_high_low$cell_type2 == 'Epithelial') ,]

lapply(as.list(1:nrow(distances_high_low_epi)), function(x) {if (distances_high_low_epi[x,2] != 'Epithelial') {distances_high_low_epi[x,1] <<- distances_high_low_epi[x,2]
distances_high_low_epi[x,2] <<- 'Epithelial'}})

distances_df <- data.frame(cell = rep(distances_high_low_epi$cell_type1, 2), condition = rep(c('High iNKT10', 'Low iNKT10'), each = nrow(distances_high_low_epi)),
                           value = c(distances_high_low_epi$mean_high, distances_high_low_epi$mean))

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/epithelial_high_low_mean.pdf', width = 5, height = 10)
ggplot(distances_df, aes(x=condition , y=cell, size = -value)) +
  geom_point(alpha=0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6, size = 15), axis.title = element_blank(), axis.text.y = element_text(size=15))
dev.off()

colnames(distances_high_low_epi) <- c('Cell type1', 'Cell type2', 'Mean distance in low iNKT10', 'Mean distance in high iNKT10')

### compute the % of ZBTB16+ T cells and of ZBTB16+IL10+ T cells over the total amount of cells per fov
num_cells_each_fov <- c(491,1079,2464,422,690,871,1240,2005,683,549,2424,855,584,711,653,853,3556,1458,2450)
names(num_cells_each_fov) <- c('21','23','24','25','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41')

fov_low_ZBTB16 = c(2, 5, 8, 156, 22, 0, 216, 9, 114)
names(fov_low_ZBTB16) <- c('21','25','29','30','37','38', '39','40','41')
fov_low_ZBTB16_IL10 = c(0, 0, 1, 6, 2, 0, 7, 0, 4)
names(fov_low_ZBTB16_IL10) <- c('21','25','29','30','37','38', '39','40','41')

fov_high_ZBTB16 = c(15,36,9,3,8,3,234,18,2,7)
names(fov_high_ZBTB16) <- c('23','24','27','28','31','32','33','34','35','36')
fov_high_ZBTB16_IL10 = c(0, 0, 1,0,0,0,6,0,0,0)
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

pdf('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/boxplot_T_ZBTB16_IL10_per_fov.pdf')
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
  if (length(readBin(file(paste0('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/number_epithelial_less_', as.character(distance), '_ZBTB16_IL10'),open="rb"),"raw", 65536)) > 0) {
    epi_less_ZBTB16_IL10 <- read.table(paste0('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/number_epithelial_less_', as.character(distance), '_ZBTB16_IL10'))
  }
  if (length(readBin(file(paste0('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/number_fibroblast_less_', as.character(distance), '_ZBTB16_IL10'),open="rb"),"raw", 65536)) > 0) {
    fibro_less_ZBTB16_IL10 <- read.table(paste0('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/number_fibroblast_less_', as.character(distance), '_ZBTB16_IL10'))
  }
  epi_less_ZBTB16 <- read.table(paste0('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/number_epithelial_less_', as.character(distance),'_ZBTB16'))
  fibro_less_ZBTB16 <- read.table(paste0('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/number_fibroblast_less_', as.character(distance),'_ZBTB16'))
  epi_less_IL10 <- read.table(paste0('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/number_epithelial_less_', as.character(distance),'_IL10'))
  fibro_less_IL10 <- read.table(paste0('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/number_fibroblast_less_', as.character(distance),'_IL10'))
  
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
  
  ggsave(paste0('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/number_epi_cells_',as.character(distance), '_T.pdf'),p, height = 5, width =3.5)
  
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
  
  ggsave(paste0('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/number_fibro_cells_', as.character(distance), '_T.pdf'),p, height = 5, width =3.5)
  
}

percentage_cells_within(100)
percentage_cells_within(500)

##### identify DEGs between epithelial/fibro closer/farther than 500 micro from T ZBTB16+IL10+ cells considering high and low cells together
##### identify DEGs between epithelial/fibro closer than 500 micro from T ZBTB16+IL10+ cells and closer than 500 micro from T ZBTB16+ cells considering high and low cells together
##### identify DEGs between epithelial/fibro closer than 500 micro from T ZBTB16+IL10+ cells and closer than 500 micro from T IL10+ cells considering high and low cells together


### DEGs between epithelial cells within 500 micro from T ZBTB16+IL10 and epithelial cells more far than 500 micro from T ZBTB16+IL10 
epi_less_500 <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/epithelial_less_500_T_ZBTB16_IL10')
epi_more_500 <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/epithelial_more_500_T_ZBTB16_IL10')

counts <- CD2@assays$Nanostring$counts[-which(grepl('NegPrb', rownames(CD2@assays$Nanostring$counts))),]
CD2 <- subset(CD2, features = rownames(counts))
                     
CD_epithelial_500 <- subset(CD2, subset= cell %in% c(unique(epi_less_500$V1), unique(epi_more_500$V1)))

a <- names(CD_epithelial_500@active.ident)
CD_epithelial_500@active.ident <- factor(unlist(lapply(as.list(CD_epithelial_500@meta.data$cell), function(x) {
  if (x %in% unique(epi_less_500$V1)) {
    '<500'
  } else {
    '>500'
  }
})))
names(CD_epithelial_500@active.ident) <- a

deg_farther_closest <- FindMarkers(CD_epithelial_500,test.use ='wilcox',ident.2 = '<500', ident.1 = '>500')
deg_farther_closest$diff <- "NO"
deg_farther_closest$diff[deg_farther_closest$avg_log2FC > 0 & deg_farther_closest$p_val_adj<0.05] <- ">500"
deg_farther_closest$diff[deg_farther_closest$avg_log2FC < 0 & deg_farther_closest$p_val_adj<0.05] <- "<500"
deg_farther_closest$diff <- factor(deg_farther_closest$diff)
write.table(deg_farther_closest, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_epithelial_less_more_500_T_ZBTB16_IL10.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/volcano_differential_epithelial_less_more_500_T_ZBTB16_IL10.pdf', height = 7, width =10)

ggplot2::ggplot(data=deg_farther_closest, aes(x=avg_log2FC, y=-log10(p_val), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('<500'="#28658f", '>500'="#ef857c")) +
  geom_point(data = subset(deg_farther_closest, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.05,rownames(deg_farther_closest),'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs between epithelial cells closer and farther than 500 µm from T ZBTB16+IL10+ cells')

dev.off()


### DEGs between fibroblast cells within 500 micro from T ZBTB16+IL10 and fibro more far than 500 micro from T ZBTB16+IL10 
fibro_less_500 <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/fibroblast_less_500_T_ZBTB16_IL10')
fibro_more_500 <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/fibroblast_more_500_T_ZBTB16_IL10')

CD_fibroblast_500 <- subset(CD2, subset= cell %in% c(unique(fibro_less_500$V1), unique(fibro_more_500$V1)))

a <- names(CD_fibroblast_500@active.ident)
CD_fibroblast_500@active.ident <- factor(unlist(lapply(as.list(CD_fibroblast_500@meta.data$cell), function(x) {
  if (x %in% unique(fibro_less_500$V1)) {
    '<500'
  } else {
    '>500'
  }
})))
names(CD_fibroblast_500@active.ident) <- a

deg_farther_closest <- FindMarkers(CD_fibroblast_500,test.use ='wilcox',ident.2 = '<500', ident.1 = '>500')
deg_farther_closest$diff <- "NO"
deg_farther_closest$diff[deg_farther_closest$avg_log2FC > 0 & deg_farther_closest$p_val_adj<0.05] <- ">500"
deg_farther_closest$diff[deg_farther_closest$avg_log2FC < 0 & deg_farther_closest$p_val_adj<0.05] <- "<500"
deg_farther_closest$diff <- factor(deg_farther_closest$diff)
write.table(deg_farther_closest, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_fibroblast_less_more_500_T_ZBTB16_IL10.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/volcano_differential_fibroblast_less_more_500_T_ZBTB16_IL10.pdf', height = 7, width =10)

ggplot2::ggplot(data=deg_farther_closest, aes(x=avg_log2FC, y=-log10(p_val), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('<500'="#28658f", '>500'="#ef857c")) +
  geom_point(data = subset(deg_farther_closest, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.05,rownames(deg_farther_closest),'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs between fibroblasts closer and farther than 500 µm from T ZBTB16+IL10+ cells')

dev.off()


### DEGs between epithelial cells within 500 micro from T ZBTB16+IL10 and epithelial cells within 500 micro from T ZBTB16
epi_less_500_ZBTB16_IL10 <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/epithelial_less_500_T_ZBTB16_IL10')
epi_less_500_ZBTB16 <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/epithelial_less_500_T_ZBTB16')

CD_epithelial_comparison2 <- subset(CD2, subset= cell %in% c(unique(epi_less_500_ZBTB16_IL10$V1), unique(epi_less_500_ZBTB16$V1)))

a <- names(CD_epithelial_comparison2@active.ident)
CD_epithelial_comparison2@active.ident <- factor(unlist(lapply(as.list(CD_epithelial_comparison2@meta.data$cell), function(x) {
  if (x %in% epi_less_500_ZBTB16_IL10$V1) {
    'ZBTB16_IL10'
  } else {
    'ZBTB16'
  }
})))
names(CD_epithelial_comparison2@active.ident) <- a

deg_farther_closest <- FindMarkers(CD_epithelial_comparison2,test.use ='wilcox',ident.2 = 'ZBTB16', ident.1 = 'ZBTB16_IL10')
deg_farther_closest$diff <- "NO"
deg_farther_closest$diff[deg_farther_closest$avg_log2FC > 0 & deg_farther_closest$p_val_adj<0.05] <- "ZBTB16_IL10"
deg_farther_closest$diff[deg_farther_closest$avg_log2FC < 0 & deg_farther_closest$p_val_adj<0.05] <- "ZBTB16"
deg_farther_closest$diff <- factor(deg_farther_closest$diff)
write.table(deg_farther_closest, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_epithelial_less_500_T_ZBTB16_IL10_less_500_T_ZBTB16.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/volcano_differential_epithelial_less_500_T_ZBTB16_IL10_less_500_T_ZBTB10.pdf', height = 7, width =10)

ggplot2::ggplot(data=deg_farther_closest, aes(x=avg_log2FC, y=-log10(p_val), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('ZBTB16'="#28658f", 'ZBTB16_IL10'="#ef857c")) +
  geom_point(data = subset(deg_farther_closest, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.05,rownames(deg_farther_closest),'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs between epithelial cells within 500 µm from T ZBTB16+IL10+ cells and\nepithelial cells within 500 µm from T ZBTB16+ cells')

dev.off()


### DEGs between fibroblast cells within 500 micro from T ZBTB16+IL10 and fibroblast cells within 500 micro from T ZBTB16
fibro_less_500_ZBTB16_IL10 <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/fibroblast_less_500_T_ZBTB16_IL10')
fibro_less_500_ZBTB16 <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/fibroblast_less_500_T_ZBTB16')

CD_fibroblast_comparison2 <- subset(CD2, subset= cell %in% c(fibro_less_500_ZBTB16_IL10$V1, fibro_less_500_ZBTB16$V1))

a <- names(CD_fibroblast_comparison2@active.ident)
CD_fibroblast_comparison2@active.ident <- factor(unlist(lapply(as.list(CD_fibroblast_comparison2@meta.data$cell), function(x) {
  if (x %in% fibro_less_500_ZBTB16_IL10$V1) {
    'ZBTB16_IL10'
  } else {
    'ZBTB16'
  }
})))
names(CD_fibroblast_comparison2@active.ident) <- a

deg_farther_closest <- FindMarkers(CD_fibroblast_comparison2,test.use ='wilcox',ident.2 = 'ZBTB16', ident.1 = 'ZBTB16_IL10')
deg_farther_closest$diff <- "NO"
deg_farther_closest$diff[deg_farther_closest$avg_log2FC > 0 & deg_farther_closest$p_val_adj<0.05] <- "ZBTB16_IL10"
deg_farther_closest$diff[deg_farther_closest$avg_log2FC < 0 & deg_farther_closest$p_val_adj<0.05] <- "ZBTB16"
deg_farther_closest$diff <- factor(deg_farther_closest$diff)
write.table(deg_farther_closest, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_fibroblast_less_500_T_ZBTB16_IL10_less_500_T_ZBTB16.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/volcano_differential_fibroblast_less_500_T_ZBTB16_IL10_less_500_T_ZBTB10.pdf', height = 7, width =10)

ggplot2::ggplot(data=deg_farther_closest, aes(x=avg_log2FC, y=-log10(p_val), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('ZBTB16'="#28658f", 'ZBTB16_IL10'="#ef857c")) +
  geom_point(data = subset(deg_farther_closest, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.05,rownames(deg_farther_closest),'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs between fibroblasts within 500 µm from T ZBTB16+IL10+ cells and\nfibroblasts within 500 µm from T ZBTB16+ cells')

dev.off()


### DEGs between epithelial cells within 500 micro from T ZBTB16+IL10 and epithelial cells within 500 micro from T IL10
epi_less_500_ZBTB16_IL10 <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/epithelial_less_500_T_ZBTB16_IL10')
epi_less_500_IL10 <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/epithelial_less_500_T_IL10')

CD_epithelial_comparison3 <- subset(CD2, subset= cell %in% c(epi_less_500_ZBTB16_IL10$V1, epi_less_500_IL10$V1))

a <- names(CD_epithelial_comparison3@active.ident)
CD_epithelial_comparison3@active.ident <- factor(unlist(lapply(as.list(CD_epithelial_comparison3@meta.data$cell), function(x) {
  if (x %in% epi_less_500_ZBTB16_IL10$V1) {
    'ZBTB16_IL10'
  } else {
    'IL10'
  }
})))
names(CD_epithelial_comparison3@active.ident) <- a

deg_farther_closest <- FindMarkers(CD_epithelial_comparison3,test.use ='wilcox',ident.2 = 'IL10', ident.1 = 'ZBTB16_IL10')
deg_farther_closest$diff <- "NO"
deg_farther_closest$diff[deg_farther_closest$avg_log2FC > 0 & deg_farther_closest$p_val_adj<0.05] <- "ZBTB16_IL10"
deg_farther_closest$diff[deg_farther_closest$avg_log2FC < 0 & deg_farther_closest$p_val_adj<0.05] <- "IL10"
deg_farther_closest$diff <- factor(deg_farther_closest$diff)
write.table(deg_farther_closest, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_epithelial_less_500_T_ZBTB16_IL10_less_500_T_IL10.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/volcano_differential_epithelial_less_500_T_ZBTB16_IL10_less_500_T_IL10.pdf', height = 7, width =10)

ggplot2::ggplot(data=deg_farther_closest, aes(x=avg_log2FC, y=-log10(p_val), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('IL10'="#28658f", 'ZBTB16_IL10'="#ef857c")) +
  geom_point(data = subset(deg_farther_closest, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.05,rownames(deg_farther_closest),'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs between epithelial cells within 500 µm from T ZBTB16+IL10+ cells and\nepithelial cells within 500 µm from T IL10+ cells')

dev.off()

### DEGs between fibroblasts within 500 nm from T ZBTB16+IL10 and fibroblasts within 500 nm from T IL10
fibro_less_500_ZBTB16_IL10 <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/fibroblast_less_500_T_ZBTB16_IL10')
fibro_less_500_IL10 <- read.table('/Users/paolamaragno/OneDrive - Università degli Studi di Milano/Policlinico/analisi/spatial/spatial/CD/results_CD/fibroblast_less_500_T_IL10')

CD_fibroblast_comparison3 <- subset(CD2, subset= cell %in% c(fibro_less_500_ZBTB16_IL10$V1, fibro_less_500_IL10$V1))

a <- names(CD_fibroblast_comparison3@active.ident)
CD_fibroblast_comparison3@active.ident <- factor(unlist(lapply(as.list(CD_fibroblast_comparison3@meta.data$cell), function(x) {
  if (x %in% fibro_less_500_ZBTB16_IL10$V1) {
    'ZBTB16_IL10'
  } else {
    'IL10'
  }
})))
names(CD_fibroblast_comparison3@active.ident) <- a

deg_farther_closest <- FindMarkers(CD_fibroblast_comparison3,test.use ='wilcox',ident.2 = 'IL10', ident.1 = 'ZBTB16_IL10')
deg_farther_closest$diff <- "NO"
deg_farther_closest$diff[deg_farther_closest$avg_log2FC > 0 & deg_farther_closest$p_val_adj<0.05] <- "ZBTB16_IL10"
deg_farther_closest$diff[deg_farther_closest$avg_log2FC < 0 & deg_farther_closest$p_val_adj<0.05] <- "IL10"
deg_farther_closest$diff <- factor(deg_farther_closest$diff)
write.table(deg_farther_closest, file='/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/DEGs_fibroblast_less_500_T_ZBTB16_IL10_less_500_T_IL10.txt',col.names =FALSE, row.names =FALSE, quote = FALSE, sep='\n')

pdf('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/volcano_differential_fibroblast_less_500_T_ZBTB16_IL10_less_500_T_IL10.pdf', height = 7, width =10)

ggplot2::ggplot(data=deg_farther_closest, aes(x=avg_log2FC, y=-log10(p_val), fill=diff)) +
  geom_point(size=2, shape=21) + 
  scale_fill_manual(values=c('IL10'="#28658f", 'ZBTB16_IL10'="#ef857c")) +
  geom_point(data = subset(deg_farther_closest, diff != 'NO'), col = "black", shape=21, size=2) + 
  theme_classic() +
  geom_text_repel(aes(label=ifelse(p_val_adj<0.05,rownames(deg_farther_closest),'')), size =4, max.overlaps = 1000)  +
  geom_hline(yintercept=-log10(0.05), col="red", linewidth=0.3) + plot_annotation(
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = 'bold'))) +
  theme(legend.title = element_blank()) +
  ggtitle('DEGs between fibroblasts within 500 µm from T ZBTB16+IL10+ cells and\nfibroblasts within 500 µm from T IL10+ cells')

dev.off()

# for cluster annotation
# gene_cells <- read.csv('gene_cells.csv', header = FALSE)
# # 3552 
# gene_names <- read.csv('gene_names.csv', header = TRUE)[,1]
# cell_names <- read.csv('cell_names.csv', header = TRUE)[,1]
# 
# colnames(gene_cells) <- gene_names
# rownames(gene_cells) <- cell_names
# 
# #unique(colnames(gene_cells)[apply(gene_cells,1,which.max)])
# 
# # CD4 memory
# CXCR6_cells <- rownames(gene_cells[gene_cells$CXCR6 >= quantile(gene_cells$CXCR6, p=0.98),])
# DPP4_cells <- rownames(gene_cells[gene_cells$DPP4 >= quantile(gene_cells$DPP4, p=0.98),])
# IL12RB2_cells <- rownames(gene_cells[gene_cells$IL12RB2 >= quantile(gene_cells$IL12RB2, p=0.98),])
# ITGAE_cells <- rownames(gene_cells[gene_cells$ITGAE >= quantile(gene_cells$ITGAE, p=0.98),])
# FLT3LG_cells <- rownames(gene_cells[gene_cells$FLT3LG >= quantile(gene_cells$FLT3LG, p=0.98),])
# TNFRSF18_cells <- rownames(gene_cells[gene_cells$TNFRSF18 >= quantile(gene_cells$TNFRSF18, p=0.98),])
# ALOX5AP_cells <- rownames(gene_cells[gene_cells$ALOX5AP >= quantile(gene_cells$ALOX5AP, p=0.98),])
# MAF_cells <- rownames(gene_cells[gene_cells$MAF >= quantile(gene_cells$MAF, p=0.98),])
# CCL20_cells <- rownames(gene_cells[gene_cells$CCL20 >= quantile(gene_cells$CCL20, p=0.98),])
# GNLY_cells <- rownames(gene_cells[gene_cells$GNLY >= quantile(gene_cells$GNLY, p=0.98),])
# ITGA1_cells <- rownames(gene_cells[gene_cells$ITGA1 >= quantile(gene_cells$ITGA1, p=0.98),])
# IL2RG_cells <- rownames(gene_cells[gene_cells$IL2RG >= quantile(gene_cells$IL2RG, p=0.98),])
# STAT1_cells <- rownames(gene_cells[gene_cells$STAT1 >= quantile(gene_cells$STAT1, p=0.98),])
# #KRT6A.B.C_cells <- rownames(gene_cells[gene_cells$KRT6A.B.C >= quantile(gene_cells$KRT6A.B.C, p=0.98),])
# TNFSF13B_cells <- rownames(gene_cells[gene_cells$TNFSF13B >= quantile(gene_cells$TNFSF13B, p=0.98),])
# BCL2_cells <- rownames(gene_cells[gene_cells$BCL2 >= quantile(gene_cells$BCL2, p=0.98),])
# #FXPB11_cells <- rownames(gene_cells[gene_cells$FXPB11 >= quantile(gene_cells$FXPB11, p=0.98),])
# ANXA1_cells <- rownames(gene_cells[gene_cells$ANXA1 >= quantile(gene_cells$ANXA1, p=0.98),])
# IFI44L_cells <- rownames(gene_cells[gene_cells$IFI44L >= quantile(gene_cells$IFI44L, p=0.98),])
# IL17A_cells <- rownames(gene_cells[gene_cells$IL17A >= quantile(gene_cells$IL17A, p=0.98),])
# PLCG1_cells <- rownames(gene_cells[gene_cells$PLCG1 >= quantile(gene_cells$PLCG1, p=0.98),])
# 
# TCD4_memory <- Reduce(union, list(CXCR6_cells,DPP4_cells,IL12RB2_cells,ITGAE_cells,FLT3LG_cells,
#                                   TNFRSF18_cells,ALOX5AP_cells,MAF_cells,CCL20_cells, GNLY_cells,
#                                   ITGA1_cells,IL2RG_cells,STAT1_cells,TNFSF13B_cells,
#                                   BCL2_cells,ANXA1_cells,IFI44L_cells,IL17A_cells,PLCG1_cells))
# 
# # CD4 naive 
# CCR7_cells <- rownames(gene_cells[gene_cells$CCR7 >= quantile(gene_cells$CCR7, p=0.98),])
# IL6ST_cells <- rownames(gene_cells[gene_cells$IL6ST >= quantile(gene_cells$IL6ST, p=0.98),])
# VSIR_cells <- rownames(gene_cells[gene_cells$VSIR >= quantile(gene_cells$VSIR, p=0.98),])
# GPR183_cells <- rownames(gene_cells[gene_cells$GPR183 >= quantile(gene_cells$GPR183, p=0.98),])
# CD48_cells <- rownames(gene_cells[gene_cells$CD48 >= quantile(gene_cells$CD48, p=0.98),])
# ATM_cells <- rownames(gene_cells[gene_cells$ATM >= quantile(gene_cells$ATM, p=0.98),])
# CTLA4_cells <- rownames(gene_cells[gene_cells$CTLA4 >= quantile(gene_cells$CTLA4, p=0.98),])
# ICA1_cells <- rownames(gene_cells[gene_cells$ICA1 >= quantile(gene_cells$ICA1, p=0.98),])
# IL16_cells <- rownames(gene_cells[gene_cells$IL16 >= quantile(gene_cells$IL16, p=0.98),])
# IL6R_cells <- rownames(gene_cells[gene_cells$IL6R >= quantile(gene_cells$IL6R, p=0.98),])
# PDCD1_cells <- rownames(gene_cells[gene_cells$PDCD1 >= quantile(gene_cells$PDCD1, p=0.98),])
# NLRP1_cells <- rownames(gene_cells[gene_cells$NLRP1 >= quantile(gene_cells$NLRP1, p=0.98),])
# PLAC8_cells <- rownames(gene_cells[gene_cells$PLAC8 >= quantile(gene_cells$PLAC8, p=0.98),])
# 
# TCD4_naive <- Reduce(union, list(CCR7_cells,IL6ST_cells,VSIR_cells,GPR183_cells,CD48_cells,
#                                  ATM_cells,CTLA4_cells,ICA1_cells,IL16_cells,IL6R_cells,PDCD1_cells,
#                                  NLRP1_cells,PLAC8_cells))
# 
# # CD8 memory 
# GZMK_cells <- rownames(gene_cells[gene_cells$GZMK >= quantile(gene_cells$GZMK, p=0.98),])
# GZMH_cells <- rownames(gene_cells[gene_cells$GZMH >= quantile(gene_cells$GZMH, p=0.98),])
# CCL4.L1.L2_cells <- rownames(gene_cells[gene_cells$CCL4.L1.L2 >= quantile(gene_cells$CCL4.L1.L2, p=0.98),])
# CCL3.L1.L3_cells <- rownames(gene_cells[gene_cells$CCL3.L1.L3 >= quantile(gene_cells$CCL3.L1.L3, p=0.98),])
# CCR5_cells <- rownames(gene_cells[gene_cells$CCR5 >= quantile(gene_cells$CCR5, p=0.98),])
# PRF1_cells <- rownames(gene_cells[gene_cells$PRF1 >= quantile(gene_cells$PRF1, p=0.98),])
# CTSW_cells <- rownames(gene_cells[gene_cells$CTSW >= quantile(gene_cells$CTSW, p=0.98),])
# DUSP2_cells <- rownames(gene_cells[gene_cells$DUSP2 >= quantile(gene_cells$DUSP2, p=0.98),])
# RUNX3_cells <- rownames(gene_cells[gene_cells$RUNX3 >= quantile(gene_cells$RUNX3, p=0.98),])
# HCST_cells <- rownames(gene_cells[gene_cells$HCST >= quantile(gene_cells$HCST, p=0.98),])
# CXCR3_cells <- rownames(gene_cells[gene_cells$CXCR3 >= quantile(gene_cells$CXCR3, p=0.98),])
# HCST_cells <- rownames(gene_cells[gene_cells$HCST >= quantile(gene_cells$HCST, p=0.98),])
# OASL_cells <- rownames(gene_cells[gene_cells$OASL >= quantile(gene_cells$OASL, p=0.98),])
# COTL1_cells <- rownames(gene_cells[gene_cells$COTL1 >= quantile(gene_cells$COTL1, p=0.98),])
# IL10RA_cells <- rownames(gene_cells[gene_cells$IL10RA >= quantile(gene_cells$IL10RA, p=0.98),])
# CYTOR_cells <- rownames(gene_cells[gene_cells$CYTOR >= quantile(gene_cells$CYTOR, p=0.98),])
# SH3BGRL3_cells <- rownames(gene_cells[gene_cells$SH3BGRL3 >= quantile(gene_cells$SH3BGRL3, p=0.98),])
# TNFRSF1B_cells <- rownames(gene_cells[gene_cells$TNFRSF1B >= quantile(gene_cells$TNFRSF1B, p=0.98),])
# TTN_cells <- rownames(gene_cells[gene_cells$TTN >= quantile(gene_cells$TTN, p=0.98),])
# AKT1_cells <- rownames(gene_cells[gene_cells$AKT1 >= quantile(gene_cells$AKT1, p=0.98),])
# 
# TCD8_memory <- Reduce(union, list(GZMK_cells,GZMH_cells,CCL4.L1.L2_cells,CCL3.L1.L3_cells, 
#                                   CCR5_cells, PRF1_cells,CTSW_cells, DUSP2_cells, RUNX3_cells, HCST_cells, CXCR3_cells,
#                                   OASL_cells, COTL1_cells, IL10RA_cells,CYTOR_cells, SH3BGRL3_cells, TNFRSF1B_cells,
#                                   TTN_cells,AKT1_cells))
# 
# # CD8 naive
# LINC02446_cells <- rownames(gene_cells[gene_cells$LINC02446 >= quantile(gene_cells$LINC02446, p=0.98),])
# IL2RB_cells <- rownames(gene_cells[gene_cells$IL2RB >= quantile(gene_cells$IL2RB, p=0.98),])
# WNT10B_cells <- rownames(gene_cells[gene_cells$WNT10B >= quantile(gene_cells$WNT10B, p=0.98),])
# CD53_cells <- rownames(gene_cells[gene_cells$CD53 >= quantile(gene_cells$CD53, p=0.98),])
# CCL17_cells <- rownames(gene_cells[gene_cells$CCL17 >= quantile(gene_cells$CCL17, p=0.98),])
# TGFBR2_cells <- rownames(gene_cells[gene_cells$TGFBR2 >= quantile(gene_cells$TGFBR2, p=0.98),])
# CD55_cells <- rownames(gene_cells[gene_cells$CD55 >= quantile(gene_cells$CD55, p=0.98),])
# KLF2_cells <- rownames(gene_cells[gene_cells$KLF2 >= quantile(gene_cells$KLF2, p=0.98),])
# IRF3_cells <- rownames(gene_cells[gene_cells$IRF3 >= quantile(gene_cells$IRF3, p=0.98),])
# CCL19_cells <- rownames(gene_cells[gene_cells$CCL19 >= quantile(gene_cells$CCL19, p=0.98),])
# OAS2_cells <- rownames(gene_cells[gene_cells$OAS2 >= quantile(gene_cells$OAS2, p=0.98),])
# CD44_cells <- rownames(gene_cells[gene_cells$CD44 >= quantile(gene_cells$CD44, p=0.98),])
# 
# TCD8_naive <- Reduce(union, list(LINC02446_cells,IL2RB_cells,WNT10B_cells,CD53_cells,CCL17_cells,TGFBR2_cells,CD55_cells,
#                                  KLF2_cells, IRF3_cells,CCL19_cells,OAS2_cells,CD44_cells))
# 
# # T reg
# FOXP3_cells <- rownames(gene_cells[gene_cells$FOXP3 >= quantile(gene_cells$FOXP3, p=0.98),])
# IL2RA_cells <- rownames(gene_cells[gene_cells$IL2RA >= quantile(gene_cells$IL2RA, p=0.98),])
# SRGN_cells <- rownames(gene_cells[gene_cells$SRGN >= quantile(gene_cells$SRGN, p=0.98),])
# RAC2_cells <- rownames(gene_cells[gene_cells$RAC2 >= quantile(gene_cells$RAC2, p=0.98),])
# RGS1_cells <- rownames(gene_cells[gene_cells$RGS1 >= quantile(gene_cells$RGS1, p=0.98),])
# LY75_cells <- rownames(gene_cells[gene_cells$LY75 >= quantile(gene_cells$LY75, p=0.98),])
# PFN1_cells <- rownames(gene_cells[gene_cells$PFN1 >= quantile(gene_cells$PFN1, p=0.98),])
# CD27_cells <- rownames(gene_cells[gene_cells$CD27 >= quantile(gene_cells$CD27, p=0.98),])
# NR3C1_cells <- rownames(gene_cells[gene_cells$NR3C1 >= quantile(gene_cells$NR3C1, p=0.98),])
# IRF4_cells <- rownames(gene_cells[gene_cells$IRF4 >= quantile(gene_cells$IRF4, p=0.98),])
# 
# Treg <- Reduce(union, list(FOXP3_cells,IL2RA_cells,SRGN_cells,RAC2_cells,RGS1_cells,LY75_cells,
#                            PFN1_cells,CD27_cells,NR3C1_cells,IRF4_cells))
# 
# 
# venn.diagram(
#   x = list(TCD4_memory, TCD4_naive, TCD8_memory, TCD8_naive, Treg),
#   category.names = c('TCD4_memory', 'TCD4_naive', 'TCD8_memory', 'TCD8_naive','Treg'),
#   filename='./venn_diagramm_0.98.png',
#   output=TRUE,
#   height = 5000, 
#   width = 5000
# )
# 
# length(Reduce(union, list(TCD4_memory, TCD4_naive, TCD8_memory, TCD8_naive)))
# #2513 di 3552
# 
# # reg
# FOXP3_cells <- rownames(gene_cells[gene_cells$FOXP3 >= quantile(gene_cells$FOXP3, p=0.98),])
# CTLA4_cells <- rownames(gene_cells[gene_cells$CTLA4 >= quantile(gene_cells$CTLA4, p=0.98),])
# ITGB8_cells <- rownames(gene_cells[gene_cells$ITGB8 >= quantile(gene_cells$ITGB8, p=0.98),])
# 
# 
# 
# 
# ##### EPITHELIAL
# gene_cells <- read.csv('gene_cells_epithelial.csv', header = FALSE)
# # 3725
# gene_names <- read.csv('gene_names_epithelial.csv', header = TRUE)[,1]
# cell_names <- read.csv('cell_names_epithelial.csv', header = TRUE)[,1]
# 
# colnames(gene_cells) <- gene_names
# rownames(gene_cells) <- cell_names
# 
# #unique(colnames(gene_cells)[apply(gene_cells,1,which.max)])
# 
# # goblet
# FCGBP_cells <- rownames(gene_cells[gene_cells$FCGBP >= quantile(gene_cells$FCGBP, p=0.8),])
# AGR2_cells <- rownames(gene_cells[gene_cells$AGR2 >= quantile(gene_cells$AGR2, p=0.8),])
# 
# goblet <- Reduce(intersect, list(FCGBP_cells,AGR2_cells))
# 
# # paneth
# REG1A_cells <- rownames(gene_cells[gene_cells$REG1A >= quantile(gene_cells$REG1A, p=0.9),])
# IFI6_cells <- rownames(gene_cells[gene_cells$IFI6 >= quantile(gene_cells$IFI6, p=0.9),])
# LYZ_cells <- rownames(gene_cells[gene_cells$LYZ >= quantile(gene_cells$LYZ, p=0.9),])
# SERPINA1_cells <- rownames(gene_cells[gene_cells$SERPINA1 >= quantile(gene_cells$SERPINA1, p=0.9),])
# LCN2_cells <- rownames(gene_cells[gene_cells$LCN2 >= quantile(gene_cells$LCN2, p=0.9),])
# 
# paneth <- Reduce(intersect, list(REG1A_cells,IFI6_cells,LYZ_cells,SERPINA1_cells,LCN2_cells))
# 
# gene_cells$cell_type <- 'epithelial'
# gene_cells[paneth,'cell_type'] <- 'paneth'
# gene_cells[goblet,'cell_type'] <- 'goblet'
# 
# write.table(cell_type, file='./cell_type_epithelial.txt', quote=FALSE, row.names = FALSE, col.names = FALSE)
# 

