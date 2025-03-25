#!/usr/bin/env python
# coding: utf-8

# In[1]:


from pathlib import Path
from copy import deepcopy

import numpy as np

import matplotlib as mp
import matplotlib.pyplot as plt
import seaborn as sns 
import scanpy as sc
import squidpy as sq
import pandas as pd
import anndata as ad
import scipy as sp
from scipy.spatial import distance_matrix
#import commot as ct
from scipy.stats import ttest_ind
from scipy.sparse import csr_array, tril

sc.logging.print_header()


# In[2]:


import gseapy
import matplotlib.pyplot as plt
import warnings
import os
import urllib.request
from matplotlib_venn import venn3_unweighted

#from nichecompass.models import NicheCompass


# In[3]:


import novae


# In[2]:


seed = 10
np.random.seed(seed)


# ### Pre-processing

# In[3]:


nanostring_dir = Path('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD').resolve()
sample_dir = nanostring_dir


# In[4]:


adata = sq.read.nanostring(
    path=sample_dir,
    counts_file="R5941_ColonTMA_exprMat_file_CD.csv",
    meta_file="R5941_ColonTMA_metadata_file_CD.csv",
    fov_file="R5941_ColonTMA_fov_positions_file_CD.csv",
)


# In[5]:


# Obtain the control probes using their names prefixed with “NegPrb-“.

adata.var["NegPrb"] = adata.var_names.str.startswith("NegPrb")
sc.pp.calculate_qc_metrics(adata, qc_vars=["NegPrb"], inplace=True)


# In[6]:


adata.obs["total_counts_NegPrb"].sum() / adata.obs["total_counts"].sum() * 100

#The percentage of unassigned “NegPrb” transcripts can be calculated from the calculated qc metrics. This can later be used to estimate 
# false discovery rate.


# In[7]:


t = ['NegPrb01',
'NegPrb02',
'NegPrb03',
'NegPrb04',
'NegPrb05',
'NegPrb06',
'NegPrb07',
'NegPrb08',
'NegPrb09',
'NegPrb10']


# In[8]:


adata = adata[:,~adata.var.index.isin(t)]


# In[9]:


exp = pd.DataFrame(adata.X.toarray(),index = adata.obs.index.values, columns = adata.var.index.values)
gene_cell = pd.DataFrame(index = exp.columns, columns = ['count'])

for g in exp.columns:
    gene_cell.loc[g] = len(exp.loc[:,g][exp.loc[:,g]!=0])
    


# In[12]:


np.min(gene_cell.loc[:,'count'])


# In[10]:


fig, axs = plt.subplots(1, 4, figsize=(17, 4))

axs[0].set_title("Total transcripts per cell")
sns.histplot(
    adata.obs["total_counts"],
    kde=False,
    ax=axs[0],
)
axs[0].set_xlim(0, 500)

axs[1].set_title("Unique transcripts per cell")
sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    ax=axs[1],
)

axs[2].set_title("Transcripts per FOV")
sns.histplot(
    adata.obs.groupby("fov").sum()["total_counts"],
    kde=False,
    bins =20,
    ax=axs[2],
)

axs[3].set_title("Number of cells expressing same transcript")
sns.histplot(
    gene_cell['count'],
    kde=False,
    bins =60,
    ax=axs[3]
)
axs[3].set(xlabel='Number of cells expressing x transcripts')
axs[3].set_xlim(0, 5000)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/quality.pdf',bbox_inches="tight")


# In[10]:


# Filter the cells based on the minimum number of counts required using scanpy.pp.filter_cells. Filter the genes based on the minimum 
# number of cells required with scanpy.pp.filter_genes. The parameters for both were specified based on the plots above. This filtering 
# is quite conservative, more relaxed settings might also be applicable. Other criteria for filtering cells could be area, 
# immunofluorescence signal or number of unique transcripts.

sc.pp.filter_cells(adata, min_counts=25)
sc.pp.filter_genes(adata, min_cells=700)


# In[11]:


# Normalize counts per cell 
sc.pp.normalize_total(adata, inplace = True)

sc.pp.log1p(adata)

# PCA
sc.pp.pca(adata, random_state=seed)

# compute a beigthorhood graph
sc.pp.neighbors(adata,random_state=seed)

# embed the neighborhood graph of the data
sc.tl.umap(adata,random_state=seed)

# cluster the cells into subgroups 
sc.tl.leiden(adata, resolution=1.2, random_state=seed)


# In[17]:


num_cells_per_patient = {}

for p in pd.unique(adata.obs.Patients):
    num_cells_per_patient[p] = adata[adata.obs.Patients == p].shape[0]


# In[85]:


num_cells_per_patient


# In[86]:


num_cells_per_fov = {}

for f in pd.unique(adata.obs.fov):
    num_cells_per_fov[f] = adata[adata.obs.fov == f].shape[0]


# In[87]:


num_cells_per_fov


# In[12]:


palette = ['#0F2080', '#006600', '#AA3377', '#BBBBBB', '#AEC7E8', '#CCBB44', '#F7B6D2','#000000',
          '#0077BB', '#EE7733', '#33BBEE', '#CC3311', '#7F7F7F', '#FFB000','#601A4A','#785EF0','#66CC00', '#CCFF99', 
          '#660000', '#E5CCFF','#FFFFCC','#CCFFE5']

sns.color_palette(palette)


# In[13]:


sc.pl.umap(
    adata, show=False,
    color= "leiden",
    palette=sns.color_palette(palette)
)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/umap_clusters_not_named.pdf',bbox_inches="tight")


# In[16]:


cat_name = "leiden"
sig_leiden = pd.DataFrame(
    columns=adata.var_names, index=adata.obs[cat_name].cat.categories
)

# for each cluster, extract the average expression of each gene in the cells of that cluster
for clust in adata.obs[cat_name].cat.categories:
    sig_leiden.loc[clust] = adata[adata.obs[cat_name].isin([clust]), :].X.mean(0)
# transpose this matrix so that the columns are the clusters and the rows are the cosmx genes
sig_leiden = sig_leiden.transpose()
sig_leiden.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/gene_counts_all_clusters_non_named_after_subclustering.xlsx')


# In[ ]:


# EPCAM (11,10,2)
# COL1A1 (4), DCN (4)
# PECAM1 (7), CD93 (7), TIE1 (7)
# IGKC (3,9), IGHG1 (3), IGHA1 -> plasma
# MS4A1 (0,6), CD19 (0,6) -> B
# CD44 (12,1,0) -> B memory/reg
# CD68 (8), CD163 (8), CLEC10A (8) -> macrophage
# CD3D (1), CD3E (1), CD3G (1,14), CD2 (1), FYN (1) -> T
# MYH11 (5), ACTG2 (5), ACTA2 (5), TAGLN (5) -> smooth
# VIM (13), S100B (13) -> enteric glia
# CPA3 (12), KIT (12), TPSAB1.B2 (12) -> mast


# In[36]:


sig_leiden.loc['CD19',].sort_values(ascending=False)


# In[28]:


num_top_genes = 50

l = list()

for inst_cluster in sig_leiden.columns.tolist():
    # for each leiden cluster (column names of sig_leiden dataframe), sort the genes in decreasing order of expression and keep the 30 genes with 
# highest expression.
    top_genes = (
        sig_leiden[inst_cluster]
        .sort_values(ascending=False)
        .index.tolist()[:num_top_genes]
    )
    l.append(inst_cluster)
    l.append((sig_leiden[inst_cluster][top_genes].to_string()))

with open('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/top_50_genes_clusters_non_named_before_subclustering.txt', 'w') as fp:
    for item in l:
        # write each item on a new line
        fp.write("%s\n" % item)
        fp.write("\n\n")


# In[14]:


map_dict = {'0':'B reg', '1':'Epithelial','2':'T', '3':'Smooth muscle cells', '4':'Plasma IgG', '5':'Plasma IgG', '6':'Fibroblast', '7':'Plasma IgA', 
            '8':'B','9':'Macrophage','10':'Endothelial', '11':'Plasma IgG','12':'Plasma IgA','13':'Plasma IgG', '14':'Smooth muscle cells', 
           '15':'Mast cells', '16':'Epithelial', '17': 'Enteric glia cells', '18':'Epithelial', '19':'Epithelial','20':'Plasma IgA','21':'Epithelial' }

adata.obs['new_leiden'] = adata.obs['leiden'].map(map_dict).astype('category')
adata.obs.new_leiden


# In[15]:


palette = ['#0F2080', '#228833', '#AA3377', '#BBBBBB', '#AEC7E8', '#CCBB44', '#F7B6D2',
          '#0077BB', '#EE7733', '#33BBEE', '#CC3311', '#7F7F7F']

sns.color_palette(palette)


# In[16]:


sc.pl.umap(
    adata, show=False,
    color= "new_leiden", palette=sns.color_palette(palette)
)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/umap_clusters_named.pdf',bbox_inches="tight")


# In[15]:


cat_name = "new_leiden"
sig_leiden = pd.DataFrame(
    columns=adata.var_names, index=adata.obs[cat_name].cat.categories
)

# for each cluster, extract the average expression of each gene in the cells of that cluster
for clust in adata.obs[cat_name].cat.categories:
    sig_leiden.loc[clust] = adata[adata.obs[cat_name].isin([clust]), :].X.mean(0)
# transpose this matrix so that the columns are the clusters and the rows are the cosmx genes
sig_leiden = sig_leiden.transpose()

sig_leiden.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/gene_counts_all_clusters_named.xlsx')


# In[62]:


num_top_genes = 50

l = list()

for inst_cluster in sig_leiden.columns.tolist():
    # for each leiden cluster (column names of sig_leiden dataframe), sort the genes in decreasing order of expression and keep the 30 genes with 
# highest expression.
    top_genes = (
        sig_leiden[inst_cluster]
        .sort_values(ascending=False)
        .index.tolist()[:num_top_genes]
    )
    l.append(inst_cluster)
    l.append((sig_leiden[inst_cluster][top_genes].to_string()))

with open('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/top_50_genes_clusters_named.txt', 'w') as fp:
    for item in l:
        # write each item on a new line
        fp.write("%s\n" % item)
        fp.write("\n\n")


# In[17]:


level = []

for p in adata.obs.Patients:
    if p in ['Patient_11','Patient_13', 'Patient_15', 'Patient_19', 'Patient_20']:
        level.append('LOW')
    else:
        level.append('HIGH')

adata.obs['level'] = level
adata.obs['level'] = adata.obs['level'].astype("category")


# In[25]:


sc.pl.umap(
    adata, show=False,
    color= "level", palette=sns.color_palette(['#0F2080', '#AA3377'], 2)
)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/umap_level.pdf',bbox_inches="tight")


# In[26]:


sc.pl.umap(
    adata, show=False,
    color= "Patients", palette=sns.color_palette(palette, 10)
)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/umap_patient.pdf',bbox_inches="tight")


# In[27]:


adata.obs.to_csv('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/metadata_adata.csv')


# In[28]:


adata.var.to_csv("/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/gene_names_adata.csv")


# In[ ]:


pd.unique(adata.obs[adata.obs['fov'] == '41'].level)


# In[ ]:


adata.uns['spatial'].keys()


# In[ ]:


fov_names = adata.obs['fov']
fov_names = adata.obs['fov']
fov_names = fov_names.astype("string")
patients_names = adata.obs['Patients']
patients_names = patients_names.astype("string")
level_names = adata.obs['level']
level_names = level_names.astype("string")


# In[ ]:


list(set('fov ' + fov_names + ' - ' + patients_names + ' - ' + level_names))


# In[18]:


t = ['fov 21 - Patient_11 - LOW',
 'fov 23 - Patient_12 - HIGH',
 'fov 24 - Patient_12 - HIGH',
 'fov 25 - Patient_13 - LOW',
 'fov 27  - Patient_14 - HIGH',
 'fov 28 - Patient_14 - HIGH',
 'fov 29 - Patient_15 - LOW',
 'fov 30 - Patient_15 - LOW',
 'fov 31 - Patient_16 - HIGH',
 'fov 32 - Patient_16 - HIGH',
 'fov 33 - Patient_17 - HIGH',
 'fov 34 - Patient_17 - HIGH',
 'fov 35 - Patient_18 - HIGH',
 'fov 36 - Patient_18 - HIGH',
 'fov 37 - Patient_19 - LOW', 
 'fov 38 - Patient_19 - LOW',
 'fov 39 - Patient_20 - LOW',
 'fov 40 - Patient_20 - LOW',
 'fov 41 - Patient_20 - LOW']


# In[38]:


sq.pl.spatial_segment(
    adata,
    color='new_leiden',
    library_key="fov",
    seg_cell_id="cell_ID",
    seg_outline=True,
    img=False,
    scalebar_dx=1.0,
    title= t,
    scalebar_kwargs={"scale_loc": "bottom", "location": "lower right"},
)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/fov_patients_level_senza_back.pdf',bbox_inches="tight")


# ### DEGs with scanpy

# In[18]:


sc.tl.rank_genes_groups(adata, 'new_leiden', method='t-test', key_added = "t-test")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key = "t-test",
                       save = 'rank_genes_groups_new_leiden_t_test.pdf')


# In[19]:


sc.tl.rank_genes_groups(adata, 'new_leiden', method='t-test_overestim_var', key_added = "t-test_ov")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key = "t-test_ov", save='rank_genes_groups_new_leiden_t_test_over_var.pdf')


# In[20]:


sc.tl.rank_genes_groups(adata, 'new_leiden', method='wilcoxon', key_added = "wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key="wilcoxon", save='rank_genes_groups_new_leiden_wilcoxon.pdf')


# In[21]:


# Take all significant DE genes for B cells with each test and compare the overlap.

wc = sc.get.rank_genes_groups_df(adata, group='B', key='wilcoxon', pval_cutoff=0.01, log2fc_min=0)['names']
tt = sc.get.rank_genes_groups_df(adata, group='B', key='t-test', pval_cutoff=0.01, log2fc_min=0)['names']
tt_ov = sc.get.rank_genes_groups_df(adata, group='B', key='t-test_ov', pval_cutoff=0.01, log2fc_min=0)['names']

venn3_unweighted([set(wc),set(tt),set(tt_ov)], ('Wilcox','T-test','T-test_ov'),)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/venn_B.pdf',bbox_inches="tight")


# In[22]:


clusters = ["T","B"   ,"B reg" ,"Plasma IgA", "Plasma IgG","Macrophage" ,"Mast cells","Smooth muscle cells" ,"Endothelial" ,"Epithelial","Fibroblast"  , 
                                     "Enteric glia cells"]
adata.obs["clusters_ordered"] = pd.Categorical(
    values=adata.obs.new_leiden, categories=clusters, ordered=True
)


# In[23]:


markers = ['CD3D',	'CD3E',	'CD3G', 'CD2', 'FYN','CD19', 'MS4A1',	'CD44',
             'IGKC','IGHA1','IGHG1', 
             'CD68', 'CD163', 'CLEC10A', 'CPA3', 'KIT',
             'TPSAB1.B2', 'MYH11', 'ACTG2', 'ACTA2', 'TAGLN',
             'PECAM1', 'CD93', 'TIE1','EPCAM','COL1A1', 'DCN', 'VIM', 'S100B']

sc.pl.dotplot(adata, markers, groupby='clusters_ordered', standard_scale="var", save=True)


# In[62]:


vp= sc.pl.stacked_violin(adata, markers, groupby='clusters_ordered', return_fig=True)
vp.add_totals().style(ylim=(0,5)).show()
vp.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/violin.pdf',bbox_inches="tight")


# In[141]:


df = sc.get.rank_genes_groups_df(adata, 'B', key = "wilcoxon")
np.array(df['names'][0:5])


# In[67]:


markers = ['PFN1', 'MARCKSL1', 'CLU', 'HMGN2', 'BASP1','MS4A1','CD19', 'CD74', 'CD37', 'CXCR4', 'HLA.DRB','CD44',
           'PECAM1', 'DUSP1', 'VIM', 'IGFBP7', 'SPARCL1','CD93','TIE1', 'CDH19', 'CRYAB', 'S100A6', 'S100B',
           'REG1A', 'OLFM4', 'PIGR', 'KRT8', 'KRT18', 'EPCAM','COL1A1', 'COL6A2', 'COL3A1', 'COL6A1', 'COL1A2',
           'DCN','SAT1', 'HLA.DRA', 'LYZ', 'CD68','CD163','CLEC10A','TPSAB1.B2', 'IL1RL1', 'KIT', 'RGS2', 'CPA3',
           'IGHA1', 'IGKC', 'JCHAIN', 'MZB1', 'XBP1','IGHG1', 'IGHG2', 'LPAR5',
           'TAGLN', 'MYL9', 'ACTA2', 'TPM2', 'MYH11','ACTG2','B2M', 'BTG1', 'IL7R', 'FYB1', 'TPT1','CD3D', 'CD3E', 'CD3G', 'CD2', 'FYN']


# In[28]:


markers = ['CD19', 'MS4A1',	'CD44','PECAM1', 'CD93', 'TIE1','VIM', 'S100B',
           'EPCAM','COL1A1', 'DCN','CD68', 'CD163', 'CLEC10A', 'CPA3', 'KIT',
             'TPSAB1.B2','IGKC','IGHA1','IGHG1','MYH11', 'ACTG2', 'ACTA2', 'TAGLN','CD3D',	'CD3E',	'CD3G', 'CD2', 'FYN']


# In[38]:


sc.pl.heatmap(adata, var_names=markers, groupby='new_leiden',show_gene_labels=True, save = 'expression_heatmap.pdf',figsize=(10,7))


# In[84]:


mp = sc.pl.matrixplot(adata, markers, 'new_leiden', return_fig=True,figsize=(12,5))
mp.add_totals().style(edge_color='black')
mp.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/matrixplot.pdf',bbox_inches="tight")


# In[35]:


sc.get.rank_genes_groups_df(adata, group="T", key='wilcoxon').head(5)


# In[37]:


dc_cluster_genes = sc.get.rank_genes_groups_df(adata, group="T", key='wilcoxon').head(5)["names"]
sc.pl.umap(
    adata,
    color=[*dc_cluster_genes, "new_leiden"],
    legend_loc="on data",
    frameon=False,
    ncols=3, save=True
)


# In[89]:


sc.settings.set_figure_params(dpi=80)


# In[91]:


sc.tl.rank_genes_groups(adata, 'new_leiden', groups=['Fibroblast'], reference='Smooth muscle cells', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groups=['Fibroblast'], n_genes=20, save='fibro_vs_smooth.pdf')

sc.pl.rank_genes_groups_violin(adata, groups='Fibroblast', n_genes=10,save='fibro_vs_smooth_violin.pdf')
mynames = [x[0] for x in adata.uns['rank_genes_groups']['names'][:10]]
sc.pl.stacked_violin(adata, mynames, groupby = 'new_leiden',save='fibro_vs_smooth_violin_all_clusters.pdf')


# In[120]:


cl2 = adata[adata.obs['new_leiden']=='B',:]
genes = sc.get.rank_genes_groups_df(cl2, group='B', key='wilcoxon')['names'][:5]
df = sc.get.obs_df(adata, genes.tolist() + ['new_leiden'], use_raw=False)
df2 = df.melt(id_vars=["new_leiden"], value_vars=genes)

g = sns.catplot(x = "new_leiden", y = "value", hue='new_leiden',  kind = 'violin', col = "variable", data = df2, inner=None)
g.tick_params(axis='x', rotation=90)
g.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/violin_B_markers.pdf')


# In[124]:


cl2 = adata[adata.obs['new_leiden']=='B',:]
cl2.obs['level'].value_counts()

sc.tl.rank_genes_groups(cl2, 'level', method='wilcoxon', key_added = "wilcoxon")
sc.pl.rank_genes_groups(cl2, n_genes=25, sharey=False, key="wilcoxon")


# In[126]:


genes1 = sc.get.rank_genes_groups_df(cl2, group='LOW', key='wilcoxon')['names'][:5]
genes2 = sc.get.rank_genes_groups_df(cl2, group='HIGH', key='wilcoxon')['names'][:5]
genes = genes1.tolist() +  genes2.tolist() 
df = sc.get.obs_df(adata, genes + ['new_leiden','level'], use_raw=False)
df2 = df.melt(id_vars=["new_leiden",'level'], value_vars=genes)

g = sns.catplot(x = "new_leiden", y = "value", hue='level',  kind = 'violin', col = "variable", data = df2, col_wrap=4,inner=None)
g.tick_params(axis='x', rotation=90)
g.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/violin_B_markers_high_low.pdf')


# In[131]:


# lets check how the top DGEs are expressed across levels
genes1 = sc.get.rank_genes_groups_df(cl2, group='LOW', key='wilcoxon')['names'][:5]
genes2 = sc.get.rank_genes_groups_df(cl2, group='HIGH', key='wilcoxon')['names'][:5]
genes = genes1.tolist() +  genes2.tolist() 

sc.pl.violin(cl2, genes1, groupby='level', save='violin_B_low.pdf')
sc.pl.violin(cl2, genes2, groupby='level', save='violin_B_high.pdf')


# In[133]:


genes1 = sc.get.rank_genes_groups_df(cl2, group='LOW', key='wilcoxon')['names'][:20]
genes2 = sc.get.rank_genes_groups_df(cl2, group='HIGH', key='wilcoxon')['names'][:20]
genes = genes1.tolist() +  genes2.tolist() 

sc.pl.dotplot(cl2,genes, groupby='level', save='dotplot_unique_high_low.pdf')


# In[142]:


df = sc.get.rank_genes_groups_df(adata, 'B', key = "wilcoxon")
np.array(df['names'][0:5])


# In[146]:


# lets check how the top DGEs in B cluster are expressed across patients
sc.pl.violin(cl2, ['PFN1', 'MARCKSL1', 'CLU', 'HMGN2', 'BASP1'], groupby='Patients', rotation=90, save = 'violin_B_patients.pdf')


# In[148]:


df = sc.get.rank_genes_groups_df(adata, 'B', key = "wilcoxon")
np.array(df['names'][0:20])


# In[149]:


sc.pl.dotplot(cl1,['PFN1', 'MARCKSL1', 'CLU', 'HMGN2', 'BASP1', 'CD74', 'MS4A1',
       'TCL1A', 'CD22', 'H4C3', 'EZR', 'P2RX5', 'ARHGDIB', 'PPIA', 'CD37',
       'PTPRCAP', 'H2AZ1', 'CD53', 'CD79A', 'SRGN'], groupby='Patients',save = 'dotplot_B_patients.pdf')


# ### Gene Set Analysis (GSA): Hypergeometric enrichment test

# In[49]:


gene_set_names = gseapy.get_library_name(organism='Human')


# In[50]:


adata_epi = adata[adata.obs.new_leiden =='Epithelial']
sc.tl.rank_genes_groups(adata_epi, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
df_adata_epi = sc.get.rank_genes_groups_df(adata_epi, group='LOW',key = "rank_genes_groups", pval_cutoff=0.1)


# In[51]:


df_adata_epi_high = df_adata_epi[df_adata_epi.logfoldchanges < 0]
df_adata_epi_low = df_adata_epi[df_adata_epi.logfoldchanges > 0]


# In[36]:


adata_epi_exp = pd.DataFrame(adata_epi.X.toarray(), index = adata_epi.obs.index, columns = adata_epi.var.index)
adata_epi_exp['level'] = adata_epi.obs.level


# In[37]:


adata_epi_exp.groupby('level', as_index=False)['CD24'].mean()


# In[52]:


# get the significant DEGs for B cluster
glist_high = df_adata_epi_high['names'].squeeze().str.strip().tolist()
glist_low = df_adata_epi_low['names'].squeeze().str.strip().tolist()
print(len(glist_high))


# In[38]:


df_adata_epi


# In[119]:


with open('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/GSA_epi_DEGs.txt', "w") as o:
    for line in glist:
        print(line, file=o)


# In[32]:


df_adata_epi.to_csv('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/GSA_epi_DEGs.txt')


# In[54]:


enr_res_high = gseapy.enrichr(gene_list=glist_high, organism='Human', gene_sets='GO_Biological_Process_2023', cutoff = 0.5)
enr_res_low = gseapy.enrichr(gene_list=glist_low, organism='Human', gene_sets='GO_Biological_Process_2023', cutoff = 0.5)


# In[55]:


df_high = pd.DataFrame(enr_res_high.results)
df_low = pd.DataFrame(enr_res_low.results)
#df.to_csv('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/GSA_epi_highDEGs_table.csv')
df_high = df_high[df_high.Term.isin(['Regulation Of Apoptotic Process (GO:0042981)','Regulation Of Fibroblast Proliferation (GO:0048145)','Epithelial To Mesenchymal Transition (GO:0001837)','Immune Response-Regulating Cell Surface Receptor Signaling Pathway (GO:0002768)','Positive Regulation Of Apoptotic Process (GO:0043065)','Negative Regulation Of Apoptotic Process (GO:0043066)','Negative Regulation Of Fibroblast Proliferation (GO:0048147)','Positive Regulation Of Fibroblast Proliferation (GO:0048146)', 'Regulation Of Oxidative Stress-Induced Intrinsic Apoptotic Signaling Pathway (GO:1902175)','Wound Healing (GO:0042060)'])]
df_low = df_low[df_low.Term.isin(['Regulation Of Apoptotic Process (GO:0042981)','Regulation Of Fibroblast Proliferation (GO:0048145)','Epithelial To Mesenchymal Transition (GO:0001837)','Immune Response-Regulating Cell Surface Receptor Signaling Pathway (GO:0002768)','Positive Regulation Of Apoptotic Process (GO:0043065)','Negative Regulation Of Apoptotic Process (GO:0043066)','Negative Regulation Of Fibroblast Proliferation (GO:0048147)','Positive Regulation Of Fibroblast Proliferation (GO:0048146)', 'Regulation Of Oxidative Stress-Induced Intrinsic Apoptotic Signaling Pathway (GO:1902175)','Wound Healing (GO:0042060)'])]


# In[58]:


df_high['-log10(Adjusted P-value)'] = -np.log10(df_high['Adjusted P-value'])
df_low['-log10(Adjusted P-value)'] = -np.log10(df_low['Adjusted P-value'])

df_high = df_high[df_high['Adjusted P-value'] < 0.1]
df_low = df_low[df_low['Adjusted P-value'] < 0.1]


# In[57]:


df_high.to_csv('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/GSA_epi_selected_table_high.csv')
df_low.to_csv('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/GSA_epi_selected_table_low.csv')


# In[59]:


plt.subplot(122)
plt.barh(width='-log10(Adjusted P-value)', y= 'Term', data = df_high)
plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/GSA_epi_selected_high.pdf',bbox_inches="tight")


# In[60]:


plt.subplot(122)
plt.barh(width='-log10(Adjusted P-value)', y= 'Term', data = df_low)
plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/GSA_epi_selected_low.pdf',bbox_inches="tight")


# In[48]:


ax= gseapy.barplot(enr_res_high.res2d,title='GO_Biological_Process_2023')
#plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/GSA_epi.pdf',bbox_inches="tight")


# In[49]:


ax= gseapy.barplot(enr_res_low.res2d,title='GO_Biological_Process_2023')


# In[58]:


ax = gseapy.dotplot(enr_res.res2d, title='GO_Biological_Process_2023',cmap='viridis_r', size=10, figsize=(3,5))
plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/GSA_dotplot_epi.pdf',bbox_inches="tight")


# In[51]:


adata_T = adata[adata.obs.new_leiden =='T']
sc.tl.rank_genes_groups(adata_T, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
df_adata_T = sc.get.rank_genes_groups_df(adata_T, group='LOW',key = "rank_genes_groups", pval_cutoff=0.05)


# In[52]:


glist = df_adata_T['names'].squeeze().str.strip().tolist()


# In[53]:


enr_res = gseapy.enrichr(gene_list=glist, organism='Human', gene_sets='GO_Biological_Process_2023', cutoff = 0.5)
enr_res.results.head()


# In[54]:


ax= gseapy.barplot(enr_res.res2d,title='GO_Biological_Process_2023')
plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/GSA_T.pdf',bbox_inches="tight")


# ### GSEA

# In[87]:


adata_epi = adata[adata.obs.new_leiden =='Epithelial']
sc.tl.rank_genes_groups(adata_epi, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
df_adata_epi = sc.get.rank_genes_groups_df(adata_epi, group='LOW',key = "rank_genes_groups", pval_cutoff=0.05)[['names','logfoldchanges']]
df_adata_epi.sort_values(by=['logfoldchanges'], inplace=True, ascending=False)

# calculate_qc_metrics will calculate number of cells per gene
sc.pp.calculate_qc_metrics(adata_epi, percent_top=None, log1p=False, inplace=True)

# filter for genes expressed in at least 30 cells.
df_adata_epi = df_adata_epi[df_adata_epi['names'].isin(adata_epi.var_names[adata_epi.var.n_cells_by_counts>30])]

df_adata_epi


# In[88]:


gene_set_names = gseapy.get_library_name(organism='Human')


# In[89]:


res = gseapy.prerank(rnk=df_adata_epi, gene_sets='GO_Biological_Process_2023', min_size =10)
a = pd.DataFrame(res.res2d)
a.sort_values(by=['NOM p-val'])[:10]


# In[90]:


a[a.Term=='Immune Response-Regulating Cell Surface Receptor Signaling Pathway (GO:0002768)']


# In[63]:


ax=gseapy.gseaplot(rank_metric=res.ranking, term=terms[0], **res.results[terms[0]])

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/GSEA_epi.pdf',bbox_inches="tight")


# In[64]:


adata_T = adata[adata.obs.new_leiden =='T']
sc.tl.rank_genes_groups(adata_T, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
df_adata_T = sc.get.rank_genes_groups_df(adata_T, group='LOW',key = "rank_genes_groups", pval_cutoff=0.05)[['names','logfoldchanges']]
df_adata_T.sort_values(by=['logfoldchanges'], inplace=True, ascending=False)

# calculate_qc_metrics will calculate number of cells per gene
sc.pp.calculate_qc_metrics(adata_T, percent_top=None, log1p=False, inplace=True)

# filter for genes expressed in at least 30 cells.
df_adata_T = df_adata_T[df_adata_T['names'].isin(adata_T.var_names[adata_T.var.n_cells_by_counts>30])]

df_adata_T


# In[27]:


gene_set_names = gseapy.get_library_name(organism='Human')


# In[65]:


res = gseapy.prerank(rnk=df_adata_T, gene_sets='GO_Biological_Process_2023', min_size =10)
a = pd.DataFrame(res.res2d)
a.sort_values(by=['NOM p-val'])


# In[66]:


terms = res.res2d.Term
terms[:10]


# In[67]:


ax=gseapy.gseaplot(rank_metric=res.ranking, term=terms[2], **res.results[terms[2]])

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/GSEA_T.pdf',bbox_inches="tight")


# ### Identify iNKT cell in fov 38, with high expression of IL10

# In[30]:


sq.pl.spatial_segment(
    adata,
    color=['IL10','new_leiden'],
    library_key="fov",
    library_id="38",
    title = ['fov 38 - IL10 level', 'fov 38 - Cells'],
    seg_cell_id="cell_ID",
    img=False, palette=mp.colors.ListedColormap(palette)
)
plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/fov38_IL10_all_cells.pdf',bbox_inches="tight")


# In[19]:


adata_T = adata[adata.obs.new_leiden == 'T']
adata_T.obs['cell_id'] = adata_T.obs.index.values
adata_T.obs['cell_id'] = adata_T.obs['cell_id'].astype("category")


# In[32]:


sq.pl.spatial_segment(
    adata_T,
    color=['IL10'],
    library_key="fov",
    library_id="38",
    seg_cell_id="cell_ID",
    img=False,
    size=60
)
plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/fov38_IL10_Tcells.pdf',bbox_inches="tight")


# In[33]:


adata_T[adata_T.obs.fov=='38'].obs['CenterY_global_px'].sort_values()


# In[20]:


adata_fov38 = adata[adata.obs['fov'] == '38']
#adata_fov38.obs['cell_id'] = adata_fov38.obs.index.values + '-' +np.array(adata_fov38.obs.new_leiden)
#adata_fov38.obs['cell_id'] = adata_fov38.obs['cell_id'].astype("category")
#adata_fov38.obs['cell_id'] = adata_fov38.obs['cell_id'].astype("category")
adata_fov38_T_epi_fibro = adata_fov38[adata_fov38.obs['new_leiden'].isin(['T','Epithelial','Fibroblast'])]


# In[21]:


adata_2274_38 = adata_T[adata_T.obs['cell_id'] == '2274_38']
sq.pl.spatial_segment(
    adata_2274_38,
    color= 'IL10',
    library_key="fov",
    library_id="38",
    seg_cell_id="cell_ID",
    img=False,
    size=60)


# In[34]:


a = pd.DataFrame(adata_2274_38.X.toarray(),index = adata_2274_38.obs.index.values, columns = adata_2274_38.var.index.values)
a['ZBTB16']


# In[36]:


pos = adata_fov38_T_epi_fibro.obs[['CenterX_global_px','CenterY_global_px']]
pos


# In[37]:


# matrix 18 x 18 with the distance between one cell and each other
dm = distance_matrix(pos, pos)


# In[38]:


dm = pd.DataFrame(dm, index = pos.index.values, columns = pos.index.values)
dm


# In[103]:


adata.obs.loc['2274_38']
# T


# In[31]:


dm_T_2274_38 = dm.loc['2274_38']
dm_T_2274_38[dm_T_2274_38 < 500].index.values


# In[34]:


distance = []

for c in adata_fov38_T_epi_fibro.obs.index.values:
    if c in dm_T_2274_38[dm_T_2274_38 < 500].index.values:
        distance.append('within_500')
    else:
        distance.append('farther_500')

adata_fov38_T_epi_fibro.obs['distance'] = distance
adata_fov38_T_epi_fibro.obs['distance'] = adata_fov38_T_epi_fibro.obs['distance'].astype("category")


# In[42]:


adata_fov38_T = adata_fov38_T_epi_fibro[adata_fov38_T_epi_fibro.obs.new_leiden == 'T']
adata_fov38_epithelial = adata_fov38_T_epi_fibro[adata_fov38_T_epi_fibro.obs.new_leiden == 'Epithelial']
adata_fov38_fibroblast = adata_fov38_T_epi_fibro[adata_fov38_T_epi_fibro.obs.new_leiden == 'Fibroblast']


# In[88]:


sc.tl.rank_genes_groups(adata_fov38_T, 'distance', groups=['within_500'], reference='farther_500', method='wilcoxon')
sc.pl.rank_genes_groups(adata_fov38_T, groups=['within_500'], n_genes=20, save='distance_wilcoxon_T_fov38_within_vs_farther.pdf')
sc.tl.rank_genes_groups(adata_fov38_T, 'distance', groups=['farther_500'], reference='within_500', method='wilcoxon')
sc.pl.rank_genes_groups(adata_fov38_T, groups=['farther_500'], n_genes=20, save='distance_wilcoxon_T_fov38_farther_vs_within.pdf')

sc.tl.rank_genes_groups(adata_fov38_epithelial, 'distance', groups=['within_500'], reference='farther_500', method='wilcoxon')
sc.pl.rank_genes_groups(adata_fov38_epithelial, groups=['within_500'], n_genes=20, save='distance_wilcoxon_epithelial_fov38_within_vs_farther.pdf')
sc.tl.rank_genes_groups(adata_fov38_epithelial, 'distance', groups=['farther_500'], reference='within_500', method='wilcoxon')
sc.pl.rank_genes_groups(adata_fov38_epithelial, groups=['farther_500'], n_genes=20, save='distance_wilcoxon_epithelial_fov38_farther_vs_within.pdf')


# In[84]:


df_T_fov38 = sc.get.rank_genes_groups_df(adata_fov38_T, group='within_500',key = "wilcoxon")
df_T_fov38.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_T_within_further_500_iNKT10_fov38.xlsx')


# In[85]:


df_epithelial_fov38 = sc.get.rank_genes_groups_df(adata_fov38_epithelial, group='within_500',key = "wilcoxon")
df_epithelial_fov38.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_epithelial_within_further_500_iNKT10_fov38.xlsx')


# ### DEGs between cells of the same cell type, high vs low iNKT10 level

# In[35]:


adata_B = adata[adata.obs.new_leiden =='B']
adata_Breg = adata[adata.obs.new_leiden =='B reg']
adata_endo = adata[adata.obs.new_leiden =='Endothelial']
adata_glia = adata[adata.obs.new_leiden =='Enteric glia cells']
adata_epi = adata[adata.obs.new_leiden =='Epithelial']
adata_fibro = adata[adata.obs.new_leiden =='Fibroblast']
adata_macro = adata[adata.obs.new_leiden =='Macrophage']
adata_mast = adata[adata.obs.new_leiden =='Mast cells']
adata_plasma_A = adata[adata.obs.new_leiden =='Plasma IgA']
adata_plasma_G = adata[adata.obs.new_leiden =='Plasma IgG']
adata_muscle = adata[adata.obs.new_leiden =='Smooth muscle cells']
adata_T = adata[adata.obs.new_leiden =='T']


# In[40]:


sc.tl.rank_genes_groups(adata_B, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
df_adata_B = sc.get.rank_genes_groups_df(adata_B, group='LOW',key = "rank_genes_groups")
df_adata_B.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_B_low_high.xlsx')

sc.tl.rank_genes_groups(adata_Breg, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
df_adata_B_reg = sc.get.rank_genes_groups_df(adata_Breg, group='LOW',key = "rank_genes_groups")
df_adata_B_reg.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_B_reg_low_high.xlsx')

sc.tl.rank_genes_groups(adata_endo, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
df_adata_endo = sc.get.rank_genes_groups_df(adata_endo, group='LOW',key = "rank_genes_groups")
df_adata_endo.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_endo_low_high.xlsx')

sc.tl.rank_genes_groups(adata_glia, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
df_adata_glia = sc.get.rank_genes_groups_df(adata_glia, group='LOW',key = "rank_genes_groups")
df_adata_glia.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_glia_low_high.xlsx')

sc.tl.rank_genes_groups(adata_epi, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
df_adata_epi = sc.get.rank_genes_groups_df(adata_epi, group='LOW',key = "rank_genes_groups")
df_adata_epi.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_epi_low_high.xlsx')

sc.tl.rank_genes_groups(adata_fibro, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
df_adata_fibro = sc.get.rank_genes_groups_df(adata_fibro, group='LOW',key = "rank_genes_groups")
df_adata_fibro.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_fibro_low_high.xlsx')

sc.tl.rank_genes_groups(adata_macro, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
df_adata_macro = sc.get.rank_genes_groups_df(adata_macro, group='LOW',key = "rank_genes_groups")
df_adata_macro.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_macro_low_high.xlsx')

sc.tl.rank_genes_groups(adata_mast, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
df_adata_mast = sc.get.rank_genes_groups_df(adata_mast, group='LOW',key = "rank_genes_groups")
df_adata_mast.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_mast_low_high.xlsx')

sc.tl.rank_genes_groups(adata_plasma_A, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
df_adata_plasma_A = sc.get.rank_genes_groups_df(adata_plasma_A, group='LOW',key = "rank_genes_groups")
df_adata_plasma_A.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_plasma_A_low_high.xlsx')

sc.tl.rank_genes_groups(adata_plasma_G, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
df_adata_plasma_G = sc.get.rank_genes_groups_df(adata_plasma_G, group='LOW',key = "rank_genes_groups")
df_adata_plasma_G.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_plasma_G_low_high.xlsx')

sc.tl.rank_genes_groups(adata_muscle, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
df_adata_muscle = sc.get.rank_genes_groups_df(adata_muscle, group='LOW',key = "rank_genes_groups")
df_adata_muscle.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_muscle_low_high.xlsx')

sc.tl.rank_genes_groups(adata_T, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
df_adata_T = sc.get.rank_genes_groups_df(adata_T, group='LOW',key = "rank_genes_groups")
df_adata_T.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_T_low_high.xlsx')


# ### Identification of T cells with high expression of ZBTB16 and IL10 
# Identification of epithelial and fibroblasts within 500 micro from these iNKT cells

# In[18]:


gene_expression_all_cells_T_cells = adata[adata.obs.new_leiden == 'T'].obs.index.values
gene_expression_all_cells = pd.DataFrame(adata.X.toarray(),index = adata.obs.index.values, columns = adata.var.index.values)

gene_expression_T_cells= gene_expression_all_cells.loc[gene_expression_all_cells_T_cells]
gene_expression_T_cells


# In[19]:


gene_expression_T_cells['IL10'][gene_expression_T_cells['IL10'] > 0.1]
gene_expression_T_cells['IL10'].max()


# In[20]:


adata_low = adata[adata.obs.level == 'LOW']
adata_high = adata[adata.obs.level == 'HIGH']


# In[21]:


# for each fov, compute how many T cells express ZBTB16 > 0.1 and how many of them express IL10 > 0.1 

fov_ZBTB16_low = {}

for fov in pd.unique(adata_low.obs.fov):
    adata_fov = adata_low[adata_low.obs.fov == fov]
    index_T_cells = adata_fov[adata_fov.obs.new_leiden == 'T'].obs.index.values
    gene_expression_fov_all_cells = pd.DataFrame(adata_fov.X.toarray(),index = adata_fov.obs.index.values, columns = adata_fov.var.index.values)
    gene_expression_fov_T = gene_expression_fov_all_cells.loc[index_T_cells]
    T_cells_ZBTB16 = gene_expression_fov_T['ZBTB16'][gene_expression_fov_T['ZBTB16'] > 0.1].index.values
    gene_expression_fov_T_ZBTB16 = gene_expression_fov_all_cells.loc[T_cells_ZBTB16]
    T_cells_ZBTB16_IL10 = gene_expression_fov_T_ZBTB16['IL10'][gene_expression_fov_T['IL10'] > 0.1].index.values

    fov_ZBTB16_low[fov] = [len(T_cells_ZBTB16), len(T_cells_ZBTB16_IL10)]


# In[55]:


fov_ZBTB16_low


# In[22]:


# for each fov, compute how many T cells express ZBTB16 > 0.1 and how many of them express IL10 > 0.1 

fov_ZBTB16_high = {}

for fov in pd.unique(adata_high.obs.fov):
    adata_fov = adata_high[adata_high.obs.fov == fov]
    index_T_cells = adata_fov[adata_fov.obs.new_leiden == 'T'].obs.index.values
    gene_expression_fov_all_cells = pd.DataFrame(adata_fov.X.toarray(),index = adata_fov.obs.index.values, columns = adata_fov.var.index.values)
    gene_expression_fov_T = gene_expression_fov_all_cells.loc[index_T_cells]
    T_cells_ZBTB16 = gene_expression_fov_T['ZBTB16'][gene_expression_fov_T['ZBTB16'] > 0.1].index.values
    gene_expression_fov_T_ZBTB16 = gene_expression_fov_all_cells.loc[T_cells_ZBTB16]
    T_cells_ZBTB16_IL10 = gene_expression_fov_T_ZBTB16['IL10'][gene_expression_fov_T['IL10'] > 0.1].index.values

    fov_ZBTB16_high[fov] = [len(T_cells_ZBTB16), len(T_cells_ZBTB16_IL10)]


# In[57]:


fov_ZBTB16_high


# In[23]:


# for each fov, high and low together, identify the percentage of epithelial/fibroblast cells closer that 500 micro from T ZBTB16+IL10+ cells/ZBTB16+/IL10+ 
# and the percentage of the epithelial/fibroblast cells farther that 500 micro

def percentage_cells(distance):
    epithelial_less_ZBTB16_IL10 = []
    fibroblast_less_ZBTB16_IL10 = []
    epithelial_more_ZBTB16_IL10 = []
    fibroblast_more_ZBTB16_IL10 = []
    
    epithelial_less_ZBTB16 = []
    fibroblast_less_ZBTB16 = []
    
    epithelial_less_IL10 = []
    fibroblast_less_IL10 = []
    
    for fov in pd.unique(adata.obs.fov):
        adata_fov = adata[adata.obs.fov == fov]
        index_T_cells = adata_fov[adata_fov.obs.new_leiden == 'T'].obs.index.values
        gene_expression_fov_all_cells = pd.DataFrame(adata_fov.X.toarray(),index = adata_fov.obs.index.values, columns = adata_fov.var.index.values)
        gene_expression_fov_T = gene_expression_fov_all_cells.loc[index_T_cells]
        T_cells_ZBTB16 = gene_expression_fov_T['ZBTB16'][gene_expression_fov_T['ZBTB16'] > 0.1].index.values
        gene_expression_fov_T_ZBTB16 = gene_expression_fov_all_cells.loc[T_cells_ZBTB16]
        T_cells_ZBTB16_IL10 = gene_expression_fov_T_ZBTB16['IL10'][gene_expression_fov_T['IL10'] > 0.1].index.values
        T_cells_IL10 = gene_expression_fov_T['IL10'][gene_expression_fov_T['IL10'] > 0.1].index.values
    
        pos = adata_fov.obs[['CenterX_global_px','CenterY_global_px']]
        dm = pd.DataFrame(distance_matrix(pos, pos), index = pos.index.values, columns = pos.index.values)

        epithelial_less_ZBTB16_IL10_fov = []
        fibroblast_less_ZBTB16_IL10_fov = []
        epithelial_less_ZBTB16_fov = []
        fibroblast_less_ZBTB16_fov = []
        epithelial_less_IL10_fov = []
        fibroblast_less_IL10_fov = []
        
        for cell_ZBTB16_IL10 in T_cells_ZBTB16_IL10:
            dm_cell = dm[cell_ZBTB16_IL10]
            dm_cell_less = dm_cell[dm_cell <= distance].index.values
            dm_cell_more = dm_cell[dm_cell > distance].index.values

            epithelial_less = adata_fov.obs.loc[dm_cell_less].new_leiden[adata_fov.obs.loc[dm_cell_less].new_leiden == 'Epithelial'].index.values
            if len(epithelial_less) > 0:
                epithelial_less_ZBTB16_IL10_fov.extend(epithelial_less)
            fibroblast_less = adata_fov.obs.loc[dm_cell_less].new_leiden[adata_fov.obs.loc[dm_cell_less].new_leiden == 'Fibroblast'].index.values
            if len(fibroblast_less) > 0:    
                fibroblast_less_ZBTB16_IL10_fov.extend(fibroblast_less)
            
        for cell_ZBTB16 in T_cells_ZBTB16:
            dm_cell = dm[cell_ZBTB16]
            dm_cell_less = dm_cell[dm_cell <= distance].index.values
    
            epithelial_less = adata_fov.obs.loc[dm_cell_less].new_leiden[adata_fov.obs.loc[dm_cell_less].new_leiden == 'Epithelial'].index.values
            if len(epithelial_less) > 0:
                epithelial_less_ZBTB16_fov.extend(epithelial_less)
            fibroblast_less = adata_fov.obs.loc[dm_cell_less].new_leiden[adata_fov.obs.loc[dm_cell_less].new_leiden == 'Fibroblast'].index.values
            if len(fibroblast_less) > 0:    
                fibroblast_less_ZBTB16_fov.extend(fibroblast_less)

        for cell_IL10 in T_cells_IL10:
            dm_cell = dm[cell_IL10]
            dm_cell_less = dm_cell[dm_cell <= distance].index.values
    
            epithelial_less = adata_fov.obs.loc[dm_cell_less].new_leiden[adata_fov.obs.loc[dm_cell_less].new_leiden == 'Epithelial'].index.values
            if len(epithelial_less) > 0:
                epithelial_less_IL10_fov.extend(epithelial_less)
            fibroblast_less = adata_fov.obs.loc[dm_cell_less].new_leiden[adata_fov.obs.loc[dm_cell_less].new_leiden == 'Fibroblast'].index.values
            if len(fibroblast_less) > 0:    
                fibroblast_less_IL10_fov.extend(fibroblast_less)

        epithelial_more_ZBTB16_IL10_fov = []
        fibroblast_more_ZBTB16_IL10_fov = []
    
        index_epithelial_cells = adata_fov[adata_fov.obs.new_leiden == 'Epithelial'].obs.index.values
        index_fibroblast_cells = adata_fov[adata_fov.obs.new_leiden == 'Fibroblast'].obs.index.values

        for c in index_epithelial_cells:
             if c not in np.unique(epithelial_less_ZBTB16_IL10_fov):
                 epithelial_more_ZBTB16_IL10_fov.append(c)
            
        for c in index_fibroblast_cells:
             if c not in np.unique(fibroblast_less_ZBTB16_IL10_fov):
                 fibroblast_more_ZBTB16_IL10_fov.append(c)
    
        if len(adata_fov[adata_fov.obs.new_leiden == 'Epithelial']) != 0:
             epithelial_less_ZBTB16_IL10.append(len(np.unique(epithelial_less_ZBTB16_IL10_fov))*100/len(adata_fov[adata_fov.obs.new_leiden == 'Epithelial']))
             epithelial_more_ZBTB16_IL10.append(len(np.unique(epithelial_more_ZBTB16_IL10_fov))*100/len(adata_fov[adata_fov.obs.new_leiden == 'Epithelial']))
             epithelial_less_ZBTB16.append(len(np.unique(epithelial_less_ZBTB16_fov))*100/len(adata_fov[adata_fov.obs.new_leiden == 'Epithelial']))
             epithelial_less_IL10.append(len(np.unique(epithelial_less_IL10_fov))*100/len(adata_fov[adata_fov.obs.new_leiden == 'Epithelial']))
        
        if len(adata_fov[adata_fov.obs.new_leiden == 'Fibroblast']) != 0:
            fibroblast_less_ZBTB16_IL10.append(len(np.unique(fibroblast_less_ZBTB16_IL10_fov))*100/len(adata_fov[adata_fov.obs.new_leiden == 'Fibroblast']))
            fibroblast_more_ZBTB16_IL10.append(len(np.unique(fibroblast_more_ZBTB16_IL10_fov))*100/len(adata_fov[adata_fov.obs.new_leiden == 'Fibroblast']))
            fibroblast_less_ZBTB16.append(len(np.unique(fibroblast_less_ZBTB16_fov))*100/len(adata_fov[adata_fov.obs.new_leiden == 'Fibroblast']))
            fibroblast_less_IL10.append(len(np.unique(fibroblast_less_IL10_fov))*100/len(adata_fov[adata_fov.obs.new_leiden == 'Fibroblast']))

    with open("/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/number_epithelial_less_" + str(distance) + "_ZBTB16_IL10", "w") as file:
        for row in epithelial_less_ZBTB16_IL10:
                s = "".join(map(str, str(row)))
                file.write(s+'\n')
    
    with open("/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/number_fibroblast_less_" + str(distance) + "_ZBTB16_IL10", "w") as file:
        for row in fibroblast_less_ZBTB16_IL10:
                s = "".join(map(str, str(row)))
                file.write(s+'\n')

    with open("/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/number_epithelial_more_" + str(distance) + "_ZBTB16_IL10", "w") as file:
        for row in epithelial_more_ZBTB16_IL10:
                s = "".join(map(str, str(row)))
                file.write(s+'\n')

    with open("/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/number_fibroblast_more_" + str(distance) + "_ZBTB16_IL10", "w") as file:
        for row in fibroblast_more_ZBTB16_IL10:
                s = "".join(map(str, str(row)))
                file.write(s+'\n')
            
    with open("/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/number_epithelial_less_" + str(distance) + "_ZBTB16", "w") as file:
        for row in epithelial_less_ZBTB16:
                s = "".join(map(str, str(row)))
                file.write(s+'\n')
            
    with open("/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/number_fibroblast_less_" + str(distance) + "_ZBTB16", "w") as file:
        for row in fibroblast_less_ZBTB16:
                s = "".join(map(str, str(row)))
                file.write(s+'\n')

    with open("/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/number_epithelial_less_" + str(distance) + "_IL10", "w") as file:
        for row in epithelial_less_IL10:
                s = "".join(map(str, str(row)))
                file.write(s+'\n')

    with open("/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/number_fibroblast_less_" + str(distance) + "_IL10", "w") as file:
        for row in fibroblast_less_IL10:
                s = "".join(map(str, str(row)))
                file.write(s+'\n')
    


# In[24]:


percentage_cells(100)
percentage_cells(500)


# In[25]:


# for each fov, high and low together, identify the id of the epithelial/fibroblast cells closer that 500 micro from T ZBTB16+IL10+ cells/ZBTB16+/IL10+ 
# and the id of the epithelial/fibroblast cells farther than 500 micro

def id_cells_within(distance):
    epithelial_less_ZBTB16_IL10 = []
    fibroblast_less_ZBTB16_IL10 = []
    
    epithelial_less_ZBTB16 = []
    fibroblast_less_ZBTB16 = []

    epithelial_less_IL10 = []
    fibroblast_less_IL10 = []
    
    for fov in pd.unique(adata.obs.fov):
        adata_fov = adata[adata.obs.fov == fov]
        index_T_cells = adata_fov[adata_fov.obs.new_leiden == 'T'].obs.index.values
        gene_expression_fov_all_cells = pd.DataFrame(adata_fov.X.toarray(),index = adata_fov.obs.index.values, columns = adata_fov.var.index.values)
        gene_expression_fov_T = gene_expression_fov_all_cells.loc[index_T_cells]
        T_cells_ZBTB16 = gene_expression_fov_T['ZBTB16'][gene_expression_fov_T['ZBTB16'] > 0.1].index.values
        gene_expression_fov_T_ZBTB16 = gene_expression_fov_all_cells.loc[T_cells_ZBTB16]
        T_cells_ZBTB16_IL10 = gene_expression_fov_T_ZBTB16['IL10'][gene_expression_fov_T['IL10'] > 0.1].index.values
        T_cells_IL10 = gene_expression_fov_T['IL10'][gene_expression_fov_T['IL10'] > 0.1].index.values

        T_cells_IL10 = [x for x in T_cells_IL10 if x not in T_cells_ZBTB16_IL10]
        T_cells_ZBTB16 = [x for x in T_cells_ZBTB16 if x not in T_cells_ZBTB16_IL10]
        
        pos = adata_fov.obs[['CenterX_global_px','CenterY_global_px']]
        dm = pd.DataFrame(distance_matrix(pos, pos), index = pos.index.values, columns = pos.index.values)

        for cell_ZBTB16_IL10 in T_cells_ZBTB16_IL10:
            adata.obs.loc[cell_ZBTB16_IL10,'within_'+str(distance)] = 'T ZBTB16+ IL10+'
            dm_cell = dm[cell_ZBTB16_IL10]
            dm_cell_less = dm_cell[dm_cell <= distance].index.values
            dm_cell_more = dm_cell[dm_cell > distance].index.values
    
            epithelial_less = adata_fov.obs.loc[dm_cell_less].new_leiden[adata_fov.obs.loc[dm_cell_less].new_leiden == 'Epithelial'].index.values
            if len(epithelial_less) > 0:
                epithelial_less_ZBTB16_IL10.extend(epithelial_less)
            fibroblast_less = adata_fov.obs.loc[dm_cell_less].new_leiden[adata_fov.obs.loc[dm_cell_less].new_leiden == 'Fibroblast'].index.values
            if len(fibroblast_less) > 0:    
                fibroblast_less_ZBTB16_IL10.extend(fibroblast_less)

        for cell_ZBTB16 in T_cells_ZBTB16:
            adata.obs.loc[cell_ZBTB16,'within_'+str(distance)] = 'T ZBTB16+'
            dm_cell = dm[cell_ZBTB16]
            dm_cell_less = dm_cell[dm_cell <= distance].index.values
    
            epithelial_less = adata_fov.obs.loc[dm_cell_less].new_leiden[adata_fov.obs.loc[dm_cell_less].new_leiden == 'Epithelial'].index.values
            if len(epithelial_less) > 0:
                epithelial_less_ZBTB16.extend(epithelial_less)
            fibroblast_less = adata_fov.obs.loc[dm_cell_less].new_leiden[adata_fov.obs.loc[dm_cell_less].new_leiden == 'Fibroblast'].index.values
            if len(fibroblast_less) > 0:    
                fibroblast_less_ZBTB16.extend(fibroblast_less)


        for cell_IL10 in T_cells_IL10:
            adata.obs.loc[cell_IL10, 'within_'+str(distance)] = 'T IL10+'
            dm_cell = dm[cell_IL10]
            dm_cell_less = dm_cell[dm_cell <= distance].index.values
    
            epithelial_less = adata_fov.obs.loc[dm_cell_less].new_leiden[adata_fov.obs.loc[dm_cell_less].new_leiden == 'Epithelial'].index.values
            if len(epithelial_less) > 0:
                epithelial_less_IL10.extend(epithelial_less)
            fibroblast_less = adata_fov.obs.loc[dm_cell_less].new_leiden[adata_fov.obs.loc[dm_cell_less].new_leiden == 'Fibroblast'].index.values
            if len(fibroblast_less) > 0:    
                fibroblast_less_IL10.extend(fibroblast_less)
            
    epithelial_less_ZBTB16 = [x for x in epithelial_less_ZBTB16 if x not in epithelial_less_ZBTB16_IL10]
    epithelial_less_IL10 = [x for x in epithelial_less_IL10 if x not in epithelial_less_ZBTB16_IL10]
    fibroblast_less_ZBTB16 = [x for x in fibroblast_less_ZBTB16 if x not in fibroblast_less_ZBTB16_IL10]
    fibroblast_less_IL10 = [x for x in fibroblast_less_IL10 if x not in fibroblast_less_ZBTB16_IL10]
    
    for c in epithelial_less_ZBTB16_IL10:
        adata.obs.loc[c,'within_'+str(distance)] = 'epithelial_within_'+str(distance)+'_T_ZBTB16_IL10'
        adata.obs.loc[c,'epithelial_'+str(distance) +'_ZBTB16_IL10'] = 'within'

    for c in epithelial_less_ZBTB16:
        adata.obs.loc[c,'within_'+str(distance)]= 'epithelial_within_'+str(distance)+'_T_ZBTB16'
    
    for c in epithelial_less_IL10:
        adata.obs.loc[c,'within_'+str(distance)] = 'epithelial_within_'+str(distance)+'_T_IL10'
    
    for c in fibroblast_less_ZBTB16_IL10:
        adata.obs.loc[c,'within_'+str(distance)] = 'fibroblast_within_'+str(distance)+'_T_ZBTB16_IL10'
        adata.obs.loc[c,'fibroblast_'+str(distance) +'_ZBTB16_IL10'] = 'within'
    
    for c in fibroblast_less_ZBTB16:
        adata.obs.loc[c,'within_'+str(distance)] = 'fibroblast_within_'+str(distance)+'_T_ZBTB16'
    
    for c in fibroblast_less_IL10:
        adata.obs.loc[c,'within_'+str(distance)] = 'fibroblast_within_'+str(distance)+'_T_IL10'
        
    epithelial_more_ZBTB16_IL10 = []
    fibroblast_more_ZBTB16_IL10 = []
    
    for fov in pd.unique(adata.obs.fov):
        adata_fov = adata[adata.obs.fov == fov]
        index_epithelial_cells = adata_fov[adata_fov.obs.new_leiden == 'Epithelial'].obs.index.values
        index_fibroblast_cells = adata_fov[adata_fov.obs.new_leiden == 'Fibroblast'].obs.index.values
        
        for c in index_epithelial_cells:
            if c not in np.unique(epithelial_less_ZBTB16_IL10):
                epithelial_more_ZBTB16_IL10.append(c)

        for c in index_fibroblast_cells:
            if c not in np.unique(fibroblast_less_ZBTB16_IL10):
                fibroblast_more_ZBTB16_IL10.append(c)    

    for c in epithelial_more_ZBTB16_IL10:
        adata.obs.loc[c,'epithelial_'+str(distance) +'_ZBTB16_IL10'] = 'farther'
    
    for c in fibroblast_more_ZBTB16_IL10:
        adata.obs.loc[c,'fibroblast_'+str(distance) +'_ZBTB16_IL10'] = 'farther'


# In[26]:


id_cells_within(100)
id_cells_within(500)


# In[62]:


adata_sub_100 = adata[adata.obs['within_100'].notna()]
adata_sub_500 = adata[adata.obs['within_500'].notna()]


# In[63]:


np.array(pd.unique(adata_sub_100.obs['fov']))


# In[64]:


t_100 = ['fov 21 - Patient_11 - LOW',
 'fov 23 - Patient_12 - HIGH',
 'fov 24 - Patient_12 - HIGH',
 'fov 25 - Patient_13 - LOW',
 'fov 27 - Patient_14 - HIGH',
 'fov 28 - Patient_14 - HIGH',
 'fov 29 - Patient_15 - LOW',
 'fov 30 - Patient_15 - LOW',
 'fov 31 - Patient_16 - HIGH',
 'fov 32 - Patient_16 - HIGH',
 'fov 33 - Patient_17 - HIGH',
 'fov 34 - Patient_17 - HIGH',
 'fov 35 - Patient_18 - HIGH',
 'fov 36 - Patient_18 - HIGH',
 'fov 37 - Patient_19 - LOW',
 'fov 38 - Patient_19 - LOW',
 'fov 39 - Patient_20 - LOW',
 'fov 40 - Patient_20 - LOW',
 'fov 41 - Patient_20 - LOW']


# In[65]:


np.array(pd.unique(adata_sub_500.obs['fov']))


# In[66]:


t_500 = ['fov 21 - Patient_11 - LOW',
 'fov 23 - Patient_12 - HIGH',
 'fov 24 - Patient_12 - HIGH',
 'fov 25 - Patient_13 - LOW',
 'fov 27 - Patient_14 - HIGH',
 'fov 28 - Patient_14 - HIGH',
 'fov 29 - Patient_15 - LOW',
 'fov 30 - Patient_15 - LOW',
 'fov 31 - Patient_16 - HIGH',
 'fov 32 - Patient_16 - HIGH',
 'fov 33 - Patient_17 - HIGH',
 'fov 34 - Patient_17 - HIGH',
 'fov 35 - Patient_18 - HIGH',
 'fov 36 - Patient_18 - HIGH',
 'fov 37 - Patient_19 - LOW',
 'fov 38 - Patient_19 - LOW',        
 'fov 39 - Patient_20 - LOW',
 'fov 40 - Patient_20 - LOW',
 'fov 41 - Patient_20 - LOW']


# In[67]:


pd.unique(adata_sub_100.obs['within_100'])


# In[68]:


palette = ['#F0E442','#000000', '#56B4E9', '#D55E00',  '#0072B2', '#CC79A7', '#F7B6D2', '#7F7F7F']

sns.color_palette(palette)


# In[69]:


sc.pl.umap(
    adata_sub_100, show=False,
    color= "within_100", palette=sns.color_palette(palette, 8)
)


# In[70]:


sq.pl.spatial_segment(
    adata_sub_100,
    color='within_100',
    library_key="fov",
    seg_cell_id="cell_ID",
    seg_outline=True,
    img=False,
    scalebar_dx=1.0,
    title = t_100, 
    scalebar_kwargs={"scale_loc": "bottom", "location": "lower right"}
)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/fov_patients_within100_senza_back.pdf',bbox_inches="tight")


# In[71]:


pd.unique(adata_sub_500.obs['within_500'])


# In[72]:


sq.pl.spatial_segment(
    adata_sub_500,
    color='within_500',
    library_key="fov",
    seg_cell_id="cell_ID",
    seg_outline=True,
    img=False,
    scalebar_dx=1.0,
    title = t_500, 
    scalebar_kwargs={"scale_loc": "bottom", "location": "lower right"}
)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/fov_patients_within500_senza_back.pdf',bbox_inches="tight")


# ### DEGs between epithelial/fibro closer/farther than 500 micro from T ZBTB16+IL10+ cells considering high and low cells together

# In[29]:


adata[adata.obs['epithelial_500_ZBTB16_IL10'] == 'farther']


# In[30]:


adata[adata.obs['fibroblast_500_ZBTB16_IL10'] == 'farther']


# In[77]:


sc.tl.rank_genes_groups(adata, 'epithelial_500_ZBTB16_IL10', groups=['within'], reference='farther', method='wilcoxon')
df_adata_epi_500 = sc.get.rank_genes_groups_df(adata, group='within',key = "rank_genes_groups")
df_adata_epi_500.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_epithelial_within_farther_500_ZBTB16_IL10.xlsx')

sc.tl.rank_genes_groups(adata, 'fibroblast_500_ZBTB16_IL10', groups=['within'], reference='farther', method='wilcoxon')
df_adata_fibro_500 = sc.get.rank_genes_groups_df(adata, group='within',key = "rank_genes_groups")
df_adata_fibro_500.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_fibroblast_within_farther_500_ZBTB16_IL10.xlsx')


# ### DEGs between epithelial cells/fibro within 500 micro from T ZBTB16+IL10 and epithelial cells within 500 micro from T ZBTB16

# In[78]:


pd.unique(adata_sub_500.obs['within_500'])


# In[92]:


adata[adata.obs['within_500'] == 'epithelial_within_500_T_ZBTB16']


# In[93]:


adata[adata.obs['within_500'] == 'fibroblast_within_500_T_ZBTB16']


# In[87]:


sc.tl.rank_genes_groups(adata, 'within_500', groups=['epithelial_within_500_T_ZBTB16_IL10'], reference='epithelial_within_500_T_ZBTB16', method='wilcoxon')
df_adata_500 = sc.get.rank_genes_groups_df(adata, group='epithelial_within_500_T_ZBTB16_IL10',key = "rank_genes_groups")
df_adata_500.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_epithelial_within_500_ZBTB16_IL10_ZBTB16.xlsx')

sc.tl.rank_genes_groups(adata, 'within_500', groups=['fibroblast_within_500_T_ZBTB16_IL10'], reference='fibroblast_within_500_T_ZBTB16', method='wilcoxon')
df_adata_500 = sc.get.rank_genes_groups_df(adata, group='fibroblast_within_500_T_ZBTB16_IL10',key = "rank_genes_groups")
df_adata_500.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_fibroblast_within_500_ZBTB16_IL10_ZBTB16.xlsx')


# ### DEGs between epithelial cells/fibroblasts within 500 micro from T ZBTB16+IL10 and epithelial cells within 500 micro from T IL10

# In[95]:


adata[adata.obs['within_500'] == 'epithelial_within_500_T_IL10']


# In[97]:


adata[adata.obs['within_500'] == 'fibroblast_within_500_T_IL10']


# In[ ]:


sc.tl.rank_genes_groups(adata, 'within_500', groups=['epithelial_within_500_T_ZBTB16_IL10'], reference='epithelial_within_500_T_IL10', method='wilcoxon')
df_adata_500 = sc.get.rank_genes_groups_df(adata, group='epithelial_within_500_T_ZBTB16_IL10',key = "rank_genes_groups")
df_adata_500.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_epithelial_within_500_ZBTB16_IL10_IL10.xlsx')

sc.tl.rank_genes_groups(adata, 'within_500', groups=['fibroblast_within_500_T_ZBTB16_IL10'], reference='fibroblast_within_500_T_IL10', method='wilcoxon')
df_adata_500 = sc.get.rank_genes_groups_df(adata, group='fibroblast_within_500_T_ZBTB16_IL10',key = "rank_genes_groups")
df_adata_500.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_fibroblast_within_500_ZBTB16_IL10_IL10.xlsx')


# ### Compare distance between cell types in high and low iNKT10 conditions

# In[31]:


sc.pl.umap(
    adata_low, show=False,
    color= "new_leiden", palette=sns.color_palette(palette),
    title = 'Low iNKT10'
)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/umap_clusters_named_low.pdf',bbox_inches="tight")


# In[32]:


sc.pl.umap(
    adata_high, show=False,
    color= "new_leiden", palette=sns.color_palette(palette),
    title = 'High iNKT10'   
)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/umap_clusters_named_high.pdf',bbox_inches="tight")


# In[16]:


sq.gr.spatial_neighbors(adata_low, coord_type="generic", spatial_key="spatial")
sq.gr.nhood_enrichment(adata_low, cluster_key="new_leiden")
sq.gr.spatial_neighbors(adata_high, coord_type="generic", spatial_key="spatial")
sq.gr.nhood_enrichment(adata_high, cluster_key="new_leiden")


# In[17]:


# plotted values are in adata_low.uns['leiden_nhood_enrichment']['zscore']

sq.pl.nhood_enrichment(
    adata_low,
    cluster_key="new_leiden",
    method="average",
    cmap='OrRd',
    figsize=(5, 5),
    title = 'Low iNKT10'
)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/neighborhood_low.pdf',bbox_inches="tight")


# In[18]:


sq.pl.nhood_enrichment(
    adata_high,
    cluster_key="new_leiden",
    method="average",
    cmap='OrRd',
    figsize=(5, 5),
     title = 'High iNKT10'
)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/neighborhood_high.pdf',bbox_inches="tight")


# In[19]:


# for each fov, for each possible comparison between cell-types report the minkovski distance 
"""
distance_per_cell_type_low_all = {}

for fov in pd.unique(adata_low.obs.fov):
    adata_fov = adata_low[adata_low.obs.fov == fov]
    pos_fov = adata_fov.obs[['CenterX_global_px','CenterY_global_px']]
    dm_fov = distance_matrix(pos_fov, pos_fov)
    dm_fov = tril(dm_fov).toarray()
    distances_fov = pd.DataFrame(dm_fov, index = adata_fov.obs.new_leiden, columns = adata_fov.obs.new_leiden)

    for r in range(1,adata_fov.shape[0]):
        if r == 1:
            c= 0
            combination1 = distances_fov.columns.values[c] + '-' + distances_fov.index.values[r]
            combination2 = distances_fov.index.values[r] + '-' + distances_fov.columns.values[c]
                
            if combination1 in distance_per_cell_type_low_all.keys() :
                distance_per_cell_type_low_all[combination1].append(distances_fov.iat[r,c])
            elif combination2 in distance_per_cell_type_low_all.keys():
                distance_per_cell_type_low_all[combination2].append(distances_fov.iat[r,c])
            else:
                distance_per_cell_type_low_all[combination1] = [distances_fov.iat[r,c]]
        
        else:
            for c in range(r):
                combination1 = distances_fov.columns.values[c] + '-' + distances_fov.index.values[r]
                combination2 = distances_fov.index.values[r] + '-' + distances_fov.columns.values[c]
                
                if combination1 in distance_per_cell_type_low_all.keys() :
                    distance_per_cell_type_low_all[combination1].append(distances_fov.iat[r,c])
                elif combination2 in distance_per_cell_type_low_all.keys():
                    distance_per_cell_type_low_all[combination2].append(distances_fov.iat[r,c])
                else:
                    distance_per_cell_type_low_all[combination1] = [distances_fov.iat[r,c]]

"""


# In[20]:


# for each fov, for each possible comparison between cell-types report the minkovski distance 
"""
distance_per_cell_type_high_all = {}

for fov in pd.unique(adata_high.obs.fov):
    adata_fov = adata_high[adata_high.obs.fov == fov]
    pos_fov = adata_fov.obs[['CenterX_global_px','CenterY_global_px']]
    dm_fov = distance_matrix(pos_fov, pos_fov)
    dm_fov = tril(dm_fov).toarray()
    distances_fov = pd.DataFrame(dm_fov, index = adata_fov.obs.new_leiden, columns = adata_fov.obs.new_leiden)

    for r in range(1,adata_fov.shape[0]):
        if r == 1:
            c= 0
            combination1 = distances_fov.columns.values[c] + '-' + distances_fov.index[r]
            combination2 = distances_fov.index[r] + '-' + distances_fov.columns.values[c]
                
            if combination1 in distance_per_cell_type_high_all.keys() :
                distance_per_cell_type_high_all[combination1].append(distances_fov.iat[r,c])
            elif combination2 in distance_per_cell_type_high_all.keys():
                distance_per_cell_type_high_all[combination2].append(distances_fov.iat[r,c])
            else:
                distance_per_cell_type_high_all[combination1] = [distances_fov.iat[r,c]]
        
        else:
            for c in range(r):
                combination1 = distances_fov.columns.values[c] + '-' + distances_fov.index[r]
                combination2 = distances_fov.index[r] + '-' + distances_fov.columns.values[c]
                
                if combination1 in distance_per_cell_type_high_all.keys() :
                    distance_per_cell_type_high_all[combination1].append(distances_fov.iat[r,c])
                elif combination2 in distance_per_cell_type_high_all.keys():
                    distance_per_cell_type_high_all[combination2].append(distances_fov.iat[r,c])
                else:
                    distance_per_cell_type_high_all[combination1] = [distances_fov.iat[r,c]]

"""


# In[131]:


# compute the t-test between the distribution of distances between two cell types in high and in low.
# bonferroni correction
"""
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests

t_test_results = pd.DataFrame(columns=['comparison', 'statistic', 'pvalue'])
# loop over column_list and execute code explained above
for combination in distance_per_cell_type_high_all.keys():
    combination2 = combination.replace('-', ' ').split(' ')[1] + '-' + combination.replace('-', ' ').split(' ')[0]
    
    if combination in distance_per_cell_type_low_all.keys():
        low = distance_per_cell_type_low_all[combination]
        high = distance_per_cell_type_high_all[combination]
        new_row = [combination, stats.ttest_ind(low,high).statistic, stats.ttest_ind(low,high).pvalue]
        t_test_results.loc[len(t_test_results)]  = new_row
    elif combination2 in distance_per_cell_type_low_all.keys():
        low = distance_per_cell_type_low_all[combination2]
        high = distance_per_cell_type_high_all[combination]
        new_row = [combination, stats.ttest_ind(low,high).statistic, stats.ttest_ind(low,high).pvalue]
        t_test_results.loc[len(t_test_results)]  = new_row


rejected, p_adjusted, _, alpha_corrected = multipletests(t_test_results['pvalue'].values, alpha=0.05, 
                               method='bonferroni', is_sorted=False, returnsorted=False)

t_test_results['p_adjusted'] = p_adjusted

t_test_results.to_csv('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/comparisons_distance_high_low.txt', index =False)
t_test_results
"""


# In[132]:


"""np.array(t_test_results['comparison'])"""


# In[153]:


"""'Epithelial-Mast' in distance_per_cell_type_low_all.keys()"""


# In[154]:


"""
epi_epi = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_all['Epithelial-Epithelial'])),
                          np.repeat('High',len(distance_per_cell_type_high_all['Epithelial-Epithelial']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_all['Epithelial-Epithelial'], distance_per_cell_type_high_all['Epithelial-Epithelial']])})
epi_epi.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_epi/Epithelial.xlsx',
           index=None)

Epithelial_Macrophage = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_all['Macrophage-Epithelial'])),
                          np.repeat('High',len(distance_per_cell_type_high_all['Epithelial-Macrophage']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_all['Macrophage-Epithelial'], distance_per_cell_type_high_all['Epithelial-Macrophage']])})
Epithelial_Macrophage.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_epi/Macrophage.xlsx',
           index=None)

Epithelial_Mast_cell = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_all['Epithelial-Mast'])),
                          np.repeat('High',len(distance_per_cell_type_high_all['Epithelial-Mast']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_all['Epithelial-Mast'], distance_per_cell_type_high_all['Epithelial-Mast']])})
Epithelial_Mast_cell.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_epi/Mast cell.xlsx',
           index=None)

Epithelial_fibro = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_all['Epithelial-Fibroblast'])),
                          np.repeat('High',len(distance_per_cell_type_high_all['Epithelial-Fibroblast']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_all['Epithelial-Fibroblast'], distance_per_cell_type_high_all['Epithelial-Fibroblast']])})
Epithelial_fibro.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_epi/Fibroblast.xlsx',
           index=None)

Epithelial_plasma_iga = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_all['Epithelial-Plasma IgA'])),
                          np.repeat('High',len(distance_per_cell_type_high_all['Epithelial-Plasma IgA']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_all['Epithelial-Plasma IgA'], distance_per_cell_type_high_all['Epithelial-Plasma IgA']])})
Epithelial_plasma_iga.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_epi/Plasma_IgA.xlsx',
           index=None)

Epithelial_plasma_igg = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_all['Epithelial-Plasma IgG'])),
                          np.repeat('High',len(distance_per_cell_type_high_all['Epithelial-Plasma IgG']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_all['Epithelial-Plasma IgG'], distance_per_cell_type_high_all['Epithelial-Plasma IgG']])})
Epithelial_plasma_igg.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_epi/Plasma_IgG.xlsx',
           index=None)

Epithelial_T = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_all['T-Epithelial'])),
                          np.repeat('High',len(distance_per_cell_type_high_all['Epithelial-T']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_all['T-Epithelial'], distance_per_cell_type_high_all['Epithelial-T']])})
Epithelial_T.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_epi/T.xlsx',
           index=None)

Epithelial_Endothelial = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_all['Endothelial-Epithelial'])),
                          np.repeat('High',len(distance_per_cell_type_high_all['Epithelial-Endothelial']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_all['Endothelial-Epithelial'], distance_per_cell_type_high_all['Epithelial-Endothelial']])})
Epithelial_Endothelial.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_epi/Endothelial.xlsx',
           index=None)

Epithelial_Bmem = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_all['Epithelial-B reg'])),
                          np.repeat('High',len(distance_per_cell_type_high_all['Epithelial-B reg']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_all['Epithelial-B reg'], distance_per_cell_type_high_all['Epithelial-B reg']])})
Epithelial_Bmem.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_epi/B memory.xlsx',
           index=None)

Epithelial_B = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_all['Epithelial-B'])),
                          np.repeat('High',len(distance_per_cell_type_high_all['Epithelial-B']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_all['Epithelial-B'], distance_per_cell_type_high_all['Epithelial-B']])})
Epithelial_B.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_epi/B.xlsx',
           index=None)

Epithelial_smooth = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_all['Epithelial-Smooth muscle cell'])),
                          np.repeat('High',len(distance_per_cell_type_high_all['Epithelial-Smooth muscle cell']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_all['Epithelial-Smooth muscle cell'], distance_per_cell_type_high_all['Epithelial-Smooth muscle cell']])})
Epithelial_smooth.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_epi/Smooth muscle cell.xlsx',
           index=None)

Epithelial_glia = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_all['Epithelial-Enteric glia cell'])),
                          np.repeat('High',len(distance_per_cell_type_high_all['Epithelial-Enteric glia cell']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_all['Epithelial-Enteric glia cell'], distance_per_cell_type_high_all['Epithelial-Enteric glia cell']])})
Epithelial_glia.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_epi/Enteric glia cell.xlsx',
           index=None)

Epithelial_undif = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_all['Epithelial-Undifferentiated'])),
                          np.repeat('High',len(distance_per_cell_type_high_all['Epithelial-Undifferentiated']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_all['Epithelial-Undifferentiated'], distance_per_cell_type_high_all['Epithelial-Undifferentiated']])})
Epithelial_undif.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_CD/distances_epi/undifferentiated.xlsx',
           index=None)
"""


# In[133]:


# compute, for each pair of cell-types, the mean distance considering all the distances from all the fovs 
"""
distance_per_cell_type_low_mean = {}

for k in distance_per_cell_type_low_all.keys():
    distance_per_cell_type_low_mean[k] = np.mean(distance_per_cell_type_low_all[k])

distance_per_cell_type_low_mean_df = pd.DataFrame(columns = ['cell_type1', 'cell_type2', 'mean'])

for i in distance_per_cell_type_low_mean.items():
    new_row = [i[0].split('-')[0], i[0].split('-')[1], i[1]]
    distance_per_cell_type_low_mean_df.loc[len(distance_per_cell_type_low_mean_df)]  = new_row

distance_per_cell_type_low_mean_df
distance_per_cell_type_low_mean_df.to_csv('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_low.txt', index =False)
"""


# In[134]:


# compute, for each pair of cell-types, the mean distance considering all the distances from all the fovs 
"""
distance_per_cell_type_high_mean = {}

for k in distance_per_cell_type_high_all.keys():
    distance_per_cell_type_high_mean[k] = np.mean(distance_per_cell_type_high_all[k])

distance_per_cell_type_high_mean_df = pd.DataFrame(columns = ['cell_type1', 'cell_type2', 'mean'])

for i in distance_per_cell_type_high_mean.items():
    new_row = [i[0].split('-')[0], i[0].split('-')[1], i[1]]
    distance_per_cell_type_high_mean_df.loc[len(distance_per_cell_type_high_mean_df)]  = new_row

distance_per_cell_type_high_mean_df
distance_per_cell_type_high_mean_df.to_csv('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_high.txt', index =False)
"""


# In[21]:


# for each fov, compute the mean distance between each pair of cell-types in that fov

distance_per_cell_type_low_fov_mean = {}

for fov in pd.unique(adata_low.obs.fov):
    adata_fov = adata_low[adata_low.obs.fov == fov]
    pos_fov = adata_fov.obs[['CenterX_global_px','CenterY_global_px']]
    dm_fov = distance_matrix(pos_fov, pos_fov)
    dm_fov = tril(dm_fov).toarray()
    distances_fov = pd.DataFrame(dm_fov, index = adata_fov.obs.new_leiden, columns = adata_fov.obs.new_leiden)
    distance_per_cell_type_low = {}
    
    for r in range(1,adata_fov.shape[0]):
        if r == 1:
            c = 0
            combination1 = distances_fov.columns.values[c] + '-' + distances_fov.index[r]
            combination2 = distances_fov.index[r] + '-' + distances_fov.columns.values[c]
                
            if combination1 in distance_per_cell_type_low.keys() :
                distance_per_cell_type_low[combination1].append(distances_fov.iat[r,c])
            elif combination2 in distance_per_cell_type_low.keys():
                distance_per_cell_type_low[combination2].append(distances_fov.iat[r,c])
            else:
                distance_per_cell_type_low[combination1] = [distances_fov.iat[r,c]]
        
        else:
            for c in range(r):
                combination1 = distances_fov.columns.values[c] + '-' + distances_fov.index[r]
                combination2 = distances_fov.index[r] + '-' + distances_fov.columns.values[c]
                
                if combination1 in distance_per_cell_type_low.keys() :
                    distance_per_cell_type_low[combination1].append(distances_fov.iat[r,c])
                elif combination2 in distance_per_cell_type_low.keys():
                    distance_per_cell_type_low[combination2].append(distances_fov.iat[r,c])
                else:
                    distance_per_cell_type_low[combination1] = [distances_fov.iat[r,c]]

    for combination in distance_per_cell_type_low.keys():
        combination2 = combination.split('-')[1] +'-'+combination.split('-')[0]
        if combination in distance_per_cell_type_low_fov_mean.keys():
            distance_per_cell_type_low_fov_mean[combination].append(np.mean(distance_per_cell_type_low[combination]))
        elif combination2 in distance_per_cell_type_low_fov_mean.keys():
            distance_per_cell_type_low_fov_mean[combination2].append(np.mean(distance_per_cell_type_low[combination]))
        else:
            distance_per_cell_type_low_fov_mean[combination] = [np.mean(distance_per_cell_type_low[combination])]


# In[22]:


# for each fov, compute the mean distance between each pair of cell-types in that fov

distance_per_cell_type_high_fov_mean = {}

for fov in pd.unique(adata_high.obs.fov):
    adata_fov = adata_high[adata_high.obs.fov == fov]
    pos_fov = adata_fov.obs[['CenterX_global_px','CenterY_global_px']]
    dm_fov = distance_matrix(pos_fov, pos_fov)
    dm_fov = tril(dm_fov).toarray()
    distances_fov = pd.DataFrame(dm_fov, index = adata_fov.obs.new_leiden, columns = adata_fov.obs.new_leiden)
    distance_per_cell_type_high = {}
    
    for r in range(1,adata_fov.shape[0]):
        if r == 1:
            c= 0
            combination1 = distances_fov.columns.values[c] + '-' + distances_fov.index[r]
            combination2 = distances_fov.index[r] + '-' + distances_fov.columns.values[c]
                
            if combination1 in distance_per_cell_type_high.keys() :
                distance_per_cell_type_high[combination1].append(distances_fov.iat[r,c])
            elif combination2 in distance_per_cell_type_high.keys():
                distance_per_cell_type_high[combination2].append(distances_fov.iat[r,c])
            else:
                distance_per_cell_type_high[combination1] = [distances_fov.iat[r,c]]
        
        else:
            for c in range(r):
                combination1 = distances_fov.columns.values[c] + '-' + distances_fov.index[r]
                combination2 = distances_fov.index[r] + '-' + distances_fov.columns.values[c]
                
                if combination1 in distance_per_cell_type_high.keys() :
                    distance_per_cell_type_high[combination1].append(distances_fov.iat[r,c])
                elif combination2 in distance_per_cell_type_high.keys():
                    distance_per_cell_type_high[combination2].append(distances_fov.iat[r,c])
                else:
                    distance_per_cell_type_high[combination1] = [distances_fov.iat[r,c]]

    for combination in distance_per_cell_type_high.keys():
        combination2 = combination.split('-')[1] +'-'+combination.split('-')[0]
        if combination in distance_per_cell_type_high_fov_mean.keys():
            distance_per_cell_type_high_fov_mean[combination].append(np.mean(distance_per_cell_type_high[combination]))
        elif combination2 in distance_per_cell_type_high_fov_mean.keys():
            distance_per_cell_type_high_fov_mean[combination2].append(np.mean(distance_per_cell_type_high[combination]))
        else:
            distance_per_cell_type_high_fov_mean[combination] = [np.mean(distance_per_cell_type_high[combination])]


# In[23]:


# compute, for each pair of cell-types, the mean distance between the mean distances of the different fovs

distance_per_cell_type_low_fov_mean2 = {}

for k in distance_per_cell_type_low_fov_mean.keys():
    distance_per_cell_type_low_fov_mean2[k] = np.mean(distance_per_cell_type_low_fov_mean[k])

distance_per_cell_type_low_fov_mean2_df = pd.DataFrame(columns = ['cell_type1', 'cell_type2', 'mean'])

for i in distance_per_cell_type_low_fov_mean2.items():
    new_row = [i[0].split('-')[0], i[0].split('-')[1], i[1]]
    distance_per_cell_type_low_fov_mean2_df.loc[len(distance_per_cell_type_low_fov_mean2_df)]  = new_row

distance_per_cell_type_low_fov_mean2_df.to_csv('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_low_mean.txt', index =False)

distance_per_cell_type_low_fov_mean2_df


# In[24]:


# compute, for each pair of cell-types, the mean distance between the mean distances of the different fovs

distance_per_cell_type_high_fov_mean2 = {}

for k in distance_per_cell_type_high_fov_mean.keys():
    distance_per_cell_type_high_fov_mean2[k] = np.mean(distance_per_cell_type_high_fov_mean[k])

distance_per_cell_type_high_fov_mean2_df = pd.DataFrame(columns = ['cell_type1', 'cell_type2', 'mean'])

for i in distance_per_cell_type_high_fov_mean2.items():
    new_row = [i[0].split('-')[0], i[0].split('-')[1], i[1]]
    distance_per_cell_type_high_fov_mean2_df.loc[len(distance_per_cell_type_high_fov_mean2_df)]  = new_row

distance_per_cell_type_high_fov_mean2_df.to_csv('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_high_mean.txt', index =False)

distance_per_cell_type_high_fov_mean2_df


# In[25]:


# for each pair of cell-types, compute the wilcoxon test between the distribution of mean distances in the fovs of high and low conditions
# fdr correction

from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests

t_test_results_mean = pd.DataFrame(columns=['comparison', 'statistic', 'pvalue'])
# loop over column_list and execute code explained above
for combination in distance_per_cell_type_high_fov_mean.keys():
    combination2 = combination.replace('-', ' ').split(' ')[1] + '-' + combination.replace('-', ' ').split(' ')[0]
    
    if combination in distance_per_cell_type_low_fov_mean.keys():
        low = distance_per_cell_type_low_fov_mean[combination]
        high = distance_per_cell_type_high_fov_mean[combination]
        new_row = [combination, stats.ttest_ind(low,high).statistic, stats.ttest_ind(low,high).pvalue]
        t_test_results_mean.loc[len(t_test_results_mean)]  = new_row
    elif combination2 in distance_per_cell_type_low_fov_mean.keys():
        low = distance_per_cell_type_low_fov_mean[combination2]
        high = distance_per_cell_type_high_fov_mean[combination]
        new_row = [combination, stats.ttest_ind(low,high).statistic, stats.ttest_ind(low,high).pvalue]
        t_test_results_mean.loc[len(t_test_results_mean)]  = new_row

rejected, p_adjusted, _, alpha_corrected = multipletests(t_test_results_mean['pvalue'].values, alpha=0.05, 
                               method='fdr_bh', is_sorted=False, returnsorted=False)

t_test_results_mean['p_adjusted'] = p_adjusted

t_test_results_mean.to_csv('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/comparisons_distance_high_low_mean.txt', index =False)
t_test_results_mean


# In[26]:


np.array(t_test_results_mean['comparison'])


# In[63]:


'Epithelial-Enteric glia cells' in distance_per_cell_type_high_fov_mean.keys()


# In[64]:


epi_epi = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_fov_mean['Epithelial-Epithelial'])),
                          np.repeat('High',len(distance_per_cell_type_high_fov_mean['Epithelial-Epithelial']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_fov_mean['Epithelial-Epithelial'], distance_per_cell_type_high_fov_mean['Epithelial-Epithelial']])})
epi_epi.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_epi_mean/Epithelial.xlsx',
           index=None)

Epithelial_Macrophage = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_fov_mean['Macrophage-Epithelial'])),
                          np.repeat('High',len(distance_per_cell_type_high_fov_mean['Epithelial-Macrophage']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_fov_mean['Macrophage-Epithelial'], distance_per_cell_type_high_fov_mean['Epithelial-Macrophage']])})
Epithelial_Macrophage.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_epi_mean/Macrophage.xlsx',
           index=None)

Epithelial_Mast_cell = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_fov_mean['Mast cells-Epithelial'])),
                          np.repeat('High',len(distance_per_cell_type_high_fov_mean['Epithelial-Mast cells']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_fov_mean['Mast cells-Epithelial'], distance_per_cell_type_high_fov_mean['Epithelial-Mast cells']])})
Epithelial_Mast_cell.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_epi_mean/Mast cell.xlsx',
           index=None)

Epithelial_fibro = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_fov_mean['Fibroblast-Epithelial'])),
                          np.repeat('High',len(distance_per_cell_type_high_fov_mean['Epithelial-Fibroblast']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_fov_mean['Fibroblast-Epithelial'], distance_per_cell_type_high_fov_mean['Epithelial-Fibroblast']])})
Epithelial_fibro.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_epi_mean/Fibroblast.xlsx',
           index=None)

Epithelial_plasma_iga = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_fov_mean['Plasma IgA-Epithelial'])),
                          np.repeat('High',len(distance_per_cell_type_high_fov_mean['Epithelial-Plasma IgA']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_fov_mean['Plasma IgA-Epithelial'], distance_per_cell_type_high_fov_mean['Epithelial-Plasma IgA']])})
Epithelial_plasma_iga.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_epi_mean/Plasma_IgA.xlsx',
           index=None)

Epithelial_plasma_igg = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_fov_mean['Plasma IgG-Epithelial'])),
                          np.repeat('High',len(distance_per_cell_type_high_fov_mean['Epithelial-Plasma IgG']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_fov_mean['Plasma IgG-Epithelial'], distance_per_cell_type_high_fov_mean['Epithelial-Plasma IgG']])})
Epithelial_plasma_igg.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_epi_mean/Plasma_IgG.xlsx',
           index=None)

Epithelial_T = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_fov_mean['T-Epithelial'])),
                          np.repeat('High',len(distance_per_cell_type_high_fov_mean['Epithelial-T']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_fov_mean['T-Epithelial'], distance_per_cell_type_high_fov_mean['Epithelial-T']])})
Epithelial_T.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_epi_mean/T.xlsx',
           index=None)

Epithelial_Endothelial = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_fov_mean['Endothelial-Epithelial'])),
                          np.repeat('High',len(distance_per_cell_type_high_fov_mean['Epithelial-Endothelial']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_fov_mean['Endothelial-Epithelial'], distance_per_cell_type_high_fov_mean['Epithelial-Endothelial']])})
Epithelial_Endothelial.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_epi_mean/Endothelial.xlsx',
           index=None)

Epithelial_Bmem = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_fov_mean['B reg-Epithelial'])),
                          np.repeat('High',len(distance_per_cell_type_high_fov_mean['Epithelial-B reg']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_fov_mean['B reg-Epithelial'], distance_per_cell_type_high_fov_mean['Epithelial-B reg']])})
Epithelial_Bmem.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_epi_mean/B memory.xlsx',
           index=None)

Epithelial_B = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_fov_mean['B-Epithelial'])),
                          np.repeat('High',len(distance_per_cell_type_high_fov_mean['Epithelial-B']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_fov_mean['B-Epithelial'], distance_per_cell_type_high_fov_mean['Epithelial-B']])})
Epithelial_B.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_epi_mean/B.xlsx',
           index=None)

Epithelial_smooth = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_fov_mean['Smooth muscle cells-Epithelial'])),
                          np.repeat('High',len(distance_per_cell_type_high_fov_mean['Epithelial-Smooth muscle cells']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_fov_mean['Smooth muscle cells-Epithelial'], distance_per_cell_type_high_fov_mean['Epithelial-Smooth muscle cells']])})
Epithelial_smooth.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_epi_mean/Smooth muscle cell.xlsx',
           index=None)

Epithelial_glia = pd.DataFrame(data = {'iNKT10 level' : np.concatenate([np.repeat('Low',len(distance_per_cell_type_low_fov_mean['Enteric glia cells-Epithelial'])),
                          np.repeat('High',len(distance_per_cell_type_high_fov_mean['Epithelial-Enteric glia cells']))]), 
                          'Distance' : np.concatenate([distance_per_cell_type_low_fov_mean['Enteric glia cells-Epithelial'], distance_per_cell_type_high_fov_mean['Epithelial-Enteric glia cells']])})
Epithelial_glia.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/distances_epi_mean/Enteric glia cell.xlsx',
           index=None)


# ### Compute Ripley’s statistics

# In[147]:


mode = "L"
sq.gr.ripley(adata_low, cluster_key="new_leiden", mode=mode)
sq.pl.ripley(
    adata_low,
    cluster_key="new_leiden",
    mode=mode
)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/ripley_low.pdf',bbox_inches="tight")


# In[158]:


mode = "L"
sq.gr.ripley(adata_high, cluster_key="new_leiden", mode=mode)
sq.pl.ripley(
    adata_high,
    cluster_key="new_leiden",
    mode=mode
)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/ripley_high.pdf',bbox_inches="tight")


# In[159]:


ripley_per_fov_low = pd.DataFrame(index = pd.unique(adata_low.obs.fov), columns = pd.unique(adata_low.obs.new_leiden))


# In[160]:


mode = "L"

for fov in pd.unique(adata_low.obs.fov):
    adata_fov = adata_low[adata_low.obs.fov == fov]
    sq.gr.ripley(adata_fov, cluster_key="new_leiden", mode=mode, max_dist = 500)

    a = adata_fov.uns['new_leiden_ripley_L']['L_stat']
    pvalues = []

    for inner_list in adata_fov.uns['new_leiden_ripley_L']['pvalues']:
        pvalues.extend(inner_list)

    a['p_value'] = pvalues

    sub_a = a[a.bins == 500]
    for c in pd.unique(adata_low.obs['new_leiden']):
        ripley_per_fov_low.loc[fov,c] = sub_a[sub_a.new_leiden==c]['stats'].values


# In[164]:


import openpyxl
ripley_per_fov_low.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/ripley_per_fov_low.xlsx')


# In[165]:


ripley_per_fov_high = pd.DataFrame(index = pd.unique(adata_high.obs.fov), columns = pd.unique(adata_high.obs.new_leiden))


# In[166]:


mode = "L"

for fov in pd.unique(adata_high.obs.fov):
    adata_fov = adata_high[adata_high.obs.fov == fov]
    sq.gr.ripley(adata_fov, cluster_key="new_leiden", mode=mode, max_dist = 500)

    a = adata_fov.uns['new_leiden_ripley_L']['L_stat']
    pvalues = []

    for inner_list in adata_fov.uns['new_leiden_ripley_L']['pvalues']:
        pvalues.extend(inner_list)

    a['p_value'] = pvalues

    sub_a = a[a.bins == 500]
    for c in pd.unique(adata_high.obs['new_leiden']):
        ripley_per_fov_high.loc[fov,c] = sub_a[sub_a.new_leiden==c]['stats'].values


# In[175]:


ripley_per_fov_high.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/ripley_per_fov_high.xlsx')


# ### Compute centrality scores

# In[168]:


degree_per_fov_high = pd.DataFrame(index = pd.unique(adata_high.obs.fov), columns = pd.unique(adata_high.obs.new_leiden))
average_per_fov_high = pd.DataFrame(index = pd.unique(adata_high.obs.fov), columns = pd.unique(adata_high.obs.new_leiden))
closeness_per_fov_high = pd.DataFrame(index = pd.unique(adata_high.obs.fov), columns = pd.unique(adata_high.obs.new_leiden))


# In[169]:


for fov in pd.unique(adata_high.obs.fov):
    adata_fov = adata_high[adata_high.obs.fov == fov]
    adata_spatial_neighbor = sq.gr.spatial_neighbors(
        adata_fov, coord_type="generic", delaunay=True)

    sq.gr.centrality_scores(adata_fov, cluster_key="new_leiden")

    a = adata_fov.uns['new_leiden_centrality_scores']
    
    for c in a.index.values:
        degree_per_fov_high.loc[fov,c] = a.loc[c]['degree_centrality']

    for c in a.index.values:
        average_per_fov_high.loc[fov,c] = a.loc[c]['average_clustering']

    for c in a.index.values:
        closeness_per_fov_high.loc[fov,c] = a.loc[c]['closeness_centrality']
        


# In[171]:


degree_per_fov_low = pd.DataFrame(index = pd.unique(adata_low.obs.fov), columns = pd.unique(adata_low.obs.new_leiden))
average_per_fov_low = pd.DataFrame(index = pd.unique(adata_low.obs.fov), columns = pd.unique(adata_low.obs.new_leiden))
closeness_per_fov_low = pd.DataFrame(index = pd.unique(adata_low.obs.fov), columns = pd.unique(adata_low.obs.new_leiden))


# In[172]:


for fov in pd.unique(adata_low.obs.fov):
    adata_fov = adata_low[adata_low.obs.fov == fov]
    adata_spatial_neighbor = sq.gr.spatial_neighbors(
        adata_fov, coord_type="generic", delaunay=True)

    sq.gr.centrality_scores(adata_fov, cluster_key="new_leiden")

    a = adata_fov.uns['new_leiden_centrality_scores']
    
    for c in a.index.values:
        degree_per_fov_low.loc[fov,c] = a.loc[c]['degree_centrality']

    for c in a.index.values:
        average_per_fov_low.loc[fov,c] = a.loc[c]['average_clustering']

    for c in a.index.values:
        closeness_per_fov_low.loc[fov,c] = a.loc[c]['closeness_centrality']


# In[170]:


degree_per_fov_high.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/degree_per_fov_high.xlsx')
average_per_fov_high.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/average_per_fov_high.xlsx')
closeness_per_fov_high.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/closeness_per_fov_high.xlsx')


# In[173]:


degree_per_fov_low.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/degree_per_fov_low.xlsx')
average_per_fov_low.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/average_per_fov_low.xlsx')
closeness_per_fov_low.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/closeness_per_fov_low.xlsx')


# ### Compute co-occurrence probability

# In[174]:


co_occurrence_per_fov_low = pd.DataFrame(index = pd.unique(adata_low.obs.fov), columns = pd.unique(adata_low.obs.new_leiden))
co_occurrence_per_fov_high = pd.DataFrame(index = pd.unique(adata_high.obs.fov), columns = pd.unique(adata_high.obs.new_leiden))


# In[78]:


adata_subset = adata[adata.obs.fov == "28"].copy()


# In[80]:


sq.gr.co_occurrence(
    adata_subset,
    cluster_key="new_leiden"
)

sq.pl.co_occurrence(
    adata_subset,
    cluster_key="new_leiden"
)
plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/co_occ_ratio_fov28_high.pdf',bbox_inches="tight")


# In[81]:


adata_subset = adata[adata.obs.fov == "38"].copy()


# In[82]:


sq.gr.co_occurrence(
    adata_subset,
    cluster_key="new_leiden"
)

sq.pl.co_occurrence(
    adata_subset,
    cluster_key="new_leiden"
)
plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/co_occ_ratio_fov38_low.pdf',bbox_inches="tight")


# ### Neighboorhod enrichment analysis 

# In[52]:


for fov in pd.unique(adata.obs.fov):
    adata_fov = adata[adata.obs.fov == fov]
    sq.gr.spatial_neighbors(
    adata_fov,
    n_neighs=10,
    coord_type="generic")

    sq.gr.nhood_enrichment(adata_fov, cluster_key="new_leiden", seed=10)

    pd.DataFrame(adata_fov.uns['new_leiden_nhood_enrichment']['zscore'],
             index = np.sort(pd.unique(adata_fov.obs.new_leiden)), columns = np.sort(pd.unique(adata_fov.obs.new_leiden))).to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/neighboorhod/neighboorhod_'+str(fov)+'.xlsx')


# # Novae
# 
# First, we need to load a pretrained Novae model.
# 
# Since we have a human colon slide, we load the "MICS-Lab/novae-human-0" model. Other model names can be found here (spot-resolution models, e.g. for Visium data, will come soon).

# In[24]:


model = novae.Novae.from_pretrained("MICS-Lab/novae-human-0")
model


# Create a Delaunay graph from the spatial coordinates of the cells (in microns). The graph is stored in adata.obsp['spatial_connectivities'] and adata.obsp['spatial_distances']. The long edges are removed from the graph according to the radius argument (if provided).
# The spatial coordinates are expected to be in microns, and stored in adata.obsm["spatial"]. If the coordinates are in pixels, set pixel_size to the size of a pixel in microns. If you don't know the pixel_size, or if you don't have adata.obsm["spatial"], you can also provide the technology argument.

# In[25]:


# this step because the spatial coordinates are not in micron but in pixel
novae.utils.spatial_neighbors(adata, slide_key='fov',technology='cosmx')


# Now, we can compute the spatial representations for each cell. In the first option below, we pass the argument zero_shot=True to run only inference (i.e., the model is not re-trained).
# 
# Instead of zero-shot, if you want to fine-tune the model, you can use the fine_tune method, and then call compute_representations without the zero_shot argument (see "option 2").

# In[26]:


# Option 1: zero-shot
model.compute_representations(adata, zero_shot=True)

# Option 2: fine-tuning
#model.fine_tune(adata)
#model.compute_representations(adata)


# To assign domains, you can use the assign_domains method, as below. By default, it creates 7 domains, but you can choose the number of domains you want with the level argument.
# 
# The function will save the domains in adata.obs, and return the name of the column in which it was saved 

# In[59]:


model.assign_domains(adata, level=10)


# In[27]:


model.assign_domains(adata, level=20)


# Then, to show the results, you can use novae.plot.domains.
# If you run model.assign_domains multiple times, you can also decide the resolution you want to show, by passing the obs_key argument to novae.plot.domains.

# In[55]:


novae.plot.domains(adata, library_id=None)


# In[94]:


palette = ['#0F2080', '#006600', '#AA3377', '#BBBBBB', '#AEC7E8', '#CCBB44', '#F7B6D2','#000000',
          '#0077BB', '#EE7733', '#33BBEE', '#CC3311', '#7F7F7F', '#FFB000','#601A4A','#785EF0','#66CC00', '#CCFF99', 
          '#660000', '#E5CCFF','#FFFFCC','#CCFFE5']

sns.color_palette(palette)


# In[95]:


sc.pl.umap(
    adata, show=False,
    color= "novae_domains_20",
    palette=sns.color_palette(palette)
)


# In[60]:


t = ['fov 21 - Patient_11 - LOW',
 'fov 23 - Patient_12 - HIGH',
 'fov 24 - Patient_12 - HIGH',
 'fov 25 - Patient_13 - LOW',
 'fov 27 - Patient_14 - HIGH',
 'fov 28 - Patient_14 - HIGH',
 'fov 29 - Patient_15 - LOW',
 'fov 30 - Patient_15 - LOW',
 'fov 31 - Patient_16 - HIGH',
 'fov 32 - Patient_16 - HIGH',
 'fov 33 - Patient_17 - HIGH',
 'fov 34 - Patient_17 - HIGH',
 'fov 35 - Patient_18 - HIGH',
 'fov 36 - Patient_18 - HIGH',
 'fov 37 - Patient_19 - LOW', 
 'fov 38 - Patient_19 - LOW',
 'fov 39 - Patient_20 - LOW',
 'fov 40 - Patient_20 - LOW',
 'fov 41 - Patient_20 - LOW']

sq.pl.spatial_segment(
    adata=adata,
    color='novae_domains_10',
    library_key="fov",
    seg_cell_id="cell_ID",
    seg_outline=True,
    img=False,
    scalebar_dx=1.0,
    title= t,
    scalebar_kwargs={"scale_loc": "bottom", "location": "lower right"},
)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/niche_all_fov_novae_10domains.pdf',bbox_inches="tight")


# In[96]:


t = ['fov 21 - Patient_11 - LOW',
 'fov 23 - Patient_12 - HIGH',
 'fov 24 - Patient_12 - HIGH',
 'fov 25 - Patient_13 - LOW',
 'fov 27 - Patient_14 - HIGH',
 'fov 28 - Patient_14 - HIGH',
 'fov 29 - Patient_15 - LOW',
 'fov 30 - Patient_15 - LOW',
 'fov 31 - Patient_16 - HIGH',
 'fov 32 - Patient_16 - HIGH',
 'fov 33 - Patient_17 - HIGH',
 'fov 34 - Patient_17 - HIGH',
 'fov 35 - Patient_18 - HIGH',
 'fov 36 - Patient_18 - HIGH',
 'fov 37 - Patient_19 - LOW', 
 'fov 38 - Patient_19 - LOW',
 'fov 39 - Patient_20 - LOW',
 'fov 40 - Patient_20 - LOW',
 'fov 41 - Patient_20 - LOW']

sq.pl.spatial_segment(
    adata=adata,
    color='novae_domains_20',
    library_key="fov",
    seg_cell_id="cell_ID",
    seg_outline=True,
    img=False,
    scalebar_dx=1.0,
    title= t,
    scalebar_kwargs={"scale_loc": "bottom", "location": "lower right"},
)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/niche_all_fov_novae_20domains.pdf',bbox_inches="tight")


# Novae has a hierarchical organization of the spatial domains. That is, if you run multiple times assign_domains with different level parameters, the domains at different resolutions will be nested inside each other.
# 
# To plot the hierarchy of the domains, you can use the plot_domains_hierarchy method of Novae, as below:

# In[85]:


model.plot_domains_hierarchy(max_level=30)


# The spatial representations of each cell are stored in adata.obsm["novae_latent"], and they are not batch-effect corrected by default. Yet, the (categorical) spatial domains are corrected. Therefore, we can use the categorical spatial domains to correct the representations, using the batch_effect_correction method, as below:

# In[63]:


adata.obs.to_csv('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/metadata_adata_niche_novae_10domains.csv')


# In[48]:


adata.obs.to_csv('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/metadata_adata_niche_novae_20domains.csv')


# ### Domains proportion per slide

# In[47]:


novae.plot.domains_proportions(adata,obs_key='novae_domains_10')


# ### Slide architecture
# 
# We run trajectory inference (PAGA) on the spatial domains to extract a graph representing the "architecture" of a slide, or the "spatial domains organization".

# In[22]:


novae.plot.paga(adata)


# ### Spatially Variable Genes (SVG)
# To extract SVG, we run DEGs on the categorical spatial domains. The function below shows the 3 most variable genes.

# In[52]:


adata_niche_D997 = adata[adata.obs.novae_domains_20 =='D997']
adata_niche_D980 = adata[adata.obs.novae_domains_20 =='D980']
adata_niche_D1002 = adata[adata.obs.novae_domains_20 =='D1002']
adata_niche_D993 = adata[adata.obs.novae_domains_20 =='D993']


# In[53]:


sc.tl.rank_genes_groups(adata_niche_D997, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
adata_niche_D997 = sc.get.rank_genes_groups_df(adata_niche_D997, group='LOW',key = "rank_genes_groups")
adata_niche_D997.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_niche_D997.xlsx')

#sc.tl.rank_genes_groups(adata_niche_D980, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
#adata_niche_D980 = sc.get.rank_genes_groups_df(adata_niche_D980, group='LOW',key = "rank_genes_groups")
#adata_niche_D980.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_niche_D980.xlsx')

sc.tl.rank_genes_groups(adata_niche_D1002, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
adata_niche_D1002 = sc.get.rank_genes_groups_df(adata_niche_D1002, group='LOW',key = "rank_genes_groups")
adata_niche_D1002.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_niche_D1002.xlsx')

sc.tl.rank_genes_groups(adata_niche_D993, 'level', groups=['LOW'], reference='HIGH', method='wilcoxon')
adata_niche_D993 = sc.get.rank_genes_groups_df(adata_niche_D993, group='LOW',key = "rank_genes_groups")
adata_niche_D993.to_excel('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/DEGs_niche_D993.xlsx')


# In[78]:


adata_niche_D980_mean_exp = pd.DataFrame(adata_niche_D980.X.mean(0).transpose(),
    index=adata_niche_D980.var_names, columns = ['mean_exp']
)


# In[79]:


adata_niche_D980_mean_exp


# In[82]:


adata_niche_D980_mean_exp['mean_exp'].mean()


# In[93]:


with open('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/higher_exp_niche_D980.txt', "w") as txt_file:
    for line in adata_niche_D980_mean_exp.sort_values(by='mean_exp', ascending=False)[1:80].index.values:
        txt_file.write(line + "\n") # works with any number of elements in a line


# In[61]:


import inspect
lines = inspect.getsource(novae.plot.pathway_scores)
print(lines)


# In[ ]:


def pathway_scores(
    adata: AnnData,
    pathways: dict[str, list[str]] | str,
    obs_key: str | None = None,
    return_df: bool = False,
    figsize: tuple[int, int] = (10, 5),
    min_pathway_size: int = 4,
    **kwargs: int,
) -> pd.DataFrame | None:
    """Show a heatmap of pathway scores for each domain.

    Info:
        Currently, this function only supports one slide per call.

    Args:
        adata: An `AnnData` object.
        pathways: Either a dictionary of pathways (keys are pathway names, values are lists of gane names), or a path to a [GSEA](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp) JSON file.
        obs_key: Key in `adata.obs` that contains the domains. By default, it will use the last available Novae domain key.
        return_df: Whether to return the DataFrame.
        figsize: Matplotlib figure size.
        min_pathway_size: Minimum number of known genes in the pathway to be considered.

    Returns:
        A DataFrame of scores per domain if `return_df` is True.
    """
    assert isinstance(adata, AnnData), f"For now, only AnnData objects are supported, received {type(adata)}"

    obs_key = utils.check_available_domains_key([adata], obs_key)

    scores = {}
    lower_var_names = adata.var_names.str.lower()

    if isinstance(pathways, str):
        pathways = _load_gsea_json(pathways)
        log.info(f"Loaded {len(pathways)} pathway(s)")

    for key, gene_names in pathways.items():
        vars = np.array([gene_name.lower() for gene_name in gene_names])
        vars = adata.var_names[np.isin(lower_var_names, vars)]
        if len(vars) >= min_pathway_size:
            sc.tl.score_genes(adata, vars, score_name="_temp")
            scores[key] = adata.obs["_temp"]
    del adata.obs["_temp"]

    assert len(scores) > 1, f"Found {len(scores)} valid pathway. Minimum 2 required."

    df = pd.DataFrame(scores)
    df[obs_key] = adata.obs[obs_key]
    df = df.groupby(obs_key).mean()
    df = df.fillna(0)

    g = sns.clustermap(df, figsize=figsize, **kwargs)
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

    if return_df:
        return df


# In[69]:


def spatially_variable_genes(
    adata,
    obs_key = None,
    top_k: int = 5,
    show: bool = False,
    cell_size: int = 10,
    min_positive_ratio: float = 0.05,
    return_list: bool = False,
    **kwargs: int,
):
    """Plot the most spatially variable genes (SVG) for a given `AnnData` object.

    !!! info
        Currently, this function only supports one slide per call.

    Args:
        adata: An `AnnData` object corresponding to one slide.
        obs_key: Key in `adata.obs` that contains the domains. By default, it will use the last available Novae domain key.
        top_k: Number of SVG to be shown.
        show: Whether to show the plot.
        cell_size: Size of the cells or spots (`spot_size` argument of `sc.pl.spatial`).
        min_positive_ratio: Genes whose "ratio of cells expressing it" is lower than this threshold are not considered.
        return_list: Whether to return the list of SVG instead of plotting them.
        **kwargs: Additional arguments for `sc.pl.spatial`.

    Returns:
        A list of SVG names if `return_list` is `True`.
    """
    
    obs_key = check_available_domains_key([adata], obs_key)

    sc.tl.rank_genes_groups(adata, groupby=obs_key)
    df = pd.concat(
        [
            sc.get.rank_genes_groups_df(adata, domain).set_index("names")["logfoldchanges"]
            for domain in adata.obs[obs_key].cat.categories
        ],
        axis=1,
    )

    where = (adata.X > 0).mean(0) > min_positive_ratio
    valid_vars = adata.var_names[where.A1 if isinstance(where, np.matrix) else where]
    assert (
        len(valid_vars) >= top_k
    ), f"Only {len(valid_vars)} genes are available. Please decrease `top_k` or `min_positive_ratio`."

    svg = df.std(1).loc[valid_vars].sort_values(ascending=False).head(top_k).index

    if return_list:
        return svg.tolist()


# In[74]:


SVGs = novae.plot.spatially_variable_genes(adata, top_k=10, library_id=None, return_list=True)


# In[75]:


adata28 = adata[adata.obs.fov =='28']


# In[76]:


sq.pl.spatial_segment(
    adata= adata28,
    color=SVGs,
    library_key="fov",
    seg_cell_id="cell_ID",
    seg_outline=True,
    img=False,
    scalebar_dx=1.0,
    scalebar_kwargs={"scale_loc": "bottom", "location": "lower right"},
)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/SVG_fov28_novae_10domains.pdf',bbox_inches="tight")


# In[85]:


SVGs = novae.plot.spatially_variable_genes(adata, top_k=10, library_id=None, return_list=True, obs_key='novae_domains_20')


# In[86]:


sq.pl.spatial_segment(
    adata= adata28,
    color=SVGs,
    library_key="fov",
    seg_cell_id="cell_ID",
    seg_outline=True,
    img=False,
    scalebar_dx=1.0,
    scalebar_kwargs={"scale_loc": "bottom", "location": "lower right"},
)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/SVG_fov28_novae_20domains.pdf',bbox_inches="tight")


# In[28]:


# Import JSON module
import json

f = open("/Users/paolamaragno/Downloads/h.all.v2024.1.Hs.json") 

data = json.load(f) 


# In[29]:


data = {key: data[key] for key in ['HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_WNT_BETA_CATENIN_SIGNALING','HALLMARK_TGF_BETA_SIGNALING','HALLMARK_IL6_JAK_STAT3_SIGNALING',
'HALLMARK_DNA_REPAIR','HALLMARK_G2M_CHECKPOINT','HALLMARK_APOPTOSIS','HALLMARK_NOTCH_SIGNALING','HALLMARK_MYOGENESIS','HALLMARK_PROTEIN_SECRETION',
'HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_APICAL_JUNCTION','HALLMARK_APICAL_SURFACE','HALLMARK_HEDGEHOG_SIGNALING',
'HALLMARK_COMPLEMENT','HALLMARK_UNFOLDED_PROTEIN_RESPONSE','HALLMARK_PI3K_AKT_MTOR_SIGNALING','HALLMARK_E2F_TARGETS','HALLMARK_MYC_TARGETS_V1',
'HALLMARK_MYC_TARGETS_V2','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION','HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_XENOBIOTIC_METABOLISM',
'HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_GLYCOLYSIS','HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY','HALLMARK_P53_PATHWAY','HALLMARK_ANGIOGENESIS',
'HALLMARK_HEME_METABOLISM','HALLMARK_COAGULATION','HALLMARK_IL2_STAT5_SIGNALING','HALLMARK_PEROXISOME','HALLMARK_KRAS_SIGNALING_UP',
'HALLMARK_KRAS_SIGNALING_DN','HALLMARK_MTORC1_SIGNALING']}


# In[30]:


data


# In[89]:


len(data)


# In[41]:


pathways = {}

for k in data.keys():
    pathways[k] = data[k]['geneSymbols']


# In[81]:


novae.plot.pathway_scores(adata, obs_key ="novae_domains_10",pathways=pathways, figsize=(15, 7))

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/pathways_novae_10domains.pdf',bbox_inches="tight")


# In[45]:


pathway_score = novae.plot.pathway_scores(adata, obs_key ="novae_domains_20",pathways=pathways, figsize=(10, 7),return_df=True)
pathway_score
#plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/pathways_novae_20domains.pdf',bbox_inches="tight")


# In[ ]:


pathway_score.to_csv('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/pathways_novae_20domains.csv') 


# In[92]:


obs_key ="novae_domains_20"

scores = {}
lower_var_names = adata.var_names.str.lower()

for key, gene_names in pathways.items():
    vars = np.array([gene_name.lower() for gene_name in gene_names])
    vars = adata.var_names[np.isin(lower_var_names, vars)]
    if len(vars) >= 4:
        sc.tl.score_genes(adata, vars, score_name="_temp")
        scores[key] = adata.obs["_temp"]
del adata.obs["_temp"]

assert len(scores) > 1, f"Found {len(scores)} valid pathway. Minimum 2 required."

df = pd.DataFrame(scores)
df[obs_key] = adata.obs[obs_key]


# In[94]:


df
df.to_csv('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/pathways_novae_20domains_score_by_cells.csv') 


# # Niche compass

# In[283]:


from nichecompass.utils import (add_gps_from_gp_dict_to_adata,
                                compute_communication_gp_network,
                                visualize_communication_gp_network,
                                create_new_color_dict,
                                extract_gp_dict_from_nichenet_lrt_interactions,
                                extract_gp_dict_from_omnipath_lr_interactions,
                                generate_enriched_gp_info_plots)


# ### Create Prior Knowledge Gene Program (GP) Mask
# 
# NicheCompass expects a prior GP mask as input, which it will use to make its latent feature space interpretable (through linear masked decoders).
# 
# The user can provide a custom GP mask to NicheCompass based on the biological question of interest.
# 
# As a default, we create a GP mask based on four databases of prior knowledge of inter- and intracellular interaction pathways:
# 
# - OmniPath (Ligand-Receptor GPs)
# 
# - MEBOCOST (Enzyme-Sensor GPs)
# 
# - CollecTRI (Transcriptional Regulation GPs)
# 
# - NicheNet (Combined Interaction GPs)
# 

# In[36]:


# Retrieve OmniPath GPs (source: ligand genes; target: receptor genes)

figure_folder_path = '/Users/paolamaragno/Downloads/NicheCompass'

omnipath_gp_dict = extract_gp_dict_from_omnipath_lr_interactions(
    species='human',
    load_from_disk=False,
    save_to_disk=True,
    lr_network_file_path=f"{figure_folder_path}" "/omnipath_lr_network.csv",
    gene_orthologs_mapping_file_path=f"{figure_folder_path}" "/human_mouse_gene_orthologs.csv",
    plot_gp_gene_count_distributions=True,
    gp_gene_count_distributions_save_path=f"{figure_folder_path}" \
                                           "/omnipath_gp_gene_count_distributions.pdf")


# In[23]:


# Display example OmniPath GP
import random

omnipath_gp_names = list(omnipath_gp_dict.keys())
random.shuffle(omnipath_gp_names)
omnipath_gp_name = omnipath_gp_names[0]
print(f"{omnipath_gp_name}: {omnipath_gp_dict[omnipath_gp_name]}")


# In[24]:


# Retrieve NicheNet GPs (source: ligand genes; target: receptor genes, target genes)
nichenet_gp_dict = extract_gp_dict_from_nichenet_lrt_interactions(
    species='human',
    version="v2",
    keep_target_genes_ratio=1.,
    max_n_target_genes_per_gp=250,
    load_from_disk=False,
    save_to_disk=True,
    lr_network_file_path=f"{figure_folder_path}" "/nichenet_lr_network_v2_human.csv",
    ligand_target_matrix_file_path=f"{figure_folder_path}" "/nichenet_ligand_target_matrix_v2_human.csv",
    gene_orthologs_mapping_file_path= f"{figure_folder_path}" "/human_mouse_gene_orthologs.csv",
    plot_gp_gene_count_distributions=True)


# In[22]:


# Display example NicheNet GP
nichenet_gp_names = list(nichenet_gp_dict.keys())
random.shuffle(nichenet_gp_names)
nichenet_gp_name = nichenet_gp_names[0]
print(f"{nichenet_gp_name}: {nichenet_gp_dict[nichenet_gp_name]}")


# In[25]:


from typing import Optional

def extract_gp_dict_from_mebocost_ms_interactions(
        species: ["mouse", "human"],
        dir_path: str="/Users/paolamaragno/Downloads/metabolite_enzyme_sensor_gps",
        plot_gp_gene_count_distributions: bool=True,
        gp_gene_count_distributions_save_path: Optional[str]=None) -> dict:
    if species == "human":
        metabolite_enzymes_df = pd.read_csv(
            dir_path + "/human_metabolite_enzymes.tsv", sep="\t")
        metabolite_sensors_df = pd.read_csv(
            dir_path + "/human_metabolite_sensors.tsv", sep="\t")

    elif species == "mouse":
        metabolite_enzymes_df = pd.read_csv(
            dir_path + "/mouse_metabolite_enzymes.tsv", sep="\t")
        metabolite_sensors_df = pd.read_csv(
            dir_path + "/mouse_metabolite_sensors.tsv", sep="\t")
    else:
        raise ValueError("Species should be either human or mouse.")

    # Retrieve metabolite names
    metabolite_names_df = (metabolite_sensors_df[["HMDB_ID",
                                                  "standard_metName"]]
                           .drop_duplicates()
                           .set_index("HMDB_ID"))

    # Keep only enzymes for which the metabolite is the product (filter enzymes
    # for which the metabolite is the substrate)
    metabolite_enzymes_df = metabolite_enzymes_df[
        metabolite_enzymes_df["direction"] == "product"]

    # Retrieve metabolite enzyme and sensor genes
    metabolite_enzymes_unrolled = []
    for _, row in metabolite_enzymes_df.iterrows():
        genes = row["gene"].split("; ")
        for gene in genes:
            tmp = row.copy()
            tmp["gene"] = gene
            metabolite_enzymes_unrolled.append(tmp)
    metabolite_enzymes_df = pd.DataFrame(metabolite_enzymes_unrolled)
    metabolite_enzymes_df["gene_name"] = metabolite_enzymes_df["gene"].apply(
        lambda x: x.split("[")[0])
    metabolite_enzymes_df = (metabolite_enzymes_df.groupby(["HMDB_ID"])
                             .agg({"gene_name": lambda x: sorted(
                                x.unique().tolist())})
                             .rename({"gene_name": "enzyme_genes"}, axis=1)
                             .reset_index()).set_index("HMDB_ID")
    metabolite_sensors_df = (metabolite_sensors_df.groupby(["HMDB_ID"])
                             .agg({"Gene_name": lambda x: sorted(
                                x.unique().tolist())})
                             .rename({"Gene_name": "sensor_genes"}, axis=1)
                             .reset_index()).set_index("HMDB_ID")

    # Combine enzyme and sensor genes based on metabolite names (sensor genes
    # are not available for most metabolites)
    metabolite_df = metabolite_enzymes_df.join(
        other=metabolite_sensors_df,
        how="inner").join(metabolite_names_df).set_index("standard_metName")

    # Convert to gene program dictionary format
    met_interaction_dict = metabolite_df.to_dict()
    gp_dict = {}
    for metabolite, enzyme_genes in met_interaction_dict["enzyme_genes"].items():
        gp_dict[metabolite + "_metabolite_enzyme_sensor_GP"] = {
            "sources": enzyme_genes,
            "sources_categories": ["enzyme"] * len(enzyme_genes)}
    for metabolite, sensor_genes in met_interaction_dict["sensor_genes"].items():
        gp_dict[metabolite + "_metabolite_enzyme_sensor_GP"][
            "targets"] = sensor_genes
        gp_dict[metabolite + "_metabolite_enzyme_sensor_GP"][
            "targets_categories"] = ["sensor"] * len(sensor_genes)

    if plot_gp_gene_count_distributions:
        create_gp_gene_count_distribution_plots(
            gp_dict=gp_dict,
            gp_plot_label="MEBOCOST",
            save_path=gp_gene_count_distributions_save_path)

    return gp_dict


# In[27]:


from anndata import AnnData
def create_gp_gene_count_distribution_plots(
        gp_dict: Optional[dict]=None,
        adata: Optional[AnnData]=None,
        gp_targets_mask_key: Optional[str]="nichecompass_gp_targets",
        gp_sources_mask_key: Optional[str]="nichecompass_gp_sources",
        gp_plot_label: str="",
        save_path: Optional[str]=None):
    """
    Create distribution plots of the gene counts for sources and targets
    of all gene programs in either a gp dict or an adata object.

    Parameters
    ----------
    gp_dict:
        A gene program dictionary.
    adata:
        An anndata object
    gp_plot_label:
        Label of the gene program plot for title.
    """
    # Get number of source and target genes for each gene program
    if gp_dict is not None:
        n_sources_list = []
        n_targets_list = []
        for _, gp_sources_targets_dict in gp_dict.items():
            n_sources_list.append(len(gp_sources_targets_dict["sources"]))
            n_targets_list.append(len(gp_sources_targets_dict["targets"]))
    elif adata is not None:
        n_targets_list = adata.varm[gp_targets_mask_key].sum(axis=0)
        n_sources_list = adata.varm[gp_sources_mask_key].sum(axis=0)

    
    # Convert the arrays to a pandas DataFrame
    targets_df = pd.DataFrame({"values": n_targets_list})
    sources_df = pd.DataFrame({"values": n_sources_list})

    # Determine plot configurations
    max_n_targets = max(n_targets_list)
    max_n_sources = max(n_sources_list)
    if max_n_targets > 200:
        targets_x_ticks_range = 100
        xticklabels_rotation = 45  
    elif max_n_targets > 100:
        targets_x_ticks_range = 20
        xticklabels_rotation = 0
    elif max_n_targets > 10:
        targets_x_ticks_range = 10
        xticklabels_rotation = 0
    else:
        targets_x_ticks_range = 1
        xticklabels_rotation = 0
    if max_n_sources > 200:
        sources_x_ticks_range = 100   
    elif max_n_sources > 100:
        sources_x_ticks_range = 20
    elif max_n_sources > 10:
        sources_x_ticks_range = 10
    else:
        sources_x_ticks_range = 1

    # Create subplot
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 5))
    plt.suptitle(
        f"{gp_plot_label} Gene Programs – Gene Count Distribution Plots")
    sns.histplot(x="values", data=targets_df, ax=ax1)
    ax1.set_title("Gene Program Targets Distribution",
                  fontsize=10)
    ax1.set(xlabel="Number of Targets",
            ylabel="Number of Gene Programs")
    ax1.set_xticks(
        np.arange(0,
                  max_n_targets + targets_x_ticks_range,
                  targets_x_ticks_range))
    ax1.set_xticklabels(
        np.arange(0,
                  max_n_targets + targets_x_ticks_range,
                  targets_x_ticks_range),
        rotation=xticklabels_rotation)
    sns.histplot(x="values", data=sources_df, ax=ax2)
    ax2.set_title("Gene Program Sources Distribution",
                  fontsize=10)
    ax2.set(xlabel="Number of Sources",
            ylabel="Number of Gene Programs")
    ax2.set_xticks(
        np.arange(0,
                  max_n_sources + sources_x_ticks_range,
                  sources_x_ticks_range))
    ax2.set_xticklabels(
        np.arange(0,
                  max_n_sources + sources_x_ticks_range,
                  sources_x_ticks_range),
        rotation=xticklabels_rotation)
    plt.subplots_adjust(wspace=0.35)
    if save_path:
        plt.savefig(save_path)
    plt.show()


# In[28]:


# Retrieve MEBOCOST GPs (source: enzyme genes; target: sensor genes)
mebocost_gp_dict = extract_gp_dict_from_mebocost_ms_interactions(
    dir_path=f"{figure_folder_path}" "/metabolite_enzyme_sensor_gps/",
    species='human',
    plot_gp_gene_count_distributions=True)


# In[27]:


# Display example MEBOCOST GP
mebocost_gp_names = list(mebocost_gp_dict.keys())
random.shuffle(mebocost_gp_names)
mebocost_gp_name = mebocost_gp_names[0]
print(f"{mebocost_gp_name}: {mebocost_gp_dict[mebocost_gp_name]}")


# In[29]:


def filter_and_combine_gp_dict_gps_v2(
        gp_dicts: list,
        overlap_thresh_target_genes: float=1.,
        verbose: bool=False) -> dict:
    """
    Combine gene program dictionaries and filter them based on gene overlaps.

    Parameters
    ----------
    gp_dicts:
        List of gene program dictionaries with keys being gene program names and
        values being dictionaries with keys ´sources´, ´targets´,
        ´sources_categories´, and ´targets_categories´, where ´targets´ contains
        a list of the names of genes in the gene program for the reconstruction
        of the gene expression of the node itself (receiving node) and ´sources´
        contains a list of the names of genes in the gene program for the
        reconstruction of the gene expression of the node's neighbors
        (transmitting nodes).
    overlap_thresh_target_genes:
        The minimum ratio of target genes that need to overlap between a GP
        without source genes and another GP for the GP to be dropped.
        Gene programs with different source genes are never combined or dropped.
    verbose:
        If `True`, print gene programs that are dropped and combined.

    Returns
    ----------
    new_gp_dict:
        Combined gene program dictionary with filtered gene programs.
    """
    # Combine gene program dictionaries
    combined_gp_dict = {}
    for i, gp_dict in enumerate(gp_dicts):
        combined_gp_dict.update(gp_dict)

    new_gp_dict = combined_gp_dict.copy()

    # Combine gene programs with overlapping genes
    all_combined = False
    while not all_combined:
        all_combined = True
        combined_gp_dict = new_gp_dict.copy()
        for i, (gp_i, gp_genes_dict_i) in enumerate(combined_gp_dict.items()):
            source_genes_i = [
                gene for gene in gp_genes_dict_i["sources"]]
            target_genes_i = [
                gene for gene in gp_genes_dict_i["targets"]]
            target_genes_categories_i = [
                target_gene_category for target_gene_category in
                gp_genes_dict_i["targets_categories"]]
            for j, (gp_j, gp_genes_dict_j) in enumerate(
                combined_gp_dict.items()):
                if j != i:
                    source_genes_j = [
                        gene for gene in gp_genes_dict_j["sources"]]
                    target_genes_j = [
                        gene for gene in gp_genes_dict_j["targets"]]
                    target_genes_categories_j = [
                        target_gene_category for target_gene_category in
                        gp_genes_dict_j["targets_categories"]]

                    if ((source_genes_i == source_genes_j) &
                        len(source_genes_i) > 0):
                        # if source genes are exactly the same, combine gene
                        # programs
                        all_combined = False
                        if verbose:
                            print(f"Combining {gp_i} and {gp_j}.")
                        source_genes = source_genes_i
                        target_genes = target_genes_i
                        target_genes_categories = target_genes_categories_i
                        for target_gene, target_gene_category in zip(
                            target_genes_j, target_genes_categories_j):
                            if target_gene not in target_genes:
                                target_genes.extend([target_gene])
                                target_genes_categories.extend(
                                    [target_gene_category])
                        new_gp_dict.pop(gp_i, None)
                        new_gp_dict.pop(gp_j, None)
                        if (gp_j.split("_")[0] + 
                            "_combined_GP") not in new_gp_dict.keys():
                            new_gp_name = gp_i.split("_")[0] + "_combined_GP"
                            new_gp_dict[new_gp_name] = {"sources": source_genes}
                            new_gp_dict[new_gp_name]["targets"] = target_genes
                            new_gp_dict[new_gp_name][
                                "sources_categories"] = gp_genes_dict_i[
                                    "sources_categories"]
                            new_gp_dict[new_gp_name][
                                "targets_categories"] = target_genes_categories
                            
                    elif len(source_genes_i) == 0:
                        target_genes_overlap = list(
                            set(target_genes_i) & set(target_genes_j))
                        n_target_gene_overlap = len(target_genes_overlap)
                        n_target_genes = len(target_genes_i)
                        ratio_shared_target_genes = (n_target_gene_overlap /
                                                     n_target_genes)
                        if ratio_shared_target_genes >= overlap_thresh_target_genes:
                            # if source genes not existent and target genes
                            # overlap more than specified, drop gene program
                            if gp_j in new_gp_dict.keys():
                                if verbose:
                                    print(f"Dropping {gp_i}.")
                                new_gp_dict.pop(gp_i, None)
                    else:
                        # otherwise do not combine or drop gene programs
                        pass

    return new_gp_dict


# In[30]:


def extract_gp_dict_from_collectri_tf_network(
        species: ["mouse", "human"],
        tf_network_file_path: Optional[str]="collectri_tf_network.csv",
        load_from_disk: bool=False,
        save_to_disk: bool=False,
        plot_gp_gene_count_distributions: bool=True,
        gp_gene_count_distributions_save_path: Optional[str]=None) -> dict:
    """
    Retrieve 1072 mouse or 1186 human transcription factor (TF) target gene gene
    programs from CollecTRI via decoupler. CollecTRI is a comprehensive resource
    containing a curated collection of TFs and their transcriptional targets
    compiled from 12 different resources. This collection provides an increased
    coverage of transcription factors and a superior performance in identifying
    perturbed TFs compared to the DoRothEA network and other literature based
    GRNs see
    https://decoupler-py.readthedocs.io/en/latest/notebooks/dorothea.html).

    Parameters
    ----------
    species:
        Species for which the gene programs will be extracted.
    load_from_disk:
        If ´True´, the CollecTRI TF network will be loaded from disk instead of
        from the decoupler library.
    save_to_disk:
        If ´True´, the CollecTRI TF network will additionally be stored on disk.
        Only applies if ´load_from_disk´ is ´False´.
    plot_gp_gene_count_distributions:
        If ´True´, display the distribution of gene programs per number of
        source and target genes.
    gp_gene_count_distributions_save_path:
        Path of the file where the gene program gene count distribution plot
        will be saved if ´plot_gp_gene_count_distributions´ is ´True´.

    Returns
    ----------
    gp_dict:
        Nested dictionary containing the CollecTRI TF target genes gene programs
        with keys being gene program names and values being dictionaries with
        keys ´sources´, ´targets´, ´sources_categories´, and
        ´targets_categories´, where ´sources´ and ´targets´ contain the
        CollecTRI TFs and target genes, and ´sources_categories´ and
        ´targets_categories´ contain the categories of all genes ('tf' or
        'target_gene').
    """
    if not load_from_disk:
        net = dc.get_collectri(organism=species, split_complexes=False)
        if save_to_disk:
            net.to_csv(tf_network_file_path, index=False)
    else:
        net = pd.read_csv(tf_network_file_path)

    tf_target_genes_df = net[["source", "target"]].groupby(
        "source")["target"].agg(list).reset_index()
    
    gp_dict = {}
    for tf, target_genes in zip(tf_target_genes_df["source"],
                                tf_target_genes_df["target"]):
        gp_dict[tf + "_TF_target_genes_GP"] = {
            "sources": [],
            "targets": [tf] + target_genes,
            "sources_categories": [],
            "targets_categories": ["tf"] + ["target_gene"] * len(target_genes)}
        
    if plot_gp_gene_count_distributions:
        create_gp_gene_count_distribution_plots(
            gp_dict=gp_dict,
            gp_plot_label="CollecTRI",
            save_path=gp_gene_count_distributions_save_path)
        
    return gp_dict


# In[31]:


# Retrieve CollecTRI GPs (source: -; target: transcription factor genes, target genes)

import decoupler as dc

collectri_gp_dict = extract_gp_dict_from_collectri_tf_network(
        species='human',
        tf_network_file_path=f"{figure_folder_path}" "/collectri_tf_network_human.csv",
        load_from_disk=False,
        save_to_disk=True,
        plot_gp_gene_count_distributions=True)


# In[31]:


# Display example CollecTRI GP
collectri_gp_names = list(collectri_gp_dict.keys())
random.shuffle(collectri_gp_names)
collectri_gp_name = collectri_gp_names[0]
print(f"{collectri_gp_name}: {collectri_gp_dict[collectri_gp_name]}")


# In[32]:


# Filter and combine GPs
gp_dicts = [omnipath_gp_dict, nichenet_gp_dict, mebocost_gp_dict,collectri_gp_dict]
combined_gp_dict = filter_and_combine_gp_dict_gps_v2(
    gp_dicts,
    verbose=True)

print(f"Number of gene programs after filtering and combining: "
      f"{len(combined_gp_dict)}.")


# ### Load Data & Compute Spatial Neighbor Graph

# In[33]:


# Compute spatial neighborhood
sq.gr.spatial_neighbors(adata,
                        coord_type="generic",
                        spatial_key='spatial')

# Make adjacency matrix symmetric
adata.obsp['spatial_connectivities'] = (
    adata.obsp['spatial_connectivities'].maximum(
        adata.obsp['spatial_connectivities'].T))


# ### Filter Genes

# In[34]:


def get_unique_genes_from_gp_dict(
        gp_dict: dict,
        retrieved_gene_entities: list=["sources", "targets"],
        retrieved_gene_categories: Optional[list]=None) -> list:
    """
    Return all unique genes of a gene program dictionary.

    Parameters
    ----------
    gp_dict:
        The gene program dictionary from which to retrieve the unique genes.
    retrieved_gene_entities:
        A list that contains all gene entities ("sources", "targets")
        for which unique genes of the gene program dictionary should be
        retrieved.
    retrieved_gene_categories:
        A list that contains all gene categories for which unique genes of the
        gene program dictionary should be retrieved. If `None`, all gene
        categories are included.

    Returns
    ----------
    unique_genes:
        A list of unique genes used in the gene program dictionary.
    """
    gene_list = []

    for _, gp in gp_dict.items():
        for gene_entity in retrieved_gene_entities:
            genes = gp[gene_entity]
            gene_categories = gp[f"{gene_entity}_categories"]
            if retrieved_gene_categories is not None:
                genes = [gene for gene, gene_category in zip(genes, gene_categories) if
                         gene_category in retrieved_gene_categories]
            gene_list.extend(genes)
    unique_genes = list(set(gene_list))
    unique_genes.sort()
    return unique_genes


# In[35]:


filter_genes = True
n_svg = 955

if filter_genes:
    print("Filtering genes...")
    # Filter genes and only keep ligand, receptor, enzyme, sensor, and
    # the 'n_svg' spatially variable genes
    gp_dict_genes = get_unique_genes_from_gp_dict(
        gp_dict=combined_gp_dict,
            retrieved_gene_entities=["sources", "targets"])
    print(f"Starting with {len(adata.var_names)} genes.")
    adata.var["sources_targets"] = adata.var_names.isin(gp_dict_genes)
    print(f"Keeping {len(adata.var_names)} ligand, receptor, enzyme, sensor genes.")
    
    # Identify spatially variable genes
    sq.gr.spatial_autocorr(adata, mode="moran", genes=adata.var_names)
    svg_genes = adata.uns["moranI"].index[:n_svg].tolist()
    adata.var["spatially_variable_all"] = adata.var_names.isin(svg_genes)
    adata = adata[:, adata.var["spatially_variable_all"] == True]
    print(f"Keeping {len(adata.var_names)} spatially variable genes.")


# ### Add GP Mask to Data

# In[36]:


# Add the GP dictionary as binary masks to the adata
add_gps_from_gp_dict_to_adata(
    gp_dict=combined_gp_dict,
    adata=adata,
    gp_targets_mask_key='nichecompass_gp_targets',
    gp_targets_categories_mask_key='nichecompass_gp_targets_categories',
    gp_sources_mask_key='nichecompass_gp_sources',
    gp_sources_categories_mask_key='nichecompass_gp_sources_categories',
    gp_names_key='nichecompass_gp_names',
    min_genes_per_gp=2,
    min_source_genes_per_gp=1,
    min_target_genes_per_gp=1,
    max_genes_per_gp=None,
    max_source_genes_per_gp=None,
    max_target_genes_per_gp=None)


# ### explore data

# In[37]:


cell_type_colors = create_new_color_dict(
    adata=adata,
    cat_key='new_leiden')

adata.layers['counts']=adata.X.copy()

print(f"Number of nodes (observations): {adata.layers['counts'].shape[0]}")
print(f"Number of gene node features: {adata.layers['counts'].shape[1]}")

# Visualize cell-level annotated data in physical space
sc.pl.spatial(adata,
              color='new_leiden',
              palette=cell_type_colors,
              spot_size=1, library_id = '21')        


# ### Model training

# In[38]:


#Initialize Model

adata.layers['counts']=adata.X.copy()
counts_key = "counts"
adj_key = "spatial_connectivities"
gp_names_key = "nichecompass_gp_names"
active_gp_names_key = "nichecompass_active_gp_names"
gp_targets_mask_key = "nichecompass_gp_targets"
gp_targets_categories_mask_key = "nichecompass_gp_targets_categories"
gp_sources_mask_key = "nichecompass_gp_sources"
gp_sources_categories_mask_key = "nichecompass_gp_sources_categories"
latent_key = "nichecompass_latent"
conv_layer_encoder = "gcnconv"
active_gp_thresh_ratio = 0.01

model = NicheCompass(adata,
                     counts_key=counts_key,
                     adj_key=adj_key,
                     gp_names_key=gp_names_key,
                     active_gp_names_key=active_gp_names_key,
                     gp_targets_mask_key=gp_targets_mask_key,
                     gp_targets_categories_mask_key=gp_targets_categories_mask_key,
                     gp_sources_mask_key=gp_sources_mask_key,
                     gp_sources_categories_mask_key=gp_sources_categories_mask_key,
                     latent_key=latent_key,
                     conv_layer_encoder=conv_layer_encoder,
                     active_gp_thresh_ratio=active_gp_thresh_ratio)


# In[39]:


adata.layers['counts'] = adata.layers['counts'].astype(np.float32)


# In[40]:


# Train model
n_epochs = 400
n_epochs_all_gps = 40
lr = 0.001
lambda_edge_recon = 500000.
lambda_gene_expr_recon = 300.
lambda_l1_masked = 0. # prior GP  regularization
lambda_l1_addon = 30. # de novo GP regularization
edge_batch_size = 1024 # increase if more memory available or decrease to save memory
n_sampled_neighbors = 10
use_cuda_if_available = True

model.train(n_epochs=n_epochs,
            n_epochs_all_gps=n_epochs_all_gps,
            lr=lr,
            lambda_edge_recon=lambda_edge_recon,
            lambda_gene_expr_recon=lambda_gene_expr_recon,
            lambda_l1_masked=lambda_l1_masked,
            lambda_l1_addon=lambda_l1_addon,
            edge_batch_size=edge_batch_size,
            use_cuda_if_available=use_cuda_if_available,
            n_sampled_neighbors=n_sampled_neighbors,
            verbose=True)


# In[41]:


# Save trained model
model.save(dir_path=f"{figure_folder_path}",
           overwrite=True,
           save_adata=True,
           adata_file_name="adata_before.h5ad")


# In[42]:


# Compute latent neighbor graph
sc.pp.neighbors(model.adata,
                use_rep=latent_key,
                key_added=latent_key,random_state=seed)

# Compute UMAP embedding
sc.tl.umap(model.adata,
           neighbors_key=latent_key,random_state=seed)


# In[43]:


# Save trained model
model.save(dir_path=f"{figure_folder_path}",
           overwrite=True,
           save_adata=True,
           adata_file_name="adata_after.h5ad")


# ### Visualize NicheCompass Latent GP Space
# 
# Nel caso delle reti neurali, i modelli  addestrati producono come risultato uno spazio vettoriale ovvero  una rappresentazione compressa dell’informazione. La versione compressa della distribuzione dei dati viene definita spazio latente

# In[313]:


model_folder_path = '/Users/paolamaragno/Downloads/NicheCompass'
gp_names_key = "nichecompass_gp_names"

model = NicheCompass.load(dir_path=model_folder_path,
                          adata=None,
                          adata_file_name="adata_after.h5ad",
                          adata_atac=None,
                          gp_names_key=gp_names_key)


# In[314]:


samples = model.adata.obs['fov'].unique().tolist()

# Note that the goal of NicheCompass is not a separation of cell types but rather to identify spatially consistent cell niches.

cell_type_colors = create_new_color_dict(
    adata=model.adata,
    cat_key='new_leiden')


# In[197]:


"""
# Create plot of cell type annotations in physical and latent space
from matplotlib import gridspec

groups = None
save_fig = True
file_path = f"{figure_folder_path}/" \
            "cell_types_latent_physical_space.pdf"
cell_type_key = 'new_leiden'
sample_key = 'fov'
spot_size = 1

fig = plt.figure(figsize=(100, 100))
title = fig.suptitle(t="Cell Types in Latent and Physical Space",
                     y=0.96,
                     x=0.55,
                     fontsize=20)
spec1 = gridspec.GridSpec(ncols=1,
                          nrows=2,
                          width_ratios=[1],
                          height_ratios=[3, 2])
spec2 = gridspec.GridSpec(ncols=len(samples),
                          nrows=2,
                          width_ratios=[1] * len(samples),
                          height_ratios=[3, 2])
axs = []
axs.append(fig.add_subplot(spec1[0]))
sc.pl.umap(adata=model.adata,
           color=[cell_type_key],
           groups=groups,palette=cell_type_colors,
           title=f"Cell Types in Latent Space",
           ax=axs[0],
           show=False)
for idx, sample in enumerate(samples):
    axs.append(fig.add_subplot(spec2[len(samples) + idx]))
    sc.pl.spatial(adata=model.adata[model.adata.obs[sample_key] == sample],
                  color=[cell_type_key],
                  groups=groups,
                  palette=cell_type_colors,
                  spot_size=spot_size,
                  title=f"Cell Types in Physical Space \n"
                        f"(Sample: {sample})",
                  legend_loc=None,
                  ax=axs[idx+1],
                  show=False, library_id=sample)

# Create and position shared legend
handles, labels = axs[0].get_legend_handles_labels()
lgd = fig.legend(handles,
                 labels,
                 loc="center left",
                 bbox_to_anchor=(0.98, 0.5))
axs[0].get_legend().remove()

# Adjust, save and display plot
plt.subplots_adjust(wspace=0.2, hspace=0.25)
if save_fig:
    fig.savefig(file_path,
                bbox_extra_artists=(lgd, title),
                bbox_inches="tight")
plt.show()
"""


# ### Identify Niches

# In[315]:


# We compute Leiden clustering of the NicheCompass latent GP space to identify spatially consistent cell niches.

latent_leiden_resolution = 0.5
latent_cluster_key = f"latent_leiden_{str(latent_leiden_resolution)}"
latent_key = "nichecompass_latent"

# Compute latent Leiden clustering
sc.tl.leiden(adata=model.adata,
             resolution=latent_leiden_resolution,
             key_added=latent_cluster_key,
             neighbors_key=latent_key,random_state=seed)


# In[316]:


model.adata.obs['latent_leiden_0.5']


# In[317]:


palette = ['#0F2080', '#006600', '#AA3377', '#BBBBBB', '#AEC7E8', '#CCBB44', '#F7B6D2','#000000',
          '#0077BB', '#EE7733', '#33BBEE', '#CC3311', '#7F7F7F', '#FFB000','#601A4A','#785EF0','#66CC00', '#CCFF99', 
          '#660000', '#E5CCFF','#FFFFCC']

sns.color_palette(palette)


# In[318]:


sc.pl.umap(adata=model.adata,
           color=[latent_cluster_key],
           palette=sns.color_palette(palette),
           title=f"Niches in Latent Space",
           show=False)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/umap_niches.pdf',bbox_inches="tight")


# In[44]:


t = ['fov 21 - Patient_11 - LOW',
 'fov 23 - Patient_12 - HIGH',
 'fov 24 - Patient_12 - HIGH',
 'fov 25 - Patient_13 - LOW',
 'fov 27 - Patient_14 - HIGH',
 'fov 28 - Patient_14 - HIGH',
 'fov 29 - Patient_15 - LOW',
 'fov 30 - Patient_15 - LOW',
 'fov 31 - Patient_16 - HIGH',
 'fov 32 - Patient_16 - HIGH',
 'fov 33 - Patient_17 - HIGH',
 'fov 34 - Patient_17 - HIGH',
 'fov 35 - Patient_18 - HIGH',
 'fov 36 - Patient_18 - HIGH',
 'fov 37 - Patient_19 - LOW', 
 'fov 38 - Patient_19 - LOW',
 'fov 39 - Patient_20 - LOW',
 'fov 40 - Patient_20 - LOW',
 'fov 41 - Patient_20 - LOW']

sq.pl.spatial_segment(
    adata=model.adata,
    color=[latent_cluster_key],
    library_key="fov",
    seg_cell_id="cell_ID",
    seg_outline=True,
    img=False,
    scalebar_dx=1.0,
    title= t,
    scalebar_kwargs={"scale_loc": "bottom", "location": "lower right"},
)

plt.savefig('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/niche_all_fov.pdf',bbox_inches="tight")


# ### Niche composition

# In[319]:


palette = ['#0F2080', '#228833', '#AA3377', '#BBBBBB', '#AEC7E8', '#CCBB44', '#F7B6D2',
          '#0077BB', '#EE7733', '#33BBEE', '#CC3311', '#7F7F7F']


# In[320]:


# We can analyze the niche composition in terms of cell type labels.

save_fig = True
file_path = f"/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/niche_composition.pdf"
cell_type_key = 'new_leiden'

df_counts = (model.adata.obs.groupby([latent_cluster_key, cell_type_key])
             .size().unstack())
df_counts.plot(kind="bar", stacked=True, figsize=(10,10), color = palette)
legend = plt.legend(bbox_to_anchor=(1, 1), loc="upper left", prop={'size': 10})
legend.set_title("Cell Type Annotations", prop={'size': 10})
plt.title("Cell Type Composition of Niches")
plt.xlabel("Niche")
plt.ylabel("Cell Type Counts")
if save_fig:
    plt.savefig(file_path,
                bbox_extra_artists=(legend,),
                bbox_inches="tight")


# ### Differential GPs

# In[304]:


# Check number of active GPs

active_gps = model.get_active_gps()
print(f"Number of total gene programs: {len(model.adata.uns[gp_names_key])}.")
print(f"Number of active gene programs: {len(active_gps)}.")


# In[135]:


model.adata.obsm['nichecompass_latent'] # for each cell, score for each active GP 
model.adata.obsm['nichecompass_latent'].shape


# In[295]:


model.adata.obsm['nichecompass_latent']


# In[305]:


# Display example active GPs
gp_summary_df = model.get_gp_summary()
gp_summary_df[gp_summary_df["gp_active"] == True].head()


# In[321]:


# Set parameters for differential gp testing
selected_cats = None
comparison_cats = "rest"
title = f"NicheCompass Strongly Enriched Niche GPs"
log_bayes_factor_thresh = 2.3
save_fig = True
file_path = f"/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/niches_enriched_gps_heatmap.pdf"


# In[322]:


# Run differential gp testing
enriched_gps = model.run_differential_gp_tests(
    cat_key=latent_cluster_key,
    selected_cats=selected_cats,
    comparison_cats=comparison_cats,
    log_bayes_factor_thresh=log_bayes_factor_thresh)


# In[308]:


# Results are stored in a df in the adata object
differential_gp_test_results_key = "nichecompass_differential_gp_test_results"

model.adata.uns[differential_gp_test_results_key]


# In[129]:


len(enriched_gps)


# In[145]:


pd.DataFrame(model.adata.obsm['nichecompass_latent'], index = model.adata.obs.index.values, columns = active_gps)['CD209_combined_GP']


# In[339]:


model.adata.obs[[latent_cluster_key] + enriched_gps].groupby(latent_cluster_key).mean()


# In[337]:


m = 
m['fov'] = model.adata.obs['fov']


# In[338]:


m.to_csv("/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/score_enriched_GP.csv")


# In[253]:


'Il10_combined_GP' in enriched_gps


# In[323]:


# Visualize GP activities of enriched GPs across niches

from sklearn.preprocessing import MinMaxScaler

# per ogni cellula vedi qual è lo score che le è stato dato relativamente ad uno dei GP arricchiti, 
# raggruppi per nicchia e calcoli la media. quindi hai per ogni nicchia la media dello score relativo ad un GP considerando tutte le 
# cellule di quella nicchia
df = model.adata.obs[[latent_cluster_key] + enriched_gps].groupby(latent_cluster_key).mean()

scaler = MinMaxScaler()
normalized_columns = scaler.fit_transform(df)
normalized_df = pd.DataFrame(normalized_columns, columns=df.columns)
normalized_df.index = df.index
# normalizzazione 

plt.figure(figsize=(20, 12))  # Set the figure size
ax = sns.clustermap(normalized_df,
            cmap='viridis')
plt.setp(ax.ax_heatmap.xaxis.get_majorticklabels(), fontsize=8) 
ax.cax.remove() 
plt.savefig(f"/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/enriched_gps_heatmap.pdf")


# In[324]:


# Store gene program summary of enriched gene programs
save_file = True
file_path = f"/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/niche_enriched_gps_summary.csv"

gp_summary_cols = ["gp_name",
                   "n_source_genes",
                   "n_non_zero_source_genes",
                   "n_target_genes",
                   "n_non_zero_target_genes",
                   "gp_source_genes",
                   "gp_target_genes",
                   "gp_source_genes_importances",
                   "gp_target_genes_importances"]

enriched_gp_summary_df = gp_summary_df[gp_summary_df["gp_name"].isin(enriched_gps)]
cat_dtype = pd.CategoricalDtype(categories=enriched_gps, ordered=True)
enriched_gp_summary_df.loc[:, "gp_name"] = enriched_gp_summary_df["gp_name"].astype(cat_dtype)
enriched_gp_summary_df = enriched_gp_summary_df.sort_values(by="gp_name")
enriched_gp_summary_df = enriched_gp_summary_df[gp_summary_cols]

if save_file:
    enriched_gp_summary_df.to_csv(f"{file_path}")
else:
    display(enriched_gp_summary_df)


# In[311]:


latent_cluster_colors = create_new_color_dict(
    adata=model.adata,
    cat_key=latent_cluster_key)


# In[62]:


# Now we will have a look at the GP activities and the log normalized counts of the most important omics features of the differential GPs.

plot_label = f"log_bayes_factor_{log_bayes_factor_thresh}_cluster_{selected_cats[0] if selected_cats else 'None'}_vs_rest"
save_figs = True
sample_key = 'fov'

generate_enriched_gp_info_plots(
    plot_label=plot_label,
    model=model,
    spot_size=0.2,
    sample_key=sample_key,
    differential_gp_test_results_key=differential_gp_test_results_key,
    cat_key=latent_cluster_key,
    cat_palette=latent_cluster_colors,
    n_top_enriched_gp_start_idx=0,
    n_top_enriched_gp_end_idx=10,
    n_top_genes_per_gp=3,
    save_figs=save_figs,
figure_folder_path=f"{figure_folder_path}/")


# ### Cell-cell Communication

# In[327]:


gp_summary_df = model.get_gp_summary()
gp = 'IL11_combined_GP'
gp_idx = model.adata.uns[model.gp_names_key_].tolist().index(gp)
active_gp_idx = model.adata.uns[model.active_gp_names_key_].tolist().index(gp)
gp_scores = model.adata.obsm[model.latent_key_][:, active_gp_idx]
gp_targets_cats = model.adata.varm[model.gp_targets_categories_mask_key_][:, gp_idx] # quali geni sono i target di questo GP
gp_sources_cats = model.adata.varm[model.gp_sources_categories_mask_key_][:, gp_idx] # quali geni sono i sources di questo GP
targets_cats_label_encoder = model.adata.uns[model.targets_categories_label_encoder_key_] 
sources_cats_label_encoder = model.adata.uns[model.sources_categories_label_encoder_key_]
sources_cat_idx_dict = {}
for source_cat, source_cat_label in sources_cats_label_encoder.items():
            sources_cat_idx_dict[source_cat] = np.where(gp_sources_cats == source_cat_label)[0]
targets_cat_idx_dict = {}
for target_cat, target_cat_label in targets_cats_label_encoder.items():
            targets_cat_idx_dict[target_cat] = np.where(gp_targets_cats == target_cat_label)[0]

# Get indices of all source and target genes
source_genes_idx = np.array([], dtype=np.int64)
for key in sources_cat_idx_dict.keys():
            source_genes_idx = np.append(source_genes_idx,
                                         sources_cat_idx_dict[key])
target_genes_idx = np.array([], dtype=np.int64)
for key in targets_cat_idx_dict.keys():
            target_genes_idx = np.append(target_genes_idx,
                                         targets_cat_idx_dict[key])
            
gp_source_scores = np.zeros((len(model.adata.obs), len(source_genes_idx)))
# array with rows equal to cells and columns equal to the number of source genes 
gp_target_scores = np.zeros((len(model.adata.obs), len(target_genes_idx)))
# array with rows equal to cells and columns equal to the number of target genes 

# per ogni cellula, per ogni source gene metti un valore pari a: 
# valore di espressione di quel gene sorgente in ogni cellula / il massimo di espressione di quel gene * il weight di quel gene * il gp score di quella cellula
for i, source_gene_idx in enumerate(source_genes_idx):
            source_gene = model.adata.var_names[source_gene_idx]
            gp_source_scores[:, i] = (
                model.adata[:, model.adata.var_names.tolist().index(source_gene)].X.toarray().flatten() / model.adata[:, model.adata.var_names.tolist().index(source_gene)].X.toarray().flatten().max() *
                gp_summary_df[gp_summary_df["gp_name"] == gp]["gp_source_genes_weights"].values[0][gp_summary_df[gp_summary_df["gp_name"] == gp]["gp_source_genes"].values[0].index(source_gene)] *
                gp_scores)

for j, target_gene_idx in enumerate(target_genes_idx):
            target_gene = model.adata.var_names[target_gene_idx]
            gp_target_scores[:, j] = (
                model.adata[:, model.adata.var_names.tolist().index(target_gene)].X.toarray().flatten() / model.adata[:, model.adata.var_names.tolist().index(target_gene)].X.toarray().flatten().max() *
                gp_summary_df[gp_summary_df["gp_name"] == gp]["gp_target_genes_weights"].values[0][gp_summary_df[gp_summary_df["gp_name"] == gp]["gp_target_genes"].values[0].index(target_gene)] *
                gp_scores)

agg_gp_source_score = gp_source_scores.mean(1).astype("float32")
agg_gp_target_score = gp_target_scores.mean(1).astype("float32")
agg_gp_source_score[agg_gp_source_score < 0] = 0.
agg_gp_target_score[agg_gp_target_score < 0] = 0.

model.adata.obs[f"{gp}_source_score"] = agg_gp_source_score
model.adata.obs[f"{gp}_target_score"] = agg_gp_target_score
        
del(gp_target_scores)
del(gp_source_scores)

agg_gp_source_score = sp.sparse.csr_matrix(agg_gp_source_score)
agg_gp_target_score = sp.sparse.csr_matrix(agg_gp_target_score)

model.adata.obsp[f"{gp}_connectivities"] = (model.adata.obsp["spatial_cci_connectivities"] > 0).multiply(
            agg_gp_source_score.T.dot(agg_gp_target_score))

group_key=latent_cluster_key
filter_key = None
gp_network_df_pivoted = aggregate_obsp_matrix_per_cell_type(
            adata=model.adata,
            obsp_key=f"{gp}_connectivities",
            cell_type_key=group_key,
            group_key=filter_key,
            agg_rows=True)

gp_network_df = gp_network_df_pivoted.melt(var_name="source", value_name="gp_score", ignore_index=False).reset_index()
gp_network_df.columns = ["source", "target", "strength"]


# In[263]:


import inspect
lines = inspect.getsource(visualize_communication_gp_network)
print(lines)


# In[297]:


from typing import Optional

def aggregate_obsp_matrix_per_cell_type(
        adata,
        obsp_key: str,
        cell_type_key: str="cell_type",
        group_key: Optional[str]=None,
        agg_rows: bool=False):
    """
    Generic function to aggregate adjacency matrices stored in
    ´adata.obsp[obsp_key]´ on cell type level. It can be used to aggregate the
    node label aggregator aggregation weights alpha or the reconstructed adjacency
    matrix of a trained NicheCompass model by neighbor cell type for downstream
    analysis.

    Parameters
    ----------
    adata:
        AnnData object which contains outputs of NicheCompass model training.
    obsp_key:
        Key in ´adata.obsp´ where the matrix to be aggregated is stored.
    cell_type_key:
        Key in ´adata.obs´ where the cell type labels are stored.
    group_key:
        Key in ´adata.obs´ where additional grouping labels are stored.    
    agg_rows:
        If ´True´, also aggregate over the observations on cell type level.

    Returns
    ----------
    cell_type_agg_df:
        Pandas DataFrame with the aggregated obsp values (dim: n_obs x
        n_cell_types if ´agg_rows == False´, else n_cell_types x n_cell_types).
    """
    n_obs = len(adata)
    n_cell_types = adata.obs[cell_type_key].nunique()
    sorted_cell_types = sorted(adata.obs[cell_type_key].unique().tolist())

    cell_type_label_encoder = {k: v for k, v in zip(
        sorted_cell_types,
        range(n_cell_types))}

    # Retrieve non zero indices and non zero values, and create row-wise
    # observation cell type index
    nz_obsp_idx = adata.obsp[obsp_key].nonzero()
    neighbor_cell_type_index = adata.obs[cell_type_key][nz_obsp_idx[1]].map(
        cell_type_label_encoder).values
    adata.obsp[obsp_key].eliminate_zeros() # In some sparse reps 0s can appear
    nz_obsp = adata.obsp[obsp_key].data

    # Use non zero indices, non zero values and row-wise observation cell type
    # index to construct new df with cell types as columns and row-wise
    # aggregated values per cell type index as values
    cell_type_agg = np.zeros((n_obs, n_cell_types))
    np.add.at(cell_type_agg,
              (nz_obsp_idx[0], neighbor_cell_type_index),
              nz_obsp)
    cell_type_agg_df = pd.DataFrame(
        cell_type_agg,
        columns=sorted_cell_types)
    
    # Add cell type labels of observations
    cell_type_agg_df[cell_type_key] = adata.obs[cell_type_key].values

    # If specified, add group label
    if group_key is not None:
        cell_type_agg_df[group_key] = adata.obs[group_key].values

    if agg_rows:
        # In addition, aggregate values across rows to get a
        # (n_cell_types x n_cell_types) df
        if group_key is not None:
                       cell_type_agg_df = cell_type_agg_df.groupby(
                [group_key, cell_type_key]).sum()
        else:
            cell_type_agg_df = cell_type_agg_df.groupby(cell_type_key).sum()

        # Sort index to have same order as columns
        cell_type_agg_df = cell_type_agg_df.loc[
            sorted(cell_type_agg_df.index.tolist()), :]
        
    return cell_type_agg_df


# In[326]:


# Now we will use the inferred activity of an enriched combined interaction GP to analyze the involved intercellular interactions.

gp_name = "IL11_combined_GP"

network_df = compute_communication_gp_network(
    gp_list=[gp_name],
    model=model,
    group_key=latent_cluster_key,
    n_neighbors=4)

visualize_communication_gp_network(
    adata=model.adata,
    network_df=network_df,
    figsize=(9, 8),
    cat_colors=latent_cluster_colors,
    edge_type_colors=["#1f77b4"], 
    cat_key=latent_cluster_key,
    save=True,
    save_path=f"/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/gp_network_{gp_name}.pdf",
    )


# In[267]:


model.adata.obs.to_csv('/Users/paolamaragno/Library/CloudStorage/OneDrive-UniversitàdegliStudidiMilano/Policlinico/analisi/spatial/spatial/CD/results_seed/metadata_adata_after_niche.csv')

