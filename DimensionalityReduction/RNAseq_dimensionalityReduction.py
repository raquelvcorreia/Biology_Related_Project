import pandas as pd
from scipy import stats
import umap
import hdbscan
import sklearn.cluster as cluster
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
import numpy as np
import seaborn as sns
from sklearn.manifold import TSNE


## just to be able to check more columns
desired_width=320
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width)
pd.set_option('display.max_columns',10)


df = pd.read_csv('GSE150316_DeseqNormCounts_final.txt.gz', compression='gzip', sep='\t', quotechar='"')


patient_data_df = pd.read_excel('2020.07.30.20165241-1_supp_table3.xlsx')
print(patient_data_df.head())

#reformat patient_data_df and create a new df with the reformated data

patient_data_df.columns = patient_data_df.columns.str.replace("\\r\\n|\\s", "_", regex = True)
print(patient_data_df.head())

patient_viral_load = patient_data_df.iloc[:,[0, 1, 2]].copy()
patient_viral_load.rename(columns ={'Case_No': 'case' }, inplace=True)
patient_viral_load['case'] = patient_viral_load['case'].astype(str)
print(patient_viral_load.head())


# filter data
# Genes where transcript levels are very low, or which don’t show a lot of
#variance, are not that informative for dimensionality reduction, and will just slow it down.
# Trim the tissue data down by removing the data on these genes
# calculate the mean across the rows, and store it in a new column mean, and remove the ones with a mean <= 0.5
# calculate the variance across the rows and add that information to a new column
# reorder the dataset by the variance level (desc)


df['mean'] = df.apply(lambda row: row.mean(), axis=1)

df = df[df['mean'] > 0.5]

df['variance'] = df.apply(lambda row: row.var(), axis=1)
df = df.sort_values(by=['variance'], ascending = False)
# converting the genes (now indexes) into a column, index to column
df.reset_index(inplace=True)
df = df.rename(columns = {'index':'gene'})

print(df.head())
print(df.tail())





##transpose the data on df to later add the information from patient_viral_load df

tissue_RNAseq = df.iloc[:, 0:89]
tissue_RNAseq: object = tissue_RNAseq.set_index('gene').T
tissue_RNAseq.insert (loc=0, column = 'details_case_tissue', value = tissue_RNAseq.index.values)
case_info = list(tissue_RNAseq.index.values)
case_info_sep = pd.DataFrame(case_info, columns = ['full_name'])

case_info_sep = case_info_sep['full_name'].str.split('-', expand= True)
case_info_sep['case'] = case_info_sep.iloc[:,0].str[4:]
case_info_sep['tissue_num'] = case_info_sep.iloc[:,1].str.extract('(\d+)', expand = False)
case_info_sep['tissue'] = case_info_sep.iloc[:,1].str.extract('([a-z]+)', expand = False)
case_info_sep.tissue_num.fillna(case_info_sep.iloc[:,2], inplace = True)
case_info_sep = case_info_sep.drop([0,1,2], axis = 1)
case_info_sep['details_case_tissue'] = pd.DataFrame(case_info, columns = ['full_name'])
case_info_sep = case_info_sep[['details_case_tissue','case', 'tissue', 'tissue_num']]


tissue_RNAseq_m = pd.merge(case_info_sep, tissue_RNAseq, on='details_case_tissue')
tissue_RNAseq_m['case'] = tissue_RNAseq_m['case'].astype(str)
tissue_RNAseq_m = tissue_RNAseq_m.merge(patient_viral_load, how = 'left', on = 'case')
# get column names from tissue_RNAseq_m
column_names = list(tissue_RNAseq_m)
# reorder
second_column = tissue_RNAseq_m.pop('Viral_load')
third_column = tissue_RNAseq_m.pop('Viral_high_vs._viral_low*')
tissue_RNAseq_m.insert(1, 'Viral_load', second_column)
tissue_RNAseq_m.insert(2, 'Viral_high_vs._viral_low*', third_column)
print(tissue_RNAseq_m.head())

# How many tissue types do we have? and how many samples from each tissue type
print('How many samples we have per tissue: ', tissue_RNAseq_m['tissue'].value_counts())

# will only be using the gene information for the dimensionality reduction
# get column names from tissue_RNAseq_m
column_names = list(tissue_RNAseq_m)
# Version of the data without the sample information
tissue_RNAseq_condensed = tissue_RNAseq_m.iloc[:,7:]
# sample information
sample_info = tissue_RNAseq_m.iloc[:,:7]
print(sample_info)


### UMAP
standard_embedding = umap.UMAP(random_state=42).fit_transform(tissue_RNAseq_condensed)


umap_df = pd.DataFrame(standard_embedding, columns = ['x_coord', 'y_coord'])

sample_umap_df = pd.concat([sample_info, umap_df], axis = 1)
sample_umap_df['tissue'].replace(to_replace=[None], value = 'negative_control', inplace = True)
sample_umap_df['Viral_high_vs._viral_low*'].fillna(value = 'DNW', inplace = True)
print(sample_umap_df)
print(sample_umap_df['tissue'].unique())

# Scatterplot
# plt.figure(figsize=(8,6))

colors = {'lung': 1, 'heart': 10, 'jejunum': 20, 'liver': 30, 'kidney': 40,'bowel': 50, 'fat': 60, 'skin':70, 'marrow': 80, 'negative_control': 90, 'placenta': 100}
colors_list = colors.keys()
sizes = {'High': 40, 'Low': 30, 'DNW' : 10}
size_list = sizes.keys()
scatter = plt.scatter(sample_umap_df.loc[:,'x_coord'], sample_umap_df.loc[:,'y_coord'], c=sample_umap_df['tissue'].map(colors),  cmap='viridis', s = sample_umap_df['Viral_high_vs._viral_low*'].map(sizes), alpha=0.5)
plt.title('UMAP plot in 2D')
plt.xlabel('umap component 1')
plt.ylabel('umap component 2')


# plt.legend(handles = scatter.legend_elements()[0], labels = colors.keys())
#
# fig, ax = plt.subplots()
# scatter = ax.scatter(sample_umap_df.loc[:,'x_coord'], sample_umap_df.loc[:,'y_coord'], c=sample_umap_df['tissue'].map(colors),  cmap='viridis', s = sample_umap_df['Viral_high_vs._viral_low*'].map(sizes), alpha=0.5)
#
# legend1 = ax.legend(*scatter.legend_elements(),
#                     loc="upper left", title="Tissues")
# ax.add_artist(legend1)
# handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
# legend2 = ax.legend(handles, labels,  loc="lower right", title="Viral Load")
# ax.set_title('UMAP plot in 2D')
# ax.set_xlabel('umap component 1')
# ax.set_ylabel('umap component 2')
#



plt.show()

### Clustering with Kmeans
kmeans_labels = cluster.KMeans(n_clusters=11).fit_predict(tissue_RNAseq_condensed)
plt.scatter(sample_umap_df.loc[:,'x_coord'], sample_umap_df.loc[:,'y_coord'], c=kmeans_labels, s= sample_umap_df['Viral_high_vs._viral_low*'].map(sizes), cmap='Spectral')
plt.show()


################## tSNE ########

tsne = TSNE(n_components=2)
rnaseq_tsne = tsne.fit_transform(tissue_RNAseq_condensed)

# Convert to data frame
tsne_df = pd.DataFrame(data = rnaseq_tsne, columns = ['tsne comp. 1', 'tsne comp. 2'])

plt.figure(figsize=(8,6))

# Scatterplot
plt.scatter(tsne_df.iloc[:,0], tsne_df.iloc[:,1], c=sample_umap_df['tissue'].map(colors),  cmap='viridis', s = sample_umap_df['Viral_high_vs._viral_low*'].map(sizes), alpha=0.5)

# Aesthetics
plt.title('t-SNE plot in 2D')
plt.xlabel('tsne component 1')
plt.ylabel('tsne component 2')
plt.show()