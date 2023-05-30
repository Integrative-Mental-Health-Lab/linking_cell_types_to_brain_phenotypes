'''
Get the dissection information of the Siletti RNAseq dataset.
The script produces two files:
    Ncells_from_Dissection_in_Cluster.csv: a matrix of size (num_clusters, num_dissections) denoting the number of cells from each dissection in a given cluster.
    dissection_info.csv: a list of size (num_clusters,) with each row containing the descending proportion of cells from each dissection in a given cluster.
'''
import h5py
import os
import numpy as np
from tqdm import tqdm

f = h5py.File('adult_human_20221007.loom', 'r')
df_cluster = f['col_attrs/Clusters'][:] 
df_roi = f['col_attrs/Roi'][:] # has 108 values but want 106

roi_list = np.unique(df_roi)
roi_list = np.delete(roi_list, [6,60])
roi_list_decode = np.array([i.decode() for i in roi_list])

n_roi = roi_list.shape[0]
n_clusters = np.max(df_cluster) + 1

clusters = np.array([("Cluster" + str(i)) for i in np.arange(n_clusters)])

df_cells = pd.DataFrame(0, index=clusters, columns=roi_list_decode)

for cluster_num in tqdm(np.arange(n_clusters)):
    for roi in roi_list:
        if roi == b'Human A25':
            df_cells.at["Cluster" + str(cluster_num), roi.decode()] = \
            np.sum(np.logical_and(df_cluster==cluster_num, np.logical_or(df_roi==roi, df_roi==b'Human A25\xe2\x80\x8b\xe2\x80\x8b')))
        elif roi == b'Human Idg':
            df_cells.at["Cluster" + str(cluster_num), roi.decode()] = \
            np.sum(np.logical_and(df_cluster==cluster_num, np.logical_or(df_roi==roi, df_roi==b'Human Idg\xe2\x80\x8b')))
        else:
            df_cells.at["Cluster" + str(cluster_num), roi.decode()] = \
            np.sum(np.logical_and(df_cluster==cluster_num, df_roi==roi))
           
df_cells.to_csv('Ncells_from_Dissection_in_Cluster.csv')
f.close()

df_d_ratio = df_dissection.div(df_dissection.sum(axis=1), axis=0)

n_clusters = df_d_ratio.shape[0]
clusters = np.array([("Cluster" + str(i)) for i in np.arange(n_clusters)])
df_info = pd.DataFrame(0, index=clusters, columns=["Dissections"])
for i in np.arange(n_clusters):
    df_curr = df_d_ratio.iloc[i,:].sort_values(ascending=False)
    info = ""
    for name, val in df_curr.iteritems():
        if val == 0: break
        info = f'{info} {name}: {val:.2%}, '
    df_info.iloc[i,0] = info.strip(" ,")
df_info.to_csv('dissection_info.csv')