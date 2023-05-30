'''
Create a level 3 (subcluster level) Siletti matrix. Expression value is not transformed.
'''
import h5py
import os
import numpy as np
import numexpr as ne # to speed up
from tqdm import tqdm

level_name = 'col_attrs/Subclusters'
new_file_name = 'Siletti_L3-subcluster_matrix.h5'

f = h5py.File('adult_human_20221007.loom', 'r')
df_cluster = f[level_name][:]
dset = f['matrix']
print(df_cluster.shape)
print(dset.shape)

n_genes = dset.shape[0]
n_clusters = np.max(df_cluster) + 1

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

f_write = h5py.File(new_file_name, 'r+')
#avg_matrix = f_write.create_dataset("matrix", (n_genes, n_clusters), '<f8')
avg_matrix = f_write['matrix']
print(avg_matrix.shape)
num_col_to_avg = [np.sum(df_cluster == i) for i in np.arange(n_clusters)]

for cluster_num in tqdm(np.arange(2990, n_clusters)):
    dset_sub = dset[:, np.where(df_cluster==cluster_num)[0]]
    avg_matrix[:,cluster_num] = ne.evaluate('sum(dset_sub, axis=1)')/dset_sub.shape[1]

f.close()
f_write.close()
