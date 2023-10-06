'''
Create the level 2 (cluster level) Siletti matrix. Cellular expression value is ln(1+p) transformed prior to aggregation.
'''
import h5py
import os
import numpy as np
import numexpr as ne # to speed up
from tqdm import tqdm

level = 'Cluster'
level_name = 'col_attrs/Clusters'
new_file_name = 'Siletti_L2-cluster-log1p_matrix.h5'

f = h5py.File('adult_human_20221007.loom', 'r')
df_cluster = f[level_name][:]
dset = f['matrix']
print(df_cluster.shape)
print(dset.shape)

n_genes = dset.shape[0]
n_clusters = np.max(df_cluster) + 1

#os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

f_write = h5py.File(new_file_name, 'w')
avg_matrix = f_write.create_dataset("matrix", (n_genes, n_clusters), '<f8')
#avg_matrix = f_write['matrix']
print(avg_matrix.shape)

for cluster_num in tqdm(np.arange(n_clusters)):
    dset_sub = dset[:, np.where(df_cluster==cluster_num)[0]]
    avg_matrix[:,cluster_num] = ne.evaluate('sum(log1p(dset_sub), axis=1)')/dset_sub.shape[1]

# Create labels
clusters = np.array([(level + str(i)).encode('UTF-8') for i in np.arange(n_clusters)])
accessions = f['row_attrs/Accession']
c_f_write = f_write.create_dataset(level, data=clusters)
a_f_write = f_write.create_dataset("Accession", data=accessions)

f.close()
f_write.close()
