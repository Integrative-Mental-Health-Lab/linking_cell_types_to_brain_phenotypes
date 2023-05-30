'''
Create the level 1 (supercluster level) Siletti matrix. Expression level is not transformed.
'''
import h5py
import os
import numpy as np
import numexpr as ne # to speed up
from tqdm import tqdm

f_l2 = h5py.File('adult_human_20221007.agg.loom', 'r')
L2_Ncells = f_l2['col_attrs/NCells']
L2_matrix = f_l2['matrix']

df_supercluster = np.array(pd.read_table("supercluster"))[:,0]
super_list = np.unique(df_supercluster)

f_write = h5py.File('Siletti_L1-supercluster_matrix.h5', 'w')
avg_matrix = f_write.create_dataset("matrix", (L2_matrix.shape[0], super_list.shape[0]), '<f8')
print(avg_matrix.shape)

for i_sc, sc in tqdm(enumerate(super_list)):
    Ncells_sc = L2_Ncells[np.where(df_supercluster == sc)[0]]
    avg_matrix[:, i_sc] = np.sum(L2_matrix[:,np.where(df_supercluster == sc)[0]] * Ncells_sc, axis=1) / sum(Ncells_sc)

f_l2.close()
f_write.close()
