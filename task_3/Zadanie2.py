import numpy as np 
from Bio.PDB import PDBList, calc_dihedral
from Bio import PDB
import matplotlib.pyplot as plt

pdb1=PDBList()
fetch_pdb = pdb1.retrieve_pdb_file('4YWO', file_format='pdb')
parser = PDB.PDBParser()
structure = parser.get_structure('4YWO', 'yw/pdb4ywo.ent')
phi_list = []
psi_list = []

for model in structure:
    for chain in model:
        res=[resid for resid in chain if resid.id[0] == ' ']
        for i in range(len(res)):
            phi = None
            psi = None
            if i>0:
                try:
                    res_p = res[i-1]
                    res_c = res[i]
                    c_p = res_p['C'].get_vector()
                    n_c = res_c['N'].get_vector()
                    ca_c = res_c['CA'].get_vector()
                    c_c = res_c['C'].get_vector()
                    phi = calc_dihedral(c_p, n_c, ca_c, c_c)
                except:
                    phi = None
            if i < len(res) - 1:
                try:
                    res_c = res[i]
                    res_n = res[i+1]
                    n_c = res_c['N'].get_vector()
                    ca_c = res_c['CA'].get_vector()
                    c_c = res_c['C'].get_vector()
                    n_n = res_n['N'].get_vector()
                    psi = calc_dihedral(n_c, ca_c, c_c, n_n)
                except:
                    psi=None
            if phi is not None and psi is not None:
                phi_list.append(np.degrees(phi))
                psi_list.append(np.degrees(psi))

plt.scatter(phi_list, psi_list, s=10, alpha=0.5)
plt.show()