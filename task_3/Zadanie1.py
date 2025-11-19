import numpy as np 
from Bio.PDB import PDBList
from Bio import PDB
import matplotlib.pyplot as plt
pdb1=PDBList()
fetch_pdb = pdb1.retrieve_pdb_file('4YWO', file_format='pdb')
parser = PDB.PDBParser()
structure = parser.get_structure('4YWO', 'yw/pdb4ywo.ent')
atoms = []
for model in structure:
    for chain in model:
        for res in chain:
            if res.id[0] == ' ':
                atoms.append(res['CA'])
n=len(atoms)
matrix = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        distance = atoms[i] - atoms[j]
        matrix[i,j] = distance
map = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        if matrix[i, j] <= 8:
            map[i, j] = 1
        else:
            map[i,j] = 0

con_x, con_y = np.where(map==1)
plt.scatter(con_x, con_y, s=1, c='red')
plt.title("Mapa kontaktów")
plt.xlabel("Numer reszty x")
plt.ylabel("Numer reszty y")
plt.show()