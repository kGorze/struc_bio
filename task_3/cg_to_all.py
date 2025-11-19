from Bio import PDB
import Bio.PDB
import sys
import copy

purines = ["G", "A"]
pyrimidines = ["C", "U"]
purines_an = ["N9", "C2", "C6"]
pyrimidines_an = ["N1", "C2", "C4"]
backbone = ["P", "C4'"]

nucleotides = ["A", "G", "C", "U"]

pdb = PDB.PDBParser(QUIET=True)

def load_templates():
    templates = {}
    for nuc in nucleotides:
        template = pdb.get_structure(nuc, f"templates/{nuc}.pdb")
        templates[nuc] = list(template.get_residues())[0]
    return templates
 
input_file = sys.argv[1]
output_file = sys.argv[2]
templates= load_templates()
cg_structure = pdb.get_structure("CG", input_file)
rec_structure = PDB.Structure.Structure("Reconstruct")
rec_model = PDB.Model.Model(0)
rec_structure.add(rec_model)
rec_chain = PDB.Chain.Chain("A")
rec_model.add(rec_chain)
sup = Bio.PDB.Superimposer()

for residue_cg in cg_structure.get_residues():
    res_name = residue_cg.get_resname().strip()
    if res_name not in templates:
        continue
    temp_copy = copy.deepcopy(templates[res_name])
    temp_copy.id = residue_cg.id
    if res_name in purines:
        anchors = purines_an + backbone
    else:
        anchors = pyrimidines_an + backbone
    moving_at = []
    fixed_at = []
    for name in anchors:
        if name in residue_cg and name in temp_copy:
            fixed_at.append(residue_cg[name])
            moving_at.append(temp_copy[name])
    sup.set_atoms(fixed_at, moving_at)
    all_at = list(temp_copy.get_atoms())
    sup.apply(all_at)
    rec_chain.add(temp_copy)
    
io = PDB.PDBIO()
io.set_structure(rec_structure)
io.save(output_file)
