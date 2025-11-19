from Bio import PDB
import sys

purines = ["G", "A"]
pyrimidines = ["C", "U"]
purine_at = ["N9", "C2", "C6"]
pyrimidines_at = ["N1", "C2", "C4"]
backbone = ["P", "C4'"]

class AtomSelector(PDB.Select):
    def accept_atom(self, atom):
        name = atom.get_name().strip()
        res = atom.get_parent()
        res_name = res.get_resname().strip()
        if name in backbone:
            return 1
        if res_name in purines:
            if name in purine_at:
                return 1
        elif res_name in pyrimidines:
            if name in pyrimidines_at:
                return 1
        return 0
if len(sys.argv) < 3:
    print("Podaj plik wejsciowy i wyjsciowy")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]
pdb = PDB.PDBParser(QUIET=True)
structure = pdb.get_structure("RNA", input_file)
io = PDB.PDBIO()
io.set_structure(structure)
io.save(output_file, AtomSelector())