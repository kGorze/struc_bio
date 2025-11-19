from Bio import PDB
import sys

purines = ["G", "A"]
pyrimidines = ["C", "U"]
purine_at = ["N9", "C2", "C6"]
pyrimidines_at = ["N1", "C2", "C4"]
backbone = ["P", "C4'"]

class AtomSelector(PDB.Select):
    def verify(self, atom):
        name = atom.get_name()
        res = atom.get_parent()
        res_name = res.get_resname().strip()
        if name in backbone:
            return True
        if res_name in purines:
            if name in purine_at:
                return True
        elif res_name in pyrimidines:
            if name in pyrimidines_at:
                return True
        return False
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