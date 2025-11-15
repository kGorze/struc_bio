import typing

def read_fasta(file_path: str) -> typing.List[str]:
    sequences = []
    current_sequence = ""
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_sequence:
                    sequences.append(current_sequence)
                current_sequence = ""
            else:
                current_sequence += line
        if current_sequence:
            sequences.append(current_sequence)
    return sequences
