def translate_gene_to_protein(sequence):
    # Mapping codons to one-letter amino acid abbreviations
    # based on the provided genetic code and abbreviation tables
    codon_map = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': 'Stop', 'UAG': 'Stop',
        'UGU': 'C', 'UGC': 'C', 'UGA': 'Stop', 'UGG': 'W',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

    # Standardize input to uppercase and convert DNA (T) to RNA (U)
    sequence = sequence.upper().replace('T', 'U')

    protein_sequence = []

    # Iterate through the sequence in triplets
    for i in range(0, len(sequence) - (len(sequence) % 3), 3):
        codon = sequence[i:i + 3]
        amino_acid = codon_map.get(codon, '?')

        if amino_acid == 'Stop':
            break

        protein_sequence.append(amino_acid)

    return "".join(protein_sequence)


# Example usage
if __name__ == "__main__":
    test = "ATGGCCATGGGCGTATAA"
    protein = translate_gene_to_protein(test)
    print(f"DNA Sequence: {test}")
    print(f"Protein Sequence: {protein}")
