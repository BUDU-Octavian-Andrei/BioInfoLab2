def translate_gene_to_protein(sequence):
    # Mapping codons to one-letter amino acid abbreviations
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

    # 1. Standardize input
    sequence = sequence.upper().replace('T', 'U')

    # 2. FIND THE START: Search for the first 'AUG'
    start_index = sequence.find('AUG')

    # If no AUG is found, there is no valid coding sequence
    if start_index == -1:
        return "No valid start codon (AUG) found."

    # Slice the sequence to start from the AUG
    coding_sequence = sequence[start_index:]

    protein_sequence = []

    # 3. TRANSLATE: Process in triplets starting from the AUG
    for i in range(0, len(coding_sequence) - (len(coding_sequence) % 3), 3):
        codon = coding_sequence[i:i + 3]
        amino_acid = codon_map.get(codon, '?')

        # STOP RULE: Stop immediately when a Stop codon is reached
        if amino_acid == 'Stop':
            break

        protein_sequence.append(amino_acid)

    # If the loop finished without hitting a 'Stop', the sequence is technically incomplete
    # but we return what we have.
    return "".join(protein_sequence)


# Example usage
if __name__ == "__main__":
    # This sequence has "junk" at the start (GGGC), then the start codon (AUG),
    # then a protein sequence, then a stop codon (UAA), then more junk (CCC).
    test_dna = "GGGCATGGCCATGGGCGTATAACCC"

    protein = translate_gene_to_protein(test_dna)
    print(f"Input Sequence: {test_dna}")
    print(f"Extracted Protein: {protein}")
    # Expected result: "MAMGV" (AUG=M, GCC=A, AUG=M, GGC=G, GUA=V)
