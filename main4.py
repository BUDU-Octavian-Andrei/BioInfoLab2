def translate_longest_orf(sequence):
    # One-letter amino acid mapping
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

    # 1. Clean sequence: remove headers, whitespace, and convert T -> U
    lines = sequence.splitlines()
    clean_seq = ""
    for line in lines:
        if not line.startswith(">"):  # Skip the FASTA header lines
            clean_seq += line.strip().upper().replace('T', 'U')

    # 2. Find ALL possible start codons (AUG)
    starts = [i for i in range(len(clean_seq)) if clean_seq.startswith('AUG', i)]

    if not starts:
        return "No start codon found."

    longest_protein = ""

    # 3. For every AUG, see how far it goes before hitting a stop
    for start_pos in starts:
        current_protein = []
        # Translate in triplets starting from this specific AUG
        for i in range(start_pos, len(clean_seq) - (len(clean_seq) % 3), 3):
            codon = clean_seq[i:i + 3]
            amino_acid = codon_map.get(codon, '?')

            if amino_acid == 'Stop':
                break  # Found a stop, this specific ORF ends here

            current_protein.append(amino_acid)

        # Join the list into a string
        protein_str = "".join(current_protein)

        # 4. "Disregard all AUG until you find a valid stop"
        # We interpret "valid" as the one that produces the longest protein.
        if len(protein_str) > len(longest_protein):
            longest_protein = protein_str

    return longest_protein


# --- Testing with the Influenza Data provided ---
influenza_data = """
>gi|58576|gb|X52226|Influenza A virus (A/FPV/Rostock/34(H7N1)) gene for neuraminidase, genomic RNA
AGCAAAAGCAGGAGTTCAAAATGAATCCAAATCAGAAAATAATAACCATTGGGTCAATCTGTATGGGGAT
CGGAATAATCAGCCTAATATTACAAATTGGAAACATAATCTCAATGTGGGTTAGTCATTCAATTCAGACT
GAAAATCAAAATCACCATGAAGCATGCAACCCAAGCATTGCTGGACAGGATGCAGCTTCAGTGGCACTAG
CAGGCAATTCCTCTCTTTGTCCCATTAGTGGGTGGGCTATATACAGTAAAGACAATGGTATAAGAATTGG
CTCCAAAGGAGACGTATTTGTCATAAGAGAGCCATTTATTTCATGCTCTCACTTGGAATGCAGGACCTTT
TTTCTGACTCAGGGCGCCTTGTTGAATGACAAACACTCAAATGGAACCGTTAAAGACAGAAGCCCTTATA
GAACCTTGATGAGCTGTCCTGTTGGTGAAGCTCCTTCTCCCTACAATTCAAGGTTCGTGTCGGTTGCATG
GTCAGCAAGTGCCTGCCATGATGGCATGGGTTGGCTAACAATCGGAATTTCTGGTCCAGATAATGGAGCA
GTGGCTGTATTAAAATACAATGGTATAATAACAGACACCATCAAAAGTTGGAAAAATAACATATTGAGAA
CGCAAGAGTCTGAATGTGCCTGTATAAATGGTTCATGTTTCACTATAATGACCGATGGCCCAAGTAATGG
GCAGGCCTCGTACAAAATTTTCAAGATAGAAAAGGGGAAGGTAGTCAAATCAAGTGAATTAAATGCACCT
AATTACCACTACGAGGAATGTTCTTGTTATCCTGATGCAGGTGAAGTAATGTGTGTATGCAGGGATAATT
GGCATGGTTCGAATCGACCATGGGTGTCTTTCAATAAAAACCTTGATTATCAAATAGGTTACATCTGCAG
TGGGGTTTTTGGTGATAATCCACGACCCAATGATGGAACAGGCAGCTGTGGTCCAGTGTCTTCTAATGGA
GCATATGGGATAAAGGGGTTCTCATTTAAGTATGGTAATGGTGTTTGGATAGGGAGAACCAAAAGCACTA
GTTCCAGAAGCGGGTTTGAGATGATTTGGGATCCCAATGGATGGACAGAAACTGACAGCAGTTTCTCTGT
GAAGCAAGATATTGTAGCAATAACTGATTGGTCGGGATATAGCGGGAGTTTCGTCCAACACCCTGAATTA
ACAGGGTTGGACTGCATGAGGCCTTGTTTCTGGGTTGAACTGATCAGGGGACGGCCCAACCACAATACGA
TCTGGACTAGTGGGAGCAGCATTCCTTCTGTGGTGTGAACAGCGATACTGTAGGTTGGTCTTGGCCAGA
CGGTGCTGAGTTGCCATTCACCATTGACAAGTAGTTTGTTCAAAAAAACTCCTTGTTTCTACT
"""

result = translate_longest_orf(influenza_data)
print(f"Longest valid protein sequence:\n{result}")
