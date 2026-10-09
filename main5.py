def parse_fasta(fasta_text):
    """
    Parses multi-record FASTA strings into a list of (header, sequence) tuples.
    """
    records = []
    current_header = None
    current_seq_chunks = []

    for line in fasta_text.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header is not None:
                records.append((current_header, "".join(current_seq_chunks)))
                current_seq_chunks = []
            current_header = line[1:]
        else:
            current_seq_chunks.append(line)

    if current_header is not None:
        records.append((current_header, "".join(current_seq_chunks)))

    return records


def extract_valid_proteins(sequence):

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

    # Standardize input by cleaning whitespace and replacing DNA 'T' with RNA 'U'
    clean_seq = sequence.upper().replace('T', 'U').replace(" ", "").replace("\n", "")
    valid_proteins = []

    # Linear scan using an explicit index cursor instead of collecting all AUGs upfront
    i = 0
    seq_len = len(clean_seq)

    while i < seq_len - 2:
        # Search for the start codon AUG
        if clean_seq[i:i + 3] == 'AUG':
            protein = []
            found_stop = False
            curr = i

            # Translate in steps of 3 from this AUG.
            # Internal AUGs are automatically read as 'M' without resetting or branching,
            # directly obeying: "disregard all AUG until you find a valid stop"
            while curr + 3 <= seq_len:
                codon = clean_seq[curr:curr + 3]
                amino_acid = codon_map.get(codon, '?')

                if amino_acid == 'Stop':
                    found_stop = True
                    # Advance the pointer past the stop codon so downstream search continues
                    i = curr + 3
                    break

                protein.append(amino_acid)
                curr += 3

            # Only keep the sequence if an actual stop codon was reached
            if found_stop:
                valid_proteins.append("".join(protein))
            else:
                # If no stop codon was reached before the end of the transcript,
                # advance by 1 nucleotide to check remaining positions
                i += 1
        else:
            i += 1

    return valid_proteins


# --- Testing with the full Influenza Data containing both FASTA records ---
if __name__ == "__main__":
    influenza_data = """
>gi|58576|gb|X52226|Influenza A virus (A/FPV/Rostock/34(H7N1)) gene for neuraminidase, genomic RNA
AGCAAAAGCAGGAGTTCAAAATGAATCCAAATCAGAAAATAATAACCATTGGGTCAATCTGTATGGGGATCGGAATAATCAGCCTAATATTACAAATTGGAAACATAATCTCAATGTGGGTTAGTCATTCAATTCAGACTGAAAATCAAAATCACCATGAAGCATGCAACCCAAGCATTGCTGGACAGGATGCAGCTTCAGTGGCACTAGCAGGCAATTCCTCTCTTTGTCCCATTAGTGGGTGGGCTATATACAGTAAAGACAATGGTATAAGAATTGGCTCCAAAGGAGACGTATTTGTCATAAGAGAGCCATTTATTTCATGCTCTCACTTGGAATGCAGGACCTTTTTTCTGACTCAGGGCGCCTTGTTGAATGACAAACACTCAAATGGAACCGTTAAAGACAGAAGCCCTTATAGAACCTTGATGAGCTGTCCTGTTGGTGAAGCTCCTTCTCCCTACAATTCAAGGTTCGTGTCGGTTGCATGGTCAGCAAGTGCCTGCCATGATGGCATGGGTTGGCTAACAATCGGAATTTCTGGTCCAGATAATGGAGCAGTGGCTGTATTAAAATACAATGGTATAATAACAGACACCATCAAAAGTTGGAAAAATAACATATTGAGAACGCAAGAGTCTGAATGTGCCTGTATAAATGGTTCATGTTTCACTATAATGACCGATGGCCCAAGTAATGGGCAGGCCTCGTACAAAATTTTCAAGATAGAAAAGGGGAAGGTAGTCAAATCAAGTGAATTAAATGCACCTAATTACCACTACGAGGAATGTTCTTGTTATCCTGATGCAGGTGAAGTAATGTGTGTATGCAGGGATAATTGGCATGGTTCGAATCGACCATGGGTGTCTTTCAATAAAAACCTTGATTATCAAATAGGTTACATCTGCAGTGGGGTTTTTGGTGATAATCCACGACCCAATGATGGAACAGGCAGCTGTGGTCCAGTGTCTTCTAATGGAGCATATGGGATAAAGGGGTTCTCATTTAAGTATGGTAATGGTGTTTGGATAGGGAGAACCAAAAGCACTAGTTCCAGAAGCGGGTTTGAGATGATTTGGGATCCCAATGGATGGACAGAAACTGACAGCAGTTTCTCTGTGAAGCAAGATATTGTAGCAATAACTGATTGGTCGGGATATAGCGGGAGTTTCGTCCAACACCCTGAATTAACAGGGTTGGACTGCATGAGGCCTTGTTTCTGGGTTGAACTGATCAGGGGACGGCCCAACCACAATACGATCTGGACTAGTGGGAGCAGCATTTCCTTCTGTGGTGTGAACAGCGATACTGTAGGTTGGTCTTGGCCAGACGGTGCTGAGTTGCCATTCACCATTGACAAGTAGTTTGTTCAAAAAAACTCCTTGTTTCTACT
>gi|59292|gb|X53029|Influenza A virus (A/USSR/90/1977(H1N1)) genes for matrix proteins 1 and 2, genomic RNA
AGCAAAAGCAGGTAGATGTTGAAAGATGAGTCTTCTAACCGAGGTCGAAACGTACGTTCTCTCTATCGTCCCGTCAGGCCCCCTCAAAGCCGAGATCGCACAGAGACTTGAAGATGTCTTTGCTGGGAAGAACACCGATCTTGAGGCTCTCATGGAATGGCTAAAGACAAGACCAATCCTGTCACCTCTGACTAAGGGGATTTTAGGATTTGTGTTCACGCTCACCGTGCCCAGTGAGCGAGGACTGCAGCGTAGACGCTTTGTCCAAAATGCCCTTAATGGGAATGGGGATCCAAATAACATGGACAGAGCAGTTAAACTGTATAGAAAGCTTAAGAGGGAGATAACATTCCATGGGGCCAAAGAAATAGCACTCAGTTATTCTGCTGGTGCACTTGCCAGTTGTATGGGCCTCATATACAACAGGATGGGGGCTGTGACCACCGAAGCGGCATTTGGCCTGATATGCGCAACCTGTGAACAGATTGCTGACTCCCAGCATAGGTCTCATAGGCAAATGGTGACAACAACCAATCCACTAATAAGACATGAGAACAGAATGGTTCTGGCCAGCACTACAGCTAAGGCTATGGAGCAAATGGCTGGATCGAGTGAGCAAGCAGCAGAGGCCATGGAGGTTGCTAGTCAGGCCAGGCAAATGGTGCAGGCAATGAGAGCCATTGGGACTCATCCTAGCTCCAGTGCTGGTCTGAAAAATGATCTTCTTGAAAATTTGCAGGCCTATCAGAAACGAATGGGGGTGCAGATGCAACGATTCAAGTGATCCTCTTGTTGTTGCCGCAAGTATCATTGGGATTTTGCACTTGATATTGTGGATTCTTGATCGTCTTTTTTTCAAATGCATTTATCGTCTCTTTAAACACGGTCTGAAAAGAGGGCCTTCTACGGAAGGAGTACCAGAGTCTATGAGGGAAGAATATCGAAAGGAACAGCAGAATGCTGTGGATGCTGACGATAGTCATTTTGTCAACATAGAGCTAGAGTAAAAAACTACCTTGTTTCTACT
"""

    # Process each FASTA entry independently instead of fusing them
    records = parse_fasta(influenza_data)
    for header, seq in records:
        print(f"Record: {header}")
        proteins = extract_valid_proteins(seq)
        for idx, p in enumerate(proteins, 1):
            print(f"  Protein {idx} (Length {len(p)} aa):\n  {p}\n")