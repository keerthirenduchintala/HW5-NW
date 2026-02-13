# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # Create NeedlemanWunsch object
    nw = NeedlemanWunsch(
        sub_matrix_file="./substitution_matrices/BLOSUM62.mat", 
        gap_open=-10,
        gap_extend=-1
    )

    # Align each species to human
    gg_score, a, b  = nw.align(hs_seq, gg_seq)
    mm_score, a, b = nw.align(hs_seq, mm_seq)
    br_score, a, b = nw.align(hs_seq, br_seq)
    tt_score, a, b = nw.align(hs_seq, tt_seq)

    # Store results
    results = [
        ("Gallus_gallus", gg_score),
        ("Mus_musculus", mm_score),
        ("Balaeniceps_rex", br_score),
        ("tursiops_truncatus", tt_score)
    ]

    # Sort by score (highest = most similar)
    results_sorted = sorted(results, key=lambda x: x[1], reverse=True)

    # Print results
    print("Species in order of similarity to human BRD2:")
    for species, score in results_sorted:
        print(f"{species}: {score}")

if __name__ == "__main__":
    main()