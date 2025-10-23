"""
Sequence Alignment Tools
- Local alignment (BLAST-like)
- Global alignment (Needleman-Wunsch)
- Similarity analysis
"""

def simple_alignment(seq1, seq2, match_score=2, mismatch_penalty=-1, gap_penalty=-1):
    """
    Simple local sequence alignment using basic scoring
    """
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    
    # Find the best local alignment
    best_score = 0
    best_position = (0, 0)
    best_length = 0
    
    # Try all possible starting positions
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            score = 0
            length = 0
            
            # Extend alignment while score is positive
            k = 0
            while (i + k < len(seq1) and j + k < len(seq2)):
                if seq1[i + k] == seq2[j + k]:
                    score += match_score
                else:
                    score += mismatch_penalty
                
                # Stop if score goes too negative
                if score < 0:
                    break
                    
                length = k + 1
                k += 1
                
                # Update best alignment
                if score > best_score:
                    best_score = score
                    best_position = (i, j)
                    best_length = length
    
    # Extract the aligned regions
    if best_length > 0:
        i, j = best_position
        aligned_seq1 = seq1[i:i + best_length]
        aligned_seq2 = seq2[j:j + best_length]
        
        return {
            'score': best_score,
            'position1': i,
            'position2': j,
            'length': best_length,
            'aligned_seq1': aligned_seq1,
            'aligned_seq2': aligned_seq2,
            'identity': calculate_identity(aligned_seq1, aligned_seq2)
        }
    else:
        return None

def calculate_identity(seq1, seq2):
    """Calculate percentage identity between two sequences"""
    if len(seq1) != len(seq2) or len(seq1) == 0:
        return 0
    
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return (matches / len(seq1)) * 100

def print_alignment(alignment, seq1_name="Sequence 1", seq2_name="Sequence 2"):
    """Pretty print alignment results"""
    if not alignment:
        print("No significant alignment found.")
        return
    
    print(f"\n🔍 SEQUENCE ALIGNMENT RESULTS:")
    print("=" * 50)
    print(f"Alignment Score: {alignment['score']}")
    print(f"Alignment Length: {alignment['length']} bp")
    print(f"Percent Identity: {alignment['identity']:.1f}%")
    print(f"Position in {seq1_name}: {alignment['position1']}-{alignment['position1'] + alignment['length']}")
    print(f"Position in {seq2_name}: {alignment['position2']}-{alignment['position2'] + alignment['length']}")
    
    print(f"\nAligned Regions:")
    print(f"{seq1_name}: {alignment['aligned_seq1']}")
    print(f"{seq2_name}: {alignment['aligned_seq2']}")
    
    # Show match/mismatch visualization
    print("\nMatch Visualization:")
    match_line = ""
    for a, b in zip(alignment['aligned_seq1'], alignment['aligned_seq2']):
        if a == b:
            match_line += "|"
        else:
            match_line += " "
    print(f"                 {match_line}")

def test_alignment():
    """Test the sequence alignment functions"""
    print("🧬 TESTING SEQUENCE ALIGNMENT")
    print("=" * 40)
    
    # Test cases
    test_cases = [
        ("ATCGATCG", "ATCGATCG", "Identical sequences"),
        ("ATCGATCG", "ATCGCTCG", "Single mismatch"),
        ("ATCGATCG", "GCTAGCTA", "Different sequences"),
        ("ATGCCCTAG", "ATGAAATGA", "Realistic gene fragments")
    ]
    
    for seq1, seq2, description in test_cases:
        print(f"\n📋 Test: {description}")
        print(f"Seq1: {seq1}")
        print(f"Seq2: {seq2}")
        
        alignment = simple_alignment(seq1, seq2)
        print_alignment(alignment)

if __name__ == "__main__":
    test_alignment()