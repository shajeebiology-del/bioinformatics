"""
Advanced Bioinformatics Tools
- Restriction enzyme mapping
- Sequence alignment
- Phylogenetic analysis
"""

# Restriction enzyme database (common enzymes with recognition sequences)
RESTRICTION_ENZYMES = {
    'EcoRI': 'GAATTC',
    'BamHI': 'GGATCC', 
    'HindIII': 'AAGCTT',
    'XbaI': 'TCTAGA',
    'NotI': 'GCGGCCGC',
    'SacI': 'GAGCTC',
    'KpnI': 'GGTACC',
    'PstI': 'CTGCAG',
    'SmaI': 'CCCGGG',
    'XhoI': 'CTCGAG',
    'SalI': 'GTCGAC',
    'NcoI': 'CCATGG',
    'AgeI': 'ACCGGT',
    'BglII': 'AGATCT',
    'EcoRV': 'GATATC'
}

def find_restriction_sites(sequence, enzyme_name=None):
    """
    Find restriction enzyme cutting sites in DNA sequence
    
    Args:
        sequence: DNA sequence to analyze
        enzyme_name: Specific enzyme to search for, or None for all enzymes
    
    Returns:
        Dictionary of enzyme: [cut_positions]
    """
    sequence = sequence.upper()
    results = {}
    
    # Determine which enzymes to search for
    enzymes_to_search = [enzyme_name] if enzyme_name else RESTRICTION_ENZYMES.keys()
    
    for enzyme in enzymes_to_search:
        if enzyme in RESTRICTION_ENZYMES:
            recognition_seq = RESTRICTION_ENZYMES[enzyme]
            cut_positions = []
            
            # Find all occurrences of recognition sequence
            for i in range(len(sequence) - len(recognition_seq) + 1):
                if sequence[i:i+len(recognition_seq)] == recognition_seq:
                    cut_positions.append(i)
            
            if cut_positions:
                results[enzyme] = cut_positions
    
    return results

def print_restriction_sites(results, sequence):
    """Pretty print restriction site results"""
    if not results:
        print("No restriction sites found.")
        return
    
    print(f"\nFound restriction sites in {len(sequence)} bp sequence:")
    print("=" * 60)
    
    for enzyme, positions in results.items():
        recognition_seq = RESTRICTION_ENZYMES[enzyme]
        print(f"{enzyme}: {recognition_seq}")
        print(f"  Cutting positions: {positions}")
        print(f"  Number of sites: {len(positions)}")
        
        # Show cutting pattern
        if positions:
            print("  Cutting pattern: ", end="")
            pattern = ["."] * len(sequence)
            for pos in positions:
                # Mark the cut site (usually after the first base in palindromic sites)
                cut_pos = pos + 1  # For most enzymes, cut after first base
                if cut_pos < len(pattern):
                    pattern[cut_pos] = "|"
            print(''.join(pattern))
        print()

def test_restriction_mapper():
    """Test the restriction enzyme mapper"""
    print("🧬 TESTING RESTRICTION ENZYME MAPPER")
    print("=" * 50)
    
    # Test sequences
    test_sequences = [
        "GAATTCGGATCCAAGCTTTCTAGAGCGGCCGC",  # Contains multiple sites
        "ATCGATCGATCGATCGATCG",  # No restriction sites
        "GAATTC" * 3  # Multiple EcoRI sites
    ]
    
    for i, seq in enumerate(test_sequences, 1):
        print(f"\nTest {i}: {seq}")
        sites = find_restriction_sites(seq)
        print_restriction_sites(sites, seq)

if __name__ == "__main__":
    test_restriction_mapper()