from advanced_tools import find_restriction_sites, print_restriction_sites

def gc_content(sequence):
    """Calculate GC content of DNA sequence"""
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

def reverse_complement(sequence):
    """Generate reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))

def validate_dna_sequence(sequence):
    """Validate if sequence contains only DNA bases"""
    valid_bases = {'A', 'T', 'G', 'C'}
    return all(base in valid_bases for base in sequence.upper())

# ===== ORF FINDER FUNCTIONS =====

def get_reading_frames(sequence):
    """
    Generate all 6 reading frames of DNA sequence
    - 3 forward frames
    - 3 reverse complement frames
    """
    sequence = sequence.upper()
    frames = []
    
    # 3 forward frames
    for i in range(3):
        frames.append(sequence[i:])
    
    # Reverse complement for 3 reverse frames
    rev_comp = reverse_complement(sequence)
    for i in range(3):
        frames.append(rev_comp[i:])
    
    return frames

def find_start_codons(sequence):
    """Find all START codon (ATG) positions in a sequence"""
    start_positions = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if codon == "ATG":
            start_positions.append(i)
    return start_positions

def find_stop_codons(sequence):
    """Find all STOP codon (TAA, TAG, TGA) positions in a sequence"""
    stop_positions = []
    stop_codons = ["TAA", "TAG", "TGA"]
    
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if codon in stop_codons:
            stop_positions.append(i)
    
    return stop_positions

def find_orfs(sequence, min_orf_length=100):
    """
    Find all Open Reading Frames in DNA sequence
    Returns: List of ORFs with start, end positions and sequence
    """
    sequence = sequence.upper()
    orfs = []
    
    # Get all 6 reading frames
    frames = get_reading_frames(sequence)
    
    for frame_num, frame in enumerate(frames):
        # Find all START and STOP codons in this frame
        starts = find_start_codons(frame)
        stops = find_stop_codons(frame)
        
        # For each START codon, find the nearest STOP codon
        for start_pos in starts:
            # Find the first STOP codon after this START
            valid_stops = [stop for stop in stops if stop > start_pos]
            if valid_stops:
                stop_pos = min(valid_stops)
                orf_length = stop_pos - start_pos + 3  # +3 to include STOP codon
                orf_sequence = frame[start_pos:stop_pos + 3]
                
                # Check if ORF meets minimum length requirement
                if orf_length >= min_orf_length:
                    orfs.append({
                        'frame': frame_num + 1,
                        'start': start_pos,
                        'end': stop_pos + 3,
                        'length': orf_length,
                        'sequence': orf_sequence
                    })
    
    return orfs

def print_orfs(orfs):
    """Pretty print ORF results"""
    if not orfs:
        print("No ORFs found meeting the criteria.")
        return
    
    print(f"\nFound {len(orfs)} ORF(s):")
    print("-" * 60)
    for i, orf in enumerate(orfs, 1):
        print(f"ORF {i}:")
        print(f"  Frame: {orf['frame']}")
        print(f"  Position: {orf['start']}-{orf['end']}")
        print(f"  Length: {orf['length']} bp")
        print(f"  Sequence: {orf['sequence'][:50]}..." if len(orf['sequence']) > 50 else f"  Sequence: {orf['sequence']}")
        print()

def test_orf_finder():
    """Test the complete ORF finder with real DNA examples"""
    
    # Test DNA with known ORFs
    test_dna = "ATGCCCTAGATGAAATGA"  # Contains ORFs
    
    print("=== Testing Complete ORF Finder ===")
    print(f"Test DNA: {test_dna}")
    
    # Find ORFs with minimum length of 1 (to find small test ORFs)
    orfs = find_orfs(test_dna, min_orf_length=1)
    
    # Print results
    print_orfs(orfs)
    
    # Test with a longer, more realistic sequence
    longer_dna = "ATG" + "GCC" * 50 + "TAA"  # 153 bp ORF
    print(f"\nTesting longer sequence (length: {len(longer_dna)} bp)")
    orfs_long = find_orfs(longer_dna, min_orf_length=50)
    print_orfs(orfs_long)
    
    print("=== ORF Finder Test Complete ===")
def test_real_genes():
    """Test ORF finder with real biological sequences"""
    
    print("\n" + "="*60)
    print("TESTING WITH REAL BIOLOGICAL SEQUENCES")
    print("="*60)
    
    # Real gene sequences (shortened for testing)
    real_genes = {
        "TP53_tumor_suppressor": "ATGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCG",
        
        "INS_insulin": "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG",
        
        "ACTB_beta_actin": "ATGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAGGCCGGCTTCGCGGGCGACGATGCTCCCCGGGCTGTATGCCTCTGGTCGTACCACTGGCATCGTGATGGACTCCGGTGACGGGGTCACCCACACTGTGCCCATCTACGAGGGCTATGCTCTCCCTCACGCCATCCTGCGTCTGGACCTGGCTGGCCGGGACCTGACTGACTACCTCATGAAGATCCTCACCGAGCGCGGCTACAGCTTCACCACCACAGCTGAGAGGGAAATCGTGCGTGACATTAAGGAGAAGCTGTGCTACGTCGCCCTGGACTTCGAGCAGGAGATGGCCACGGCTGCTTCCAGCTCCTCCCTGGAGAAGAGCTACGAGCTGCCTGACGGCCAGGTCATCACCATTGGCAATGAGCGGTTCCGCTGCCCTGAGGCACTCTTCCAGCCTTCCTTCCTGGGCATGGAGTCCTGTGGCATCCACGAAACTACCTTCAACTCCATCATGAAGTGTGACGTGGACATCCGCAAAGACCTGTACGCCAACACAGTGCTGTCTGGTGGTCCCTCCGTCGCCCTCCCCCTCCCTCATCCTCTCGTCGCATGGAGTCCTGTGGCA"
    }
    
    for gene_name, sequence in real_genes.items():
        print(f"\n🔬 Analyzing: {gene_name}")
        print(f"Sequence length: {len(sequence)} bp")
        
        # Find ORFs with realistic minimum length (typically > 300bp for proteins)
        orfs = find_orfs(sequence, min_orf_length=50)
        
        # Print results
        print_orfs(orfs)
        
        # Show GC content for this gene
        print(f"GC Content: {gc_content(sequence):.2f}%")
        print("-" * 50)
def translate_dna_to_protein(dna_sequence):
    """
    Translate DNA sequence to protein sequence using standard genetic code
    """
    genetic_code = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
    }
    
    protein = ""
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3].upper()
        if len(codon) == 3:
            protein += genetic_code.get(codon, '?')
    return protein

def analyze_orfs_with_translation(sequence, min_orf_length=100):
    """
    Enhanced ORF finder that also translates to protein
    """
    orfs = find_orfs(sequence, min_orf_length)
    
    for orf in orfs:
        # Translate the ORF to protein
        protein_sequence = translate_dna_to_protein(orf['sequence'])
        orf['protein_sequence'] = protein_sequence
        orf['protein_length'] = len(protein_sequence)
    
    return orfs

def print_detailed_orfs(orfs):
    """Enhanced printing with protein translation"""
    if not orfs:
        print("No ORFs found meeting the criteria.")
        return
    
    print(f"\nFound {len(orfs)} ORF(s) with Protein Translation:")
    print("=" * 70)
    for i, orf in enumerate(orfs, 1):
        print(f"ORF {i}:")
        print(f"  Frame: {orf['frame']}")
        print(f"  DNA Position: {orf['start']}-{orf['end']}")
        print(f"  DNA Length: {orf['length']} bp")
        print(f"  Protein Length: {orf['protein_length']} aa")
        print(f"  DNA Sequence: {orf['sequence'][:30]}..." if len(orf['sequence']) > 30 else f"  DNA Sequence: {orf['sequence']}")
        print(f"  Protein: {orf['protein_sequence'][:30]}..." if len(orf['protein_sequence']) > 30 else f"  Protein: {orf['protein_sequence']}")
        print()

def test_real_genes():
    """Test ORF finder with real biological sequences and translation"""
    
    print("\n" + "="*60)
    print("TESTING WITH REAL BIOLOGICAL SEQUENCES + PROTEIN TRANSLATION")
    print("="*60)
    
    # Real gene sequences
    real_genes = {
        "TP53_tumor_suppressor": "ATGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCG",
        
        "INS_insulin": "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG",
        
        "ACTB_beta_actin": "ATGGATGATGATATCGCCGCGCTCGTCGTCGACAACGGCTCCGGCATGTGCAAGGCCGGCTTCGCGGGCGACGATGCTCCCCGGGCTGTATGCCTCTGGTCGTACCACTGGCATCGTGATGGACTCCGGTGACGGGGTCACCCACACTGTGCCCATCTACGAGGGCTATGCTCTCCCTCACGCCATCCTGCGTCTGGACCTGGCTGGCCGGGACCTGACTGACTACCTCATGAAGATCCTCACCGAGCGCGGCTACAGCTTCACCACCACAGCTGAGAGGGAAATCGTGCGTGACATTAAGGAGAAGCTGTGCTACGTCGCCCTGGACTTCGAGCAGGAGATGGCCACGGCTGCTTCCAGCTCCTCCCTGGAGAAGAGCTACGAGCTGCCTGACGGCCAGGTCATCACCATTGGCAATGAGCGGTTCCGCTGCCCTGAGGCACTCTTCCAGCCTTCCTTCCTGGGCATGGAGTCCTGTGGCATCCACGAAACTACCTTCAACTCCATCATGAAGTGTGACGTGGACATCCGCAAAGACCTGTACGCCAACACAGTGCTGTCTGGTGGTCCCTCCGTCGCCCTCCCCCTCCCTCATCCTCTCGTCGCATGGAGTCCTGTGGCA"
    }
    
    for gene_name, sequence in real_genes.items():
        print(f"\n🔬 Analyzing: {gene_name}")
        print(f"Sequence length: {len(sequence)} bp")
        print(f"GC Content: {gc_content(sequence):.2f}%")
        
        # Find ORFs with protein translation
        orfs = analyze_orfs_with_translation(sequence, min_orf_length=50)
        
        # Print detailed results
        print_detailed_orfs(orfs)
        print("-" * 50)

def test_integrated_analysis():
    """Test combined ORF finding and restriction mapping"""
    print("\n" + "="*60)
    print("INTEGRATED ANALYSIS: ORF FINDING + RESTRICTION MAPPING")
    print("="*60)
    
    # Test with a gene that has both ORFs and restriction sites
    test_gene = "ATGGAATTCGCCGGATCCATGAAATGA"  # Contains ORF and restriction sites
    
    print(f"Test Sequence: {test_gene}")
    print(f"Length: {len(test_gene)} bp")
    print(f"GC Content: {gc_content(test_gene):.2f}%")
    
    # ORF Analysis
    print("\n--- ORF ANALYSIS ---")
    orfs = find_orfs(test_gene, min_orf_length=1)
    print_orfs(orfs)
    
    # Restriction Analysis
    print("--- RESTRICTION ANALYSIS ---")
    restriction_sites = find_restriction_sites(test_gene)
    print_restriction_sites(restriction_sites, test_gene)
    
    # Combined insights
    if orfs and restriction_sites:
        print("--- MOLECULAR CLONING INSIGHTS ---")
        print("This sequence contains both protein-coding regions and")
        print("restriction sites useful for genetic engineering.")
from advanced_tools import find_restriction_sites, print_restriction_sites

print("🧪 TESTING RESTRICTION MAPPER DIRECTLY")
print("=" * 50)

# Test sequence with known restriction sites
test_seq = "GAATTCGGATCCAAGCTT"  # EcoRI, BamHI, HindIII

print(f"Sequence: {test_seq}")
print(f"Length: {len(test_seq)} bp")

# Call restriction functions directly
sites = find_restriction_sites(test_seq)
print(f"Raw result: {sites}")

if sites:
    print("\n📊 RESTRICTION SITES:")
    print_restriction_sites(sites, test_seq)
else:
    print("❌ No restriction sites found")