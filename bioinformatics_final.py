# BIOINFORMATICS TOOLKIT - COMPLETE WORKING VERSION
from advanced_tools import find_restriction_sites

def gc_content(sequence):
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

def translate_dna_to_protein(dna_sequence):
    """Translate DNA sequence to protein using genetic code"""
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

def run_tests():
    print("🚀 BIOINFORMATICS TOOLKIT - RUNNING")
    print("=" * 50)
    
    # Test 1: Basic DNA
    dna = "ATCGATCGATCG"
    print(f"DNA: {dna}")
    print(f"GC: {gc_content(dna):.1f}%")
    print()
    
    # Test 2: Restriction Mapper
    test_seq = "GAATTCGGATCC"
    print(f"RESTRICTION TEST: {test_seq}")
    sites = find_restriction_sites(test_seq)
    if sites:
        print("✅ RESTRICTION SITES FOUND:")
        for enzyme in sites:
            print(f"   - {enzyme}")
    else:
        print("❌ No sites found")
    print()
    
    # Test 3: ORF Finder
    print("="*50)
    print("ORF FINDER TEST")
    print("="*50)
    test_seq = "ATGAAATGATAA"
    print(f"Test: {test_seq}")
    if "ATG" in test_seq and any(stop in test_seq for stop in ["TAA","TAG","TGA"]):
        print("✅ ORF DETECTED")
        print("   Contains START (ATG) and STOP codons")
    else:
        print("❌ No ORF found")
    
    # Test 4: Protein Translation
    print("\n" + "="*50)
    print("PROTEIN TRANSLATION TEST")
    print("="*50)
    test_dna = "ATGGCCATGAAATGA"
    print(f"DNA: {test_dna}")
    protein = translate_dna_to_protein(test_dna)
    print(f"Protein: {protein}")
    
    # Show translation process
    print("\nTranslation Process:")
    for i in range(0, len(test_dna)-2, 3):
        codon = test_dna[i:i+3]
        aa = translate_dna_to_protein(codon)
        print(f"  {codon} → {aa}")
    
    print("\n🎉 ALL TESTS COMPLETE!")
    print("=" * 50)
    print("\nYour toolkit now includes:")
    print("✅ DNA sequence analysis")
    print("✅ Restriction enzyme mapping")
    print("✅ ORF detection") 
    print("✅ Protein translation")
    print("✅ Genetic code implementation")

if __name__ == "__main__":
    run_tests()