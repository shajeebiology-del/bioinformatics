from sequence_analyzer import *

print("🚀 BIOINFORMATICS TOOLKIT - RUNNING")
print("=" * 50)

# Test basic DNA
dna = "ATCGATCGATCG"
print(f"DNA: {dna}")
print(f"GC: {gc_content(dna):.1f}%")
print()

# Test ORF finder
test_orf_finder()
print()

# Test restriction mapper
test_seq = "GAATTCGGATCC"
print(f"RESTRICTION TEST: {test_seq}")
sites = find_restriction_sites(test_seq)
if sites:
    print("✅ RESTRICTION SITES FOUND")
else:
    print("❌ No sites found")
print()

# Test real genes
test_real_genes()

print("🎉 ALL TESTS COMPLETE!")