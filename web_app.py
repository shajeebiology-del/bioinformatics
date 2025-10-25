from flask import Flask, render_template, request, jsonify
import sys
import os

# Import your existing bioinformatics modules
from sequence_analyzer import (
    gc_content, 
    reverse_complement, 
    validate_dna_sequence,
    find_orfs,
    translate_dna_to_protein,
    analyze_orfs_with_translation
)

from advanced_tools import find_restriction_sites

app = Flask(__name__)

# Recognition sequences mapping
RECOGNITION_SEQUENCES = {
    'EcoRI': 'GAATTC', 'BamHI': 'GGATCC', 'HindIII': 'AAGCTT',
    'XhoI': 'CTCGAG', 'NotI': 'GCGGCCGC', 'SacI': 'GAGCTC',
    'KpnI': 'GGTACC', 'PstI': 'CTGCAG', 'SmaI': 'CCCGGG',
    'XbaI': 'TCTAGA', 'SalI': 'GTCGAC', 'NcoI': 'CCATGG',
    'AgeI': 'ACCGGT', 'BglII': 'AGATCT', 'NdeI': 'CATATG',
    'SpeI': 'ACTAGT', 'EcoRV': 'GATATC', 'ScaI': 'AGTACT',
    'PvuI': 'CGATCG', 'SphI': 'GCATGC'
}

def calculate_sequence_stats(sequence):
    """Calculate basic statistics for DNA sequence"""
    sequence = sequence.upper().replace('\n', '').replace(' ', '')
    stats = {
        'length': len(sequence),
        'gc_content': round(gc_content(sequence), 2),
        'at_content': round(100 - gc_content(sequence), 2),
        'a_count': sequence.count('A'),
        't_count': sequence.count('T'),
        'g_count': sequence.count('G'),
        'c_count': sequence.count('C'),
        'is_valid': validate_dna_sequence(sequence)
    }
    return stats

def convert_restriction_sites(restriction_dict, sequence):
    """Convert restriction sites from dictionary to list format"""
    restriction_sites = []
    for enzyme, positions in restriction_dict.items():
        for position in positions:
            # Get the recognition sequence
            seq_pattern = RECOGNITION_SEQUENCES.get(enzyme, 'Unknown')
            restriction_sites.append({
                'enzyme': enzyme,
                'sequence': seq_pattern,
                'position': position
            })
    return restriction_sites

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/analyze', methods=['POST'])
def analyze_sequence():
    try:
        data = request.get_json()
        sequence = data.get('sequence', '').upper().strip()
        analysis_type = data.get('analysis_type', 'all')
        
        # Clean the sequence
        sequence = sequence.replace('\n', '').replace(' ', '')
        
        if not sequence:
            return jsonify({
                'success': False,
                'error': 'Please enter a DNA sequence'
            })
        
        if not validate_dna_sequence(sequence):
            return jsonify({
                'success': False,
                'error': 'Invalid DNA sequence. Only A, T, C, G characters allowed.'
            })
        
        results = {}
        stats = calculate_sequence_stats(sequence)
        
        # Basic sequence info always included
        results['stats'] = stats
        results['reverse_complement'] = reverse_complement(sequence)
        
        # Perform requested analyses
        if analysis_type in ['all', 'orfs']:
            orfs = analyze_orfs_with_translation(sequence, min_orf_length=30)
            results['orfs'] = orfs
        
        if analysis_type in ['all', 'translation']:
            results['protein_translation'] = translate_dna_to_protein(sequence)
        
        if analysis_type in ['all', 'restriction']:
            restriction_sites_dict = find_restriction_sites(sequence)
            restriction_sites = convert_restriction_sites(restriction_sites_dict, sequence)
            results['restriction_sites'] = restriction_sites
        
        return jsonify({
            'success': True,
            'results': results
        })
    
    except Exception as e:
        return jsonify({
            'success': False,
            'error': f'Analysis error: {str(e)}'
        })

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)