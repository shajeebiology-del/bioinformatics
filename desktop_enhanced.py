import sys
from PyQt5.QtWidgets import (QApplication, QMainWindow, QTextEdit, 
                             QVBoxLayout, QWidget, QPushButton, QLabel,
                             QMessageBox)

class EnhancedBioApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle('Bioinformatics Toolkit - Enhanced Version')
        self.setGeometry(100, 100, 900, 700)
        
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)
        
        layout.addWidget(QLabel("Enter DNA Sequence:"))
        self.dna_input = QTextEdit()
        # Insulin gene - WILL find restriction sites
        self.dna_input.setPlainText("""ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG""")
        layout.addWidget(self.dna_input)
        
        # Test buttons
        self.gc_btn = QPushButton("1. Analyze GC Content")
        self.gc_btn.clicked.connect(self.simple_gc)
        layout.addWidget(self.gc_btn)
        
        self.restriction_btn = QPushButton("2. Find Restriction Sites (Enhanced)")
        self.restriction_btn.clicked.connect(self.enhanced_restriction_test)
        layout.addWidget(self.restriction_btn)
        
        self.orf_btn = QPushButton("3. Detect ORFs")
        self.orf_btn.clicked.connect(self.simple_orf)
        layout.addWidget(self.orf_btn)
        
        # Results area
        self.result_display = QTextEdit()
        self.result_display.setReadOnly(True)
        layout.addWidget(self.result_display)
    
    def simple_gc(self):
        try:
            sequence = self.dna_input.toPlainText().strip()
            if not sequence:
                self.result_display.setText("Please enter a DNA sequence")
                return
                
            gc_count = sequence.upper().count('G') + sequence.upper().count('C')
            gc_percent = (gc_count / len(sequence)) * 100
            
            result = f"🧬 GC CONTENT ANALYSIS\n"
            result += f"Sequence length: {len(sequence)} bp\n"
            result += f"GC Content: {gc_percent:.1f}%\n"
            result += f"G+C bases: {gc_count}\n"
            result += f"A+T bases: {len(sequence) - gc_count}\n"
            
            # GC content interpretation
            if gc_percent < 40:
                result += "Classification: AT-rich\n"
            elif gc_percent > 60:
                result += "Classification: GC-rich\n"
            else:
                result += "Classification: Balanced\n"
            
            self.result_display.setText(result)
        except Exception as e:
            self.result_display.setText(f"Error in GC analysis: {str(e)}")
    
    def enhanced_restriction_test(self):
        try:
            sequence = self.dna_input.toPlainText().strip().upper()
            if not sequence:
                self.result_display.setText("Please enter a DNA sequence")
                return
            
            # ENHANCED restriction enzyme database (20+ enzymes)
            restriction_enzymes = {
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
                'EcoRV': 'GATATC',
                'SpeI': 'ACTAGT',
                'NheI': 'GCTAGC',
                'AvrII': 'CCTAGG',
                'NdeI': 'CATATG',
                'SphI': 'GCATGC'
            }
            
            result = f"🔪 ENHANCED RESTRICTION SITE ANALYSIS\n"
            result += f"Sequence length: {len(sequence)} bp\n"
            result += f"Testing {len(restriction_enzymes)} enzymes\n\n"
            
            found_sites = {}
            for enzyme, site in restriction_enzymes.items():
                positions = []
                for i in range(len(sequence) - len(site) + 1):
                    if sequence[i:i+len(site)] == site:
                        positions.append(i)
                if positions:
                    found_sites[enzyme] = {
                        'positions': positions,
                        'recognition_site': site
                    }
            
            if found_sites:
                result += f"✅ Found {len(found_sites)} restriction enzymes:\n\n"
                for enzyme, data in found_sites.items():
                    result += f"🔹 {enzyme} ({data['recognition_site']}):\n"
                    result += f"   Sites: {len(data['positions'])}\n"
                    result += f"   Positions: {data['positions']}\n\n"
                    
                    # Show cutting pattern
                    if data['positions']:
                        result += "   Cutting pattern: "
                        pattern = ["."] * len(sequence)
                        for pos in data['positions']:
                            if pos < len(pattern):
                                pattern[pos] = "|"
                        result += ''.join(pattern) + "\n\n"
            else:
                result += "❌ No restriction sites found with these enzymes.\n"
                result += "Try a different sequence or add more enzymes to the database.\n"
            
            self.result_display.setText(result)
            
        except Exception as e:
            self.result_display.setText(f"Error in restriction analysis: {str(e)}")
            QMessageBox.critical(self, "Error", f"Restriction analysis failed: {str(e)}")
    
    def simple_orf(self):
        try:
            sequence = self.dna_input.toPlainText().strip().upper()
            if not sequence:
                self.result_display.setText("Please enter a DNA sequence")
                return
            
            result = f"🧬 ORF ANALYSIS\n"
            result += f"Sequence length: {len(sequence)} bp\n\n"
            
            # Find START codons
            start_positions = []
            for i in range(len(sequence)-2):
                if sequence[i:i+3] == "ATG":
                    start_positions.append(i)
            
            if start_positions:
                result += f"✅ START codons (ATG) found: {len(start_positions)}\n"
                result += f"   Positions: {start_positions}\n\n"
            else:
                result += "❌ No START codons found\n\n"
            
            # Find STOP codons
            stop_codons = {"TAA": [], "TAG": [], "TGA": []}
            for stop in stop_codons.keys():
                positions = []
                for i in range(len(sequence)-2):
                    if sequence[i:i+3] == stop:
                        positions.append(i)
                stop_codons[stop] = positions
            
            total_stops = sum(len(positions) for positions in stop_codons.values())
            if total_stops > 0:
                result += f"✅ STOP codons found: {total_stops}\n"
                for stop, positions in stop_codons.items():
                    if positions:
                        result += f"   {stop}: {len(positions)} at positions {positions}\n"
            else:
                result += "❌ No STOP codons found\n"
            
            # ORF detection
            if start_positions and total_stops > 0:
                result += f"\n💡 This sequence has protein-coding potential!\n"
                result += f"   Contains both START and STOP codons\n"
            elif start_positions:
                result += f"\n⚠️  Has START but no STOP codons\n"
            elif total_stops > 0:
                result += f"\n⚠️  Has STOP but no START codons\n"
            else:
                result += f"\n❌ No protein-coding regions detected\n"
            
            self.result_display.setText(result)
            
        except Exception as e:
            self.result_display.setText(f"Error in ORF analysis: {str(e)}")

def main():
    app = QApplication(sys.argv)
    window = EnhancedBioApp()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()