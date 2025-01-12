import streamlit as st
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import io

class DNASequenceAnalyzer:
    def __init__(self, sequence):
        self.sequence = Seq(sequence)

    def validate_sequence(self):
        """Check if the sequence contains only valid DNA bases."""
        valid_bases = {"A", "T", "G", "C"}
        return set(self.sequence.upper()).issubset(valid_bases)

    def gc_content(self):
        gc_count = self.sequence.count('G') + self.sequence.count('C')
        return (gc_count / len(self.sequence)) * 100

    def transcribe(self):
        return self.sequence.transcribe()

    def translate(self):
        return self.sequence.translate()
    
    def find_motif(self, motif):
        locations = []
        start = 0
        while True:
            start = self.sequence.find(motif, start)
            if start == -1:
                break
            locations.append(start)
            start += len(motif)
        return locations

    def reverse_complement(self):
        return self.sequence.reverse_complement()

    def align_sequences(self, seq2):
        alignments = pairwise2.align.globalxx(self.sequence, Seq(seq2))
        return alignments

    @staticmethod
    def build_phylogenetic_tree(sequences):
        alignment = MultipleSeqAlignment(
            [SeqRecord(Seq(seq), id=name) for name, seq in sequences.items()]
        )
        calculator = DistanceCalculator("identity")
        distance_matrix = calculator.get_distance(alignment)
        constructor = DistanceTreeConstructor(calculator)
        tree = constructor.nj(distance_matrix)
        return tree

# Streamlit Interface
st.title("DNA Sequence Analyzer and Phylogenetic Tree Builder")

st.sidebar.header("Input Sequences")
sequences = {}

# Input multiple DNA sequences
sequence_count = st.sidebar.number_input("Number of Sequences", min_value=1, step=1, value=2)
for i in range(sequence_count):
    species_name = st.sidebar.text_input(f"Species {i + 1} Name", f"Species_{i + 1}")
    dna_sequence = st.sidebar.text_area(f"DNA Sequence for {species_name}", key=f"sequence_{i + 1}")
    if species_name and dna_sequence:
        sequences[species_name] = dna_sequence

# DNA Analysis for Individual Sequences
st.header("Individual Sequence Analysis")
selected_species = st.selectbox("Select a Species for Analysis", list(sequences.keys()))
if selected_species:
    dna_analyzer = DNASequenceAnalyzer(sequences[selected_species])

    if not dna_analyzer.validate_sequence():
        st.error("Invalid DNA sequence. Please ensure the sequence contains only A, T, G, and C.")
    else:
        st.subheader("GC Content")
        st.write(f"GC Content: {dna_analyzer.gc_content():.2f}%")

        st.subheader("Transcription")
        st.write(dna_analyzer.transcribe())

        st.subheader("Translation")
        st.write(dna_analyzer.translate())

        st.subheader("Find Motif")
        motif = st.text_input("Enter Motif to Search", "AGC")
        if motif:
            locations = dna_analyzer.find_motif(motif)
            st.write(f"Motif '{motif}' found at positions: {locations}")
        
        # Export Results
        if st.button("Export Results"):
            result = f"""Species: {species_name}
GC Content: {dna_analyzer.gc_content():.2f}%
Transcription: {dna_analyzer.transcribe()}
Translation: {dna_analyzer.translate()}
Reverse Complement: {dna_analyzer.reverse_complement()}"""
            st.download_button("Download Results", result, file_name="sequence_analysis.txt")

# Reverse Complement Analysis
st.header("Reverse Complement Analysis")
if selected_species:
    reverse_complement = dna_analyzer.reverse_complement()
    reverse_analyzer = DNASequenceAnalyzer(str(reverse_complement))

    st.subheader("Reverse Complement")
    st.write(reverse_complement)

    st.subheader("Transcription")
    st.write(reverse_analyzer.transcribe())

    st.subheader("Translation")
    st.write(reverse_analyzer.translate())

    st.subheader("Find Motif")
    motif_rc = st.text_input("Enter Motif to Search in Reverse Complement", "AGC")
    if motif_rc:
        locations_rc = reverse_analyzer.find_motif(motif_rc)
        st.write(f"Motif '{motif_rc}' found at positions: {locations_rc}")

    # Export Reverse Complement Results
    if st.button("Export Reverse Complement Results"):
        result_rc = f"""Species: {species_name}
Reverse Complement: {reverse_complement}
GC Content (Reverse Complement): {reverse_analyzer.gc_content():.2f}%
Transcription (Reverse Complement): {reverse_analyzer.transcribe()}
Translation (Reverse Complement): {reverse_analyzer.translate()}
Motif '{motif_rc}' found at positions (Reverse Complement): {locations_rc if motif_rc else 'Not searched'}"""
        st.download_button("Download Reverse Complement Results", result_rc, file_name="reverse_complement_analysis.txt")

# Pairwise Sequence Alignment
st.header("Pairwise Sequence Alignment")
if len(sequences) > 1:
    seq1 = st.selectbox("Select Sequence 1", list(sequences.keys()), key="align_seq1")
    seq2 = st.selectbox("Select Sequence 2", list(sequences.keys()), key="align_seq2")
    if seq1 != seq2:
        analyzer1 = DNASequenceAnalyzer(sequences[seq1])
        if analyzer1.validate_sequence():
            alignments = analyzer1.align_sequences(sequences[seq2])
            st.subheader("Alignment")
            st.write(alignments[0])

# Build Phylogenetic Tree
st.header("Phylogenetic Tree")
if st.button("Build Phylogenetic Tree"):
    if len(sequences) > 1:
        valid_sequences = all(DNASequenceAnalyzer(seq).validate_sequence() for seq in sequences.values())
        if valid_sequences:
            tree = DNASequenceAnalyzer.build_phylogenetic_tree(sequences)
            fig = plt.figure(figsize=(10, 6))
            ax = fig.add_subplot(1, 1, 1)
            Phylo.draw(tree, do_show=False, axes=ax)
            st.pyplot(fig)

            # Save Tree
            with io.StringIO() as output:
                Phylo.write(tree, output, "newick")
                tree_newick = output.getvalue()
            st.download_button("Download Tree", tree_newick, file_name="phylogenetic_tree.newick")
        else:
            st.warning("Invalid sequences detected. Ensure all sequences contain only A, T, G, and C.")
    else:
        st.warning("Please input at least two sequences to build a phylogenetic tree.")

# Load Tree
st.subheader("Load and Display Saved Tree")
uploaded_tree = st.file_uploader("Upload a Newick Tree File", type=["newick"])
if uploaded_tree:
    try:
        # Decode the uploaded file to a string
        tree_data = uploaded_tree.getvalue().decode("utf-8")
        
        # Read the tree using Phylo.read from the string
        tree = Phylo.read(io.StringIO(tree_data), "newick")
        
        # Plot the tree
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, do_show=False, axes=ax)
        st.pyplot(fig)

    except Exception as e:
        st.error(f"An error occurred while loading the tree: {e}")