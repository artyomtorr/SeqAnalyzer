import os
import unittest

from seqanalyzer import (AminoAcidSequence, DNASequence, InvalidAlphabetError,
                         NucleicAcidSequence, RNASequence, filter_fastq)


class TestBiologicalSequences(unittest.TestCase):

    def test_nucleic_acid_sequence(self):
        seq = NucleicAcidSequence("ATCG")
        self.assertEqual(len(seq), 4)
        self.assertEqual(seq[0], "A")

    def test_dna_sequence(self):
        seq = DNASequence("ATCG")
        self.assertEqual(seq.complement().seq, "TAGC")
        self.assertEqual(seq.reverse_complement().seq, "CGAT")
        self.assertAlmostEqual(seq.gc_content(), 50.0)
        self.assertEqual(seq.transcribe().seq, "AUCG")

    def test_rna_sequence(self):
        seq = RNASequence("AUCG")
        self.assertEqual(seq.complement().seq, "UAGC")
        self.assertEqual(seq.reverse_complement().seq, "CGAU")
        self.assertAlmostEqual(seq.gc_content(), 50.0)

    def test_amino_acid_sequence(self):
        seq = AminoAcidSequence("ACDEFGHIKLMNPQRSTVWY")
        self.assertAlmostEqual(seq.compute_hydrophobicity(), 40.0)

    def test_invalid_sequence(self):
        with self.assertRaises(InvalidAlphabetError):
            DNASequence("XYZ").check_alphabet()


class TestFilterFastq(unittest.TestCase):

    def test_output_file_created(self):
        input_path = "test_input.fastq"
        with open(input_path, "w") as file:
            file.write("@read1\nACGT\n+\nHHHH\n")
            file.write("@read2\nAAAA\n+\nHHHH\n")
        filter_fastq(input_path)
        output_path = "fastq_filtrator_results/test_input_filtered.fastq"
        self.assertTrue(os.path.exists(output_path))
        os.remove(input_path)
        os.remove(output_path)

    def test_gc_content_filter(self):
        input_path = "test_input.fastq"
        with open(input_path, "w") as file:
            file.write("@read1\nGCGG\n+\nHHHH\n")
            file.write("@read2\nAAGG\n+\nHHHH\n")
        filter_fastq(input_path, gc_bounds=(40, 60))
        output_path = "fastq_filtrator_results/test_input_filtered.fastq"
        with open(output_path, "r") as file:
            lines = file.readlines()
            self.assertEqual(len(lines), 4)  # Only the second read passes GC-content filter
        os.remove(input_path)
        os.remove(output_path)

    def test_length_filter(self):
        input_path = "test_input.fastq"
        with open(input_path, "w") as file:
            file.write("@read1\nACGT\n+\nHHHH\n")
            file.write("@read2\nAAAAAAA\n+\nHHHHHHH\n")
        filter_fastq(input_path, length_bounds=(3, 4))
        output_path = "fastq_filtrator_results/test_input_filtered.fastq"
        with open(output_path, "r") as file:
            lines = file.readlines()
            self.assertEqual(len(lines), 4)  # Only the first read passes length filter
        os.remove(input_path)
        os.remove(output_path)

    def test_quality_filter(self):
        input_path = "test_input.fastq"
        with open(input_path, "w") as file:
            file.write("@read1\nACGT\n+\n!!!!\n")
            file.write("@read2\nAAAA\n+\n!!!!\n")
        filter_fastq(input_path, quality_threshold=30)
        output_path = "fastq_filtrator_results/test_input_filtered.fastq"
        with open(output_path, "r") as file:
            lines = file.readlines()
            self.assertEqual(len(lines), 0)  # No reads pass quality filter
        os.remove(input_path)
        os.remove(output_path)


if __name__ == '__main__':
    unittest.main()