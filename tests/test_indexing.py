import unittest
from autoseq.tools.indexing import *


class TestIndexing(unittest.TestCase):
    def test_bwa_index(self):
        bwa_index = BwaIndex()
        bwa_index.input_fasta = "input.fasta"
        bwa_index.output = "test_output"
        cmd = bwa_index.command()
        self.assertIn('input.fasta', cmd)

    def test_samtools_faidx(self):
        samtools_faidx = SamtoolsFaidx()
        samtools_faidx.input_fasta = "input.fasta"
        samtools_faidx.output = "test_output"
        cmd = samtools_faidx.command()
        self.assertIn('input.fasta', cmd)

    def test_generate_chr_sizes(self):
        samtools_faidx = GenerateChrSizes()
        samtools_faidx.input_fai = "input.fai"
        samtools_faidx.output = "test_output"
        cmd = samtools_faidx.command()
        self.assertIn('input.fai', cmd)
        self.assertIn('test_output', cmd)
