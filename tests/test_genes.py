import unittest
from autoseq.tools.genes import *


class TestGenes(unittest.TestCase):
    def test_filter_gtf_chromosomes(self):
        filter_gtf_chromosomes = FilterGTFChromosomes()
        filter_gtf_chromosomes.input = "test_input"
        filter_gtf_chromosomes.output = "test_output"
        cmd = filter_gtf_chromosomes.command()
        self.assertIn('test_input', cmd)
        self.assertIn('test_output', cmd)

    def test_filter_gtf_genes(self):
        filter_gtf_genes = FilterGTFGenes()
        filter_gtf_genes.input = "test_input"
        filter_gtf_genes.output = "test_output"
        cmd = filter_gtf_genes.command()
        self.assertIn('test_input', cmd)
        self.assertIn('test_output', cmd)

    def test_gtf_2_gene_pred(self):
        filter_gtf_2_gene_pred = GTF2GenePred()
        filter_gtf_2_gene_pred.input = "test_input"
        filter_gtf_2_gene_pred.output = "test_output"
        cmd = filter_gtf_2_gene_pred.command()
        self.assertIn('test_input', cmd)
        self.assertIn('test_output', cmd)
