import unittest
from autoseq.tools.contamination import *


class TestContamination(unittest.TestCase):
    def test_create_contest_vcfs(self):
        create_contest_vcfs = CreateContestVCFs()
        create_contest_vcfs.input_target_regions_bed_1 = "test1.bed"
        create_contest_vcfs.input_target_regions_bed_2 = "test2.bed"
        create_contest_vcfs.input_population_vcf = "test.vcf"
        create_contest_vcfs.output = "output.txt"
        cmd = create_contest_vcfs.command()
        self.assertIn('test1.bed', cmd)
        self.assertIn('test2.bed', cmd)
        self.assertIn('test.vcf', cmd)
        self.assertIn('output.txt', cmd)

    def test_contest(self):
        contest = ContEst()
        contest.reference_genome = "dummy.fasta"
        contest.input_eval_bam = "test_eval.bam"
        contest.input_genotype_bam = "test_genotype.bam"
        contest.input_population_af_vcf = "test.vcf"
        contest.output = "output.txt"
        cmd = contest.command()
        self.assertIn('dummy.fasta', cmd)
        self.assertIn('test_eval.bam', cmd)
        self.assertIn('test_genotype.bam', cmd)
        self.assertIn('test.vcf', cmd)
        self.assertIn('output.txt', cmd)

    def test_contest_to_contam_caveat(self):
        contest_to_contam_caveat = ContEstToContamCaveat()
        contest_to_contam_caveat.input_contest_results = "input.txt"
        contest_to_contam_caveat.output = "output.txt"
        cmd = contest_to_contam_caveat.command()
        self.assertIn('input.txt', cmd)
        self.assertIn('output.txt', cmd)
