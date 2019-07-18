import unittest
from autoseq.tools.qc import *


class TestQC(unittest.TestCase):
    def test_heterozygote_concordance(self):
        test_job = HeterzygoteConcordance()
        test_job.input_vcf = "input.vcf"
        test_job.input_bam = "input.bam"
        test_job.reference_sequence = "dummy_reference"
        test_job.target_regions = "target.bed"
        test_job.normalid = "normalid"
        test_job.output = "test_output"
        cmd = test_job.command()
        self.assertIn('input.vcf', cmd)
        self.assertIn('input.bam', cmd)
        self.assertIn('dummy_reference', cmd)
        self.assertIn('target.bed', cmd)
        self.assertIn('normalid', cmd)
        self.assertIn('test_output', cmd)

    def test_fast_qc(self):
        test_job = FastQC()
        test_job.input = "test_input"
        test_job.output = "test_output"
        test_job.outdir = "test_outdir"
        cmd = test_job.command()
        self.assertIn('test_input', cmd)
        self.assertIn('test_outdir', cmd)

    def test_multi_qc(self):
        test_job = MultiQC()
        test_job.input_files = ["test_input1", "test_input2"]
        test_job.search_dir = "dummy_search_dir"
        test_job.output = "test_output"
        test_job.report_title = "dummy_report_title"
        cmd = test_job.command()
        self.assertIn('dummy_search_dir', cmd)
        self.assertIn('dummy_report_title', cmd)

    def test_sambamba_depth(self):
        test_job = SambambaDepth()
        test_job.input = "test_input"
        test_job.targets_bed = "targets.bed"
        test_job.output = "test_output"
        cmd = test_job.command()
        self.assertIn('test_input', cmd)
        self.assertIn('targets.bed', cmd)
        self.assertIn('test_output', cmd)

    def test_bedtools_coverage_histogram(self):
        test_job = BedtoolsCoverageHistogram()
        test_job.input_bam = "input.bam"
        test_job.input_bed = "input.bed"
        test_job.tag = "test_tag"
        test_job.output = "test_output"
        cmd = test_job.command()
        self.assertIn('input.bam', cmd)
        self.assertIn('input.bed', cmd)
        self.assertIn('test_tag', cmd)
        self.assertIn('test_output', cmd)

    def test_coverage_histogram(self):
        test_job = CoverageHistogram()
        test_job.input_bam = "input.bam"
        test_job.input_bed = "input.bed"
        test_job.output = "test_output"
        cmd = test_job.command()
        self.assertIn('input.bam', cmd)
        self.assertIn('input.bed', cmd)
        self.assertIn('test_output', cmd)

    def test_coverage_caveat(self):
        test_job = CoverageCaveat()
        test_job.input_histogram = "dummy_input"
        test_job.output = "test_output"
        cmd = test_job.command()
        self.assertIn('dummy_input', cmd)
        self.assertIn('test_output', cmd)
