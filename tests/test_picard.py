import unittest
from autoseq.tools.picard import *


class TestPicard(unittest.TestCase):
    def test_picard_collect_insert_size_metrics(self):
        test_job = PicardCollectInsertSizeMetrics()
        test_job.input = "test_input"
        test_job.output_metrics = "test_output"
        cmd = test_job.command()
        self.assertIn('test_input', cmd)
        self.assertIn('test_output', cmd)

    def test_picard_collect_gc_bias_metrics(self):
        test_job = PicardCollectGcBiasMetrics()
        test_job.input = "test_input"
        test_job.reference_sequence = "dummy_reference"
        test_job.output_metrics = "test_output_metrics"
        test_job.output_summary = "test_output_summary"
        test_job.stop_after = "test_stop_after_value"
        cmd = test_job.command()
        self.assertIn('test_input', cmd)
        self.assertIn('dummy_reference', cmd)
        self.assertIn('test_output_metrics', cmd)
        self.assertIn('test_output_summary', cmd)
        self.assertIn('test_stop_after_value', cmd)

    def test_picard_collect_oxo_g_metrics(self):
        test_job = PicardCollectOxoGMetrics()
        test_job.input = "test_input"
        test_job.reference_sequence = "dummy_reference"
        test_job.output_metrics = "test_output_metrics"
        cmd = test_job.command()
        self.assertIn('test_input', cmd)
        self.assertIn('dummy_reference', cmd)
        self.assertIn('test_output_metrics', cmd)

    def test_picard_collect_hs_metrics(self):
        test_job = PicardCollectHsMetrics()
        test_job.input = "test_input"
        test_job.reference_sequence = "dummy_reference"
        test_job.target_regions = "target.bed"
        test_job.bait_regions = "bait.bed"
        test_job.output_metrics = "test_output_metrics"
        cmd = test_job.command()
        self.assertIn('test_input', cmd)
        self.assertIn('dummy_reference', cmd)
        self.assertIn('test_output_metrics', cmd)
        self.assertIn('target.bed', cmd)
        self.assertIn('bait.bed', cmd)

    def test_picard_collect_wgs_metrics(self):
        test_job = PicardCollectWgsMetrics()
        test_job.input = "test_input"
        test_job.reference_sequence = "dummy_reference"
        test_job.minimum_mapping_quality = "dummy_min_mapping"
        test_job.minimum_base_quality = "dummy_min_base"
        test_job.coverage_cap = "dummy_coverage_cap"
        test_job.output_metrics = "test_output_metrics"
        cmd = test_job.command()
        self.assertIn('test_input', cmd)
        self.assertIn('dummy_reference', cmd)
        self.assertIn('dummy_min_mapping', cmd)
        self.assertIn('dummy_min_base', cmd)
        self.assertIn('dummy_coverage_cap', cmd)
        self.assertIn('test_output_metrics', cmd)

    def test_picard_create_sequence_dictionary(self):
        test_job = PicardCreateSequenceDictionary()
        test_job.input = "test_input"
        test_job.output_dict = "test_output"
        cmd = test_job.command()
        self.assertIn('test_input', cmd)
        self.assertIn('test_output', cmd)

    def test_picard_bed_to_interval_list(self):
        test_job = PicardBedToIntervalList()
        test_job.input = "test_input"
        test_job.reference_dict = "dummy_reference_dict"
        test_job.output = "test_output"
        cmd = test_job.command()
        self.assertIn('test_input', cmd)
        self.assertIn('dummy_reference_dict', cmd)
        self.assertIn('test_output', cmd)

    def test_picard_merge_sam_files(self):
        test_job = PicardMergeSamFiles(["input.bam"], "output.bam")
        cmd = test_job.command()
        self.assertIn('input.bam', cmd)
        self.assertIn('output.bam', cmd)

    def test_picard_mark_duplicates(self):
        test_job = PicardMarkDuplicates("input.bam", "output.bam", "dummy_output_metrics")
        cmd = test_job.command()
        self.assertIn('input.bam', cmd)
        self.assertIn('output.bam', cmd)
        self.assertIn('dummy_output_metrics', cmd)
