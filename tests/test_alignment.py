import unittest

from mock import patch

from autoseq.tools.alignment import *
from autoseq.pipeline.clinseq import ClinseqPipeline


class TestAlignment(unittest.TestCase):
    def setUp(self):
        sample_data = {
            "sdid": "P-NA12877",
            "T": ["AL-P-NA12877-T-03098849-TD1-TT1", "AL-P-NA12877-T-03098849-TD1-WGS"],
            "N": ["AL-P-NA12877-N-03098121-TD1-TT1", "AL-P-NA12877-N-03098121-TD1-WGS"],
            "CFDNA": ["LB-P-NA12877-CFDNA-03098850-TD1-TT1", "LB-P-NA12877-CFDNA-03098850-TD1-TT2",
                      "LB-P-NA12877-CFDNA-03098850-TD1-WGS"]
        }
        ref_data = {
            "bwaIndex": "bwa/test-genome-masked.fasta",
            "chrsizes": "genome/test-genome-masked.chrsizes.txt",
            "clinvar": "variants/clinvar_20160203.vcf.gz",
            "cosmic": "variants/CosmicCodingMuts_v71.vcf.gz",
            "dbSNP": "variants/dbsnp142-germline-only.vcf.gz",
            "exac": "variants/ExAC.r0.3.1.sites.vep.vcf.gz",
            "icgc": "variants/icgc_release_20_simple_somatic_mutation.aggregated.vcf.gz",
            "reference_dict": "genome/test-genome-masked.dict",
            "reference_genome": "genome/test-genome-masked.fasta",
            "swegene_common": "variants/swegen_common.vcf.gz",
            "targets": {
                "test-regions": {
                    "cnvkit-ref": None,
                    "msisites": "intervals/targets/test-regions.msisites.tsv",
                    "targets-bed-slopped20": "intervals/targets/test-regions-GRCh37.slopped20.bed",
                    "targets-interval_list": "intervals/targets/test-regions-GRCh37.slopped20.interval_list",
                    "targets-interval_list-slopped20": "intervals/targets/test-regions-GRCh37.slopped20.interval_list"
                }
            },
            "contest_vcfs": {
                "test-regions": "test_contest.vcf"
            },
            "vep_dir": None
        }
        self.test_clinseq_pipeline = ClinseqPipeline(sample_data, ref_data, {}, "/tmp", "/nfs/LIQBIO/INBOX/exomes")

    def test_bwa_pe(self):
        """
        test bwa paired end, should contain string "foo_1.fq  foo_2.fq" with a double-space between the file names
        """
        bwa = Bwa()
        bwa.input_fastq1 = "foo_1.fq"
        bwa.input_fastq2 = "foo_2.fq"
        bwa.input_reference_sequence = "ref.fasta"
        bwa.readgroup = "__readgroup__"
        bwa.output = "out.bam"
        cmd = bwa.command()
        self.assertIn('foo_1.fq  foo_2.fq', cmd)

    def test_bwa_se(self):
        """
        test bwa single end, should not include any reference to foo_2.fq
        """
        bwa = Bwa()
        bwa.input_fastq1 = "foo_1.fq"
        bwa.input_fastq2 = None
        bwa.input_reference_sequence = "ref.fasta"
        bwa.readgroup = "__readgroup__"
        bwa.output = "out.bam"
        cmd = bwa.command()
        self.assertIn('foo_1.fq', cmd)
        self.assertNotIn('foo_2.fq', cmd)

    def test_bwa_removes_temporary_logs(self):
        """
        test that the command includes a removal of the temporary log files
        """
        bwa = Bwa()
        bwa.input_fastq1 = "foo_1.fq"
        bwa.input_fastq2 = None
        bwa.input_reference_sequence = "ref.fasta"
        bwa.readgroup = "__readgroup__"
        bwa.output = "out.bam"
        cmd = bwa.command()
        self.assertTrue(cmd.endswith('rm out.bam.bwa.log out.bam.samblaster.log'),
                        msg="rm temp logs must be the final part of the command")

    @patch('uuid.uuid4')
    def test_skewer_deletes_tmp(self, mock_uuid):
        """
        test that skewer includes commands to remove the temp dir it uses
        """
        mock_uuid.return_value = "foo"

        skewer = Skewer()
        skewer.input1 = "in_1.fq.gz"
        skewer.input2 = "in_2.fq.gz"
        skewer.output1 = "out_1.fq.gz"
        skewer.output2 = "out_2.fq.gz"
        skewer.stats = "stats.txt"
        skewer.scratch = "/scratch"
        cmd = skewer.command()

        self.assertTrue(cmd.endswith('rm -r /scratch/skewer-foo'))

    @patch('uuid.uuid4')
    def test_skewer_pe_uses_both_fqs(self, mock_uuid):
        """
        test that both input1 and input2 are used in the command line
        and that both output1 and output2 are copied back
        """
        mock_uuid.return_value = "foo"

        skewer = Skewer()
        skewer.input1 = "in_1.fq.gz"
        skewer.input2 = "in_2.fq.gz"
        skewer.output1 = "/path/to/out_1.fq.gz"
        skewer.output2 = "/path/to/out_2.fq.gz"
        skewer.stats = "stats.txt"
        skewer.scratch = "/scratch"
        cmd = skewer.command()
        self.assertIn("in_1.fq.gz  in_2.fq.gz", cmd)
        self.assertIn("cp /scratch/skewer-foo/skewer-trimmed-pair1.fastq.gz /path/to/out_1.fq.gz", cmd)
        self.assertIn("cp /scratch/skewer-foo/skewer-trimmed-pair2.fastq.gz /path/to/out_2.fq.gz", cmd)

    @patch('uuid.uuid4')
    def test_skewer_se_only_uses_input1(self, mock_uuid):
        """
        test that skewer doesn't try to copy back fastq2 in single-end mode
        """
        mock_uuid.return_value = "foo"
        skewer = Skewer()
        skewer.input1 = "in_1.fq.gz"
        skewer.input2 = None
        skewer.output1 = "/path/to/out_1.fq.gz"
        skewer.output2 = "/path/to/out_2.fq.gz"
        skewer.stats = "stats.txt"
        skewer.scratch = "/scratch"
        cmd = skewer.command()

        self.assertNotIn("/path/to/out_2.fq.gz", cmd)

    def test_skewer_raises_valuerror_if_output_not_gzipped(self):
        """
        test that skewer raises a ValueError if the user tries to use uncompressed output files
        """
        skewer = Skewer()
        skewer.input1 = "in_1.fq.gz"
        skewer.input2 = None
        skewer.output1 = "/path/to/out_1.fq"
        skewer.output2 = "/path/to/out_2.fq.gz"
        skewer.stats = "stats.txt"
        skewer.scratch = "/scratch"
        with self.assertRaises(ValueError):
            skewer.command()

        skewer.output1 = "/path/to/out_1.fq.gz"
        skewer.output2 = "/path/to/out_2.fq"
        with self.assertRaises(ValueError):
            skewer.command()

    @patch('autoseq.tools.alignment.align_se')
    def test_align_library_se(self, mock_align_se):
        mock_align_se.return_value = "dummy_bwa_output_se"
        bwa_output = align_library(self.test_clinseq_pipeline, ["test.fq.gz"], [], "AL-P-NA12877-T-03098849-TD1-TT1",
                                   "dummy_reference.fasta", "dummy_output_dir")
        self.assertTrue(mock_align_se.called)
        self.assertEquals(bwa_output, "dummy_bwa_output_se")

    @patch('autoseq.tools.alignment.align_pe')
    def test_align_library_pe(self, mock_align_pe):
        mock_align_pe.return_value = "dummy_bwa_output_pe"
        bwa_output = align_library(self.test_clinseq_pipeline, ["test1.fq.gz"], ["test2.fq.gz"], "AL-P-NA12877-T-03098849-TD1-TT1",
                                   "dummy_reference.fasta", "dummy_output_dir")
        self.assertTrue(mock_align_pe.called)
        self.assertEquals(bwa_output, "dummy_bwa_output_pe")

    def test_align_se(self):
        bwa_output = align_se(self.test_clinseq_pipeline, ["test.fq.gz"], "AL-P-NA12877-T-03098849-TD1-TT1",
                              "dummy_reference.fasta", "dummy_output_dir", 1)
        self.assertEquals(len(self.test_clinseq_pipeline.graph.nodes()), 3)
        self.assertEquals(bwa_output.split(".")[-1], "bam")

    def test_align_pe(self):
        bwa_output = align_pe(self.test_clinseq_pipeline, ["test1.fq.gz"], ["test2.fq.gz"],
                              "AL-P-NA12877-T-03098849-TD1-TT1", "dummy_reference.fasta",
                              "dummy_output_dir", 1)
        self.assertEquals(len(self.test_clinseq_pipeline.graph.nodes()), 4)
        self.assertEquals(bwa_output.split(".")[-1], "bam")
