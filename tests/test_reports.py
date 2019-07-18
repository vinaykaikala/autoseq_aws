import unittest
from autoseq.tools.reports import *


class TestReports(unittest.TestCase):
    def test_compile_metadata(self):
        test_job = CompileMetadata("dummy.json", "1234567", "7654321", "dummy2.json", "dummy_addresses.csv")
        cmd = test_job.command()
        self.assertIn('dummy.json', cmd)
        self.assertIn('1234567', cmd)
        self.assertIn('7654321', cmd)
        self.assertIn('dummy2.json', cmd)
        self.assertIn('dummy_addresses.csv', cmd)

    def test_compile_alascca_genomic_json(self):
        test_job = CompileAlasccaGenomicJson("dummy.vcf", "dummy_cnv.txt", "dummy_msi_sensor.txt",
                                             "dummy_purity.json", "dummy_contam.json", "dummy_tcov.json",
                                             "dummy_ncov.json", "output.json")
        cmd = test_job.command()
        self.assertIn('dummy.vcf', cmd)
        self.assertIn('dummy_cnv.txt', cmd)
        self.assertIn('dummy_msi_sensor.txt', cmd)
        self.assertIn('dummy_purity.json', cmd)
        self.assertIn('dummy_contam.json', cmd)
        self.assertIn('dummy_tcov.json', cmd)
        self.assertIn('dummy_ncov.json', cmd)
        self.assertIn('output.json', cmd)

    def test_write_alascca_report_full(self):
        test_job = WriteAlasccaReport("dummy_genomic.json", "dummy_metadata.json", "output.pdf")
        cmd = test_job.command()
        self.assertIn('dummy_genomic.json', cmd)
        self.assertIn('dummy_metadata.json', cmd)
        self.assertIn('output.pdf', cmd)
        self.assertNotIn("--alascca_only", cmd)

    def test_write_alascca_report_alascca_only(self):
        test_job = WriteAlasccaReport("dummy_genomic.json", "dummy_metadata.json", "output.pdf", True)
        cmd = test_job.command()
        self.assertIn('dummy_genomic.json', cmd)
        self.assertIn('dummy_metadata.json', cmd)
        self.assertIn('output.pdf', cmd)
        self.assertIn("--alascca_only", cmd)
