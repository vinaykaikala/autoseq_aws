import unittest
from autoseq.tools.purity import *


class TestPurity(unittest.TestCase):
    def test_purecn(self):
        purecn = PureCN("input.seg", "input.vcf", "tumorID", "dummy-dir", gcgene_file = "dummy-file.txt",
                        minpurity = 0.03, postopt = True)
        cmd = purecn.command()
        self.assertIn('input.seg', cmd)
        self.assertIn('input.vcf', cmd)
        self.assertIn('dummy-dir', cmd)
        self.assertIn('postoptimize', cmd)
        self.assertIn('0.03', cmd)
        self.assertIn('source activate purecn-env', cmd)