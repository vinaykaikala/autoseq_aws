import unittest
from autoseq.tools.cnvcalling import *


class TestCNVCalling(unittest.TestCase):
    def test_qdna_seq(self):
        qdna_seq = QDNASeq("dummy.bam", "output_segments.txt")
        cmd = qdna_seq.command()
        self.assertIn('dummy.bam', cmd)
        self.assertIn('output_segments.txt', cmd)
        self.assertIn('qdnaseq.R', cmd)

    def test_qdnaseq2bed(self):
        qdnaseq2bed = QDNASeq2Bed("segments.txt", "output.bed", "genes.gtf")
        cmd = qdnaseq2bed.command()
        self.assertIn('genes.gtf', cmd)
        self.assertIn('segments.txt', cmd)
        self.assertIn('output.bed', cmd)
        self.assertIn('qdnaseq2bed.py', cmd)

    def test_alascca_cna_plot(self):
        alascca_cna_plot = AlasccaCNAPlot()
        alascca_cna_plot.input_cnr = "input.cnr"
        alascca_cna_plot.input_cns = "input.cns"
        alascca_cna_plot.input_germline_vcf = "input_germline.vcf"
        alascca_cna_plot.input_somatic_vcf = "input_somatic.vcf"
        alascca_cna_plot.chrsizes = "chrsizes.txt"
        alascca_cna_plot.output_png = "output.png"
        alascca_cna_plot.output_cna = "output.cna"
        alascca_cna_plot.output_purity = "purity.json"

        cmd = alascca_cna_plot.command()

        self.assertIn("input.cnr", cmd)
        self.assertIn("input.cns", cmd)
        self.assertIn("input_germline.vcf", cmd)
        self.assertIn("input_somatic.vcf", cmd)
        self.assertIn("chrsizes.txt", cmd)
        self.assertIn("output.png", cmd)
        self.assertIn("output.cna", cmd)
        self.assertIn("purity.json", cmd)

    def test_cnv_kit_reference(self):
        cnvkit = CNVkit("input.bam", "output.cns", "output.cnr", reference="dummy_reference.txt")
        cmd = cnvkit.command()
        self.assertIn('input.cns', cmd)
        self.assertIn('output.cns', cmd)
        self.assertIn('output.cnr', cmd)

    def test_cnv_kit_targets(self):
        cnvkit = CNVkit("input.bam", "output.cns", "output.cnr", targets_bed="dummy_targets.bed", fasta="dummy_fasta.fa")
        cmd = cnvkit.command()
        self.assertIn('input.cns', cmd)
        self.assertIn('output.cns', cmd)
        self.assertIn('output.cnr', cmd)
        self.assertIn('dummy_fasta.fa', cmd)
        self.assertIn('dummy_targets.bed', cmd)

    def test_cnv_kit_both(self):
        cnvkit = CNVkit("input.bam", "output.cns", "output.cnr", reference="dummy_reference.txt",
                        targets_bed="dummy_targets.bed")
        self.assertRaises(ValueError,
                          lambda: cnvkit.command())

    def test_cnv_kit_neither(self):
        cnvkit = CNVkit("input.bam", "output.cns", "output.cnr")
        self.assertRaises(ValueError,
                          lambda: cnvkit.command())

    def test_cns2seg(self):
        cns2seg = Cns2Seg("input.cns", "output.seg")
        cmd = cns2seg.command()
        self.assertIn('input.cns', cmd)
        self.assertIn('output.seg', cmd)
