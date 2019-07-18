import json
import os
import unittest

import subprocess
from genomicassertions.readassertions import ReadAssertions
from genomicassertions.variantassertions import VariantAssertions

from autoseq.tests import alascca_test_outdir, alascca_purity_test_outdir


class TestAlasccaPurity(unittest.TestCase, VariantAssertions, ReadAssertions):
    returncode = None
    tmpdir = None
    somatic_vcf = None
    blood_barcode = '03098849'
    tumor_barcode = '03098121'
    outdir = alascca_purity_test_outdir
    jobdb = os.path.join(outdir, "jobdb.json")

    @classmethod
    def setUpClass(cls):
        cmd_str = "autoseq " + \
                  " --ref /tmp/test-genome/autoseq-genome.json " + \
                  " --outdir {} ".format(cls.outdir) + \
                  " --libdir /tmp/libraries " + \
                  " --scratch /scratch/tmp/autoseq-integration-tests/alascca --jobdb {} --cores 2 alascca ".format(cls.jobdb) + \
                  " tests/alascca-test-sample_low_purity.json"
        subprocess.check_call(cmd_str, shell=True)

    def test_purity_estimate_file(self):
        purity_json_fn = "{}/variants/AL-P-NA12877-T-03098849-TD-TT-alascca-purity.json".format(self.outdir)
        with open(purity_json_fn, 'r') as fh:
            # NOTE: The purity call is no longer "FAIL" for this example, as the number of SNPs
            # is insufficient to estimate purity => It should be "OK" now. Need to run a larger
            # test (with more SNPs) or change the minimum SNP threshold in alasccaCNA.R in order
            # to test a purity "FAIL" example.
            purity_json = json.load(fh)
            self.assertEqual(purity_json['CALL'], 'OK')


class TestAlascca(unittest.TestCase, VariantAssertions, ReadAssertions):
    returncode = None
    tmpdir = None
    somatic_vcf = None
    blood_barcode = '03098849'
    tumor_barcode = '03098121'
    outdir = alascca_test_outdir
    jobdb = os.path.join(outdir, "jobdb.json")

    @classmethod
    def setUpClass(cls):
        cmd_str = "autoseq " + \
                  " --ref /tmp/test-genome/autoseq-genome.json " + \
                  " --outdir {} ".format(cls.outdir) + \
                  " --job-params tests/alascca-test-sample-params.json " + \
                  " --libdir /tmp/libraries " + \
                  " --scratch /scratch/tmp/autoseq-integration-tests/alascca --jobdb {} --cores 2 alascca ".format(cls.jobdb) + \
                  " tests/alascca-test-sample.json"
        subprocess.check_call(cmd_str, shell=True)

    def test_jobdb(self):
        jobdb = json.load(open(self.jobdb, 'r'))
        self.assertEqual(set([job['status'] for job in jobdb['jobs']]),
                         {'COMPLETED'})

    def test_vardict_somatic(self):
        vcf = os.path.join(self.outdir, "variants",
                           "AL-P-NA12877-T-03098849-TD-TT-AL-P-NA12877-N-03098121-TD-TT.vardict-somatic.vcf.gz")

        # TP53 insertion: MU2185182, chr17:g.7578475->G
        # TP53 deletion: MU25947, chr17:g.7577558G>-
        # TP53 DNV: MU52971976, chr17:g.7574003GG>AA
        # PIK3CA hotspot E545K, MU5219, chr3:g.178936091G>A
        # PTEN hotspot R130Q, MU29098, chr10:g.89692905G>A
        # PTEN hotspot R233*, MU589331, chr10:g.89717672C>T
        # AR intron variant, MU50988553, chrX:g.66788924G>A

        self.assertVcfHasSample(vcf, 'AL-P-NA12877-N-03098121')  # N lib id is set
        self.assertVcfHasSample(vcf, 'AL-P-NA12877-T-03098849')  # T lib id is set

        # deletion is called
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, 17, 7577557, 'AG', 'A')

        # insertion is called
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, 17, 7578474, 'C', 'CG')

        # dnv (GG>AA) is called as a single event
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, 17, 7574003, 'GG', 'AA')

        # PTEN hotspots are called
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, 10, 89692905, 'G', 'A')
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, 10, 89717672, 'C', 'T')

        # PIK3CA hotspot is called
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, 3, 178936091, 'G', 'A')

    def test_panel_bam_readgroups(self):
        bam = os.path.join(self.outdir, "bams", "TT", "AL-P-NA12877-T-03098849-TD1-TT1.bam")
        self.assertBamHeaderElementEquals(bam, 'RG', [{'ID': 'AL-P-NA12877-T-03098849-TD1-TT1',
                                                       'SM': 'AL-P-NA12877-T-03098849',
                                                       'LB': 'TD1',
                                                       'PL': 'ILLUMINA'}])

        bam = os.path.join(self.outdir, "bams", "TT", "AL-P-NA12877-N-03098121-TD1-TT1.bam")
        self.assertBamHeaderElementEquals(bam, 'RG', [{'ID': 'AL-P-NA12877-N-03098121-TD1-TT1',
                                                       'SM': 'AL-P-NA12877-N-03098121',
                                                       'LB': 'TD1',
                                                       'PL': 'ILLUMINA'}])

    def test_germline_vcf(self):
        vcf = os.path.join(self.outdir, "variants",
                           "AL-P-NA12877-N-03098121-TD-TT.freebayes-germline.vcf.gz")

        self.assertVcfHasSample(vcf, 'AL-P-NA12877-N-03098121')  # sample is the name of the normal lib id
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, '3', 178925677, 'G', 'A')  # SNP
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, '17', 7579643, 'CCCCCAGCCCTCCAGGT', 'C')  # deletion

    def test_msisensor(self):
        with open(os.path.join(self.outdir, "msisensor-AL-P-NA12877-N-03098121-TD-TT-AL-P-NA12877-T-03098849-TD-TT.tsv")) as fh:
            header = fh.readline().strip()
            self.assertTrue(header.startswith("Total_Number_of_Sites"))
            dataln = fh.readline().strip()
            self.assertTrue(dataln.startswith("20"))

    def test_report_metadata(self):
        metadata_json_fn = "{}/report/{}-{}.metadata.json".format(self.outdir, self.tumor_barcode, self.blood_barcode)
        with open(metadata_json_fn, 'r') as fh:
            metadata_json = json.load(fh)
            self.assertEqual(metadata_json['blood_referral_ID'], 159725)
            self.assertEqual(metadata_json['tumor_referral_ID'], 159977)
            self.assertEqual(len(metadata_json['return_addresses']), 2)

    def test_genomic_json(self):
        genomic_json_fn = "{}/report/{}-{}.genomic.json".format(self.outdir, self.tumor_barcode, self.blood_barcode)
        with open(genomic_json_fn, 'r') as fh:
            genomic_json = json.load(fh)

            # sample is of class A
            self.assertEqual(genomic_json['alascca_class_report']['alascca_class'], 'Mutation class A')

            # MSI status is not determined since only 20 markers are in there
            self.assertEqual(genomic_json["msi_report"]["msi_status"], "Not determined")

            # no mutations in BRAF, KRAS or NRAS
            self.assertDictEqual(genomic_json["simple_somatic_mutations_report"],
                                 {
                                     "BRAF": {
                                         "alterations": [],
                                         "status": "Not determined"
                                     },
                                     "KRAS": {
                                         "alterations": [],
                                         "status": "Not determined"
                                     },
                                     "NRAS": {
                                         "alterations": [],
                                         "status": "Not determined"
                                     }
                                 }
                                 )

    def test_purity_estimate_file(self):
        purity_json_fn = "{}/variants/AL-P-NA12877-T-03098849-TD-TT-alascca-purity.json".format(self.outdir)
        with open(purity_json_fn, 'r') as fh:
            purity_json = json.load(fh)
            self.assertEqual(purity_json['CALL'], 'OK')
