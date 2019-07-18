import unittest
from autoseq.tools.variantcalling import *
from autoseq.pipeline.clinseq import ClinseqPipeline


class TestVariantCalling(unittest.TestCase):
    def setUp(self):
        sample_data = {
            "sdid": "P-NA12877",
            "T": ["AL-P-NA12877-T-03098849-TD1-TT1", "AL-P-NA12877-T-03098849-TD1-WGS"],
            "N": ["AL-P-NA12877-N-03098121-TD1-TT1", "AL-P-NA12877-N-03098121-TD1-WGS"],
            "CFDNA": ["LB-P-NA12877-CFDNA-03098850-TD1-TT1", "LB-P-NA12877-CFDNA-03098850-TD1-TT2",
                      "LB-P-NA12877-CFDNA-03098850-TD1-WGS"]
        }
        self.ref_data = {
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
            "ar_regions": "intervals/ar_regions.bed",
            "ts_regions": "intervals/ts_regions.bed",
            "fusion_regions": "intervals/fusion_regions.bed",
            "targets": {
                "test-regions": {
                    "cnvkit-ref": {
                        "THRUPLEX_PLASMASEQ": {
                            "CFDNA": "intervals/targets/progression.THRUPLEX_PLASMASEQ.CFDNA.cnn",
                            "N": "intervals/targets/progression.THRUPLEX_PLASMASEQ.N.cnn"
                        }
                    },
                    "cnvkit-fix": {
                        "THRUPLEX_PLASMASEQ": {
                            "CFDNA": "intervals/targets/progression.THRUPLEX_PLASMASEQ.CFDNA.cnvkit-fix.tsv"
                        }
                    },
                    "msisites": "intervals/targets/test-regions.msisites.tsv",
                    "targets-bed-slopped20": "intervals/targets/test-regions-GRCh37.slopped20.bed",
                    "targets-interval_list": "intervals/targets/test-regions-GRCh37.slopped20.interval_list",
                    "targets-interval_list-slopped20": "intervals/targets/test-regions-GRCh37.slopped20.interval_list",
                    "blacklist-bed": None,
                    "purecn_targets": "intervals/targets/purecn.bed",
                }
            },
            "contest_vcfs": {
                "test-regions": "test_contest.vcf"
            },
            "vep_dir": "dummy_vep_dir"
        }
        self.test_cancer_capture = UniqueCapture("AL", "P-NA12877", "CFDNA", "03098850", "TD", "TT")
        self.test_normal_capture = UniqueCapture("AL", "P-NA12877", "N", "03098121", "TD", "TT")
        self.test_clinseq_pipeline = ClinseqPipeline(sample_data, self.ref_data, {}, "/tmp", "/nfs/LIQBIO/INBOX/exomes")

    def test_freebayes(self):
        freebayes = Freebayes()
        freebayes.input_bams = ["input.bam"]
        freebayes.reference_sequence = "dummy.fasta"
        freebayes.target_bed = "dummy_targets.bed"
        freebayes.output = "output.txt"
        cmd = freebayes.command()
        self.assertIn('input.bam', cmd)
        self.assertIn('dummy.fasta', cmd)
        self.assertIn('dummy_targets.bed', cmd)
        self.assertIn('output.txt', cmd)

    def test_vardict(self):
        vardict = VarDict()
        vardict.input_tumor = "input_tumor.bam"
        vardict.input_normal = "input_normal.bam"
        vardict.tumorid = "tumor_id"
        vardict.normalid = "normal_id"
        vardict.reference_sequence = "dummy.fasta"
        vardict.reference_dict = {}
        vardict.target_bed = "dummy_targets.bed"
        vardict.output = "output.txt"
        cmd = vardict.command()
        self.assertIn('input_tumor.bam', cmd)
        self.assertIn('input_normal.bam', cmd)
        self.assertIn('tumor_id', cmd)
        self.assertIn('normal_id', cmd)
        self.assertIn('dummy.fasta', cmd)
        self.assertIn('dummy_targets.bed', cmd)
        self.assertIn('output.txt', cmd)

    def test_vardictForPureCN(self):
        vardict = VarDictForPureCN()
        vardict.input_tumor = "input_tumor.bam"
        vardict.input_normal = "input_normal.bam"
        vardict.tumorid = "tumor_id"
        vardict.normalid = "normal_id"
        vardict.reference_sequence = "dummy.fasta"
        vardict.reference_dict = {}
        vardict.target_bed = "dummy_targets.bed"
        vardict.dbsnp = "dummy_dbsnp.vcf.gz"
        vardict.output = "output.txt"
        cmd = vardict.command()
        self.assertIn('input_tumor.bam', cmd)
        self.assertIn('input_normal.bam', cmd)
        self.assertIn('tumor_id', cmd)
        self.assertIn('normal_id', cmd)
        self.assertIn('dummy.fasta', cmd)
        self.assertIn('dummy_targets.bed', cmd)
        self.assertIn('dummy_dbsnp.vcf.gz', cmd)
        self.assertIn('output.txt', cmd)

    def test_vep(self):
        vep = VEP()
        vep.input_vcf = "input.vcf"
        vep.output_vcf = "output.vcf"
        vep.reference_sequence = "dummy.fasta"
        vep.vep_dir = "dummy_dir"
        cmd = vep.command()
        self.assertIn('input.vcf', cmd)
        self.assertIn('dummy.fasta', cmd)
        self.assertIn('dummy_dir', cmd)
        self.assertIn('output.vcf', cmd)

    def test_vep_fork(self):
        vep = VEP()
        vep.input_vcf = "input.vcf"
        vep.output_vcf = "output.vcf"
        vep.reference_sequence = "dummy.fasta"
        vep.vep_dir = "dummy_dir"
        vep.threads = 2
        cmd = vep.command()
        self.assertIn('fork', cmd)

    def test_vep_gz(self):
        vep = VEP()
        vep.input_vcf = "input.vcf"
        vep.output_vcf = "output.vcf.gz"
        vep.reference_sequence = "dummy.fasta"
        vep.vep_dir = "dummy_dir"
        cmd = vep.command()
        self.assertIn('bgzip', cmd)

    def test_vcf_add_sample(self):
        vcf_add_sample = VcfAddSample()
        vcf_add_sample.input_vcf = "input.vcf"
        vcf_add_sample.input_bam = "input.bam"
        vcf_add_sample.samplename = "dummy_name"
        vcf_add_sample.output = "output.vcf"
        cmd = vcf_add_sample.command()
        self.assertIn('input.vcf', cmd)
        self.assertIn('input.bam', cmd)
        self.assertIn('dummy_name', cmd)
        self.assertIn('output.vcf', cmd)

    def test_vcf_add_sample_gz(self):
        vcf_add_sample = VcfAddSample()
        vcf_add_sample.input_vcf = "input.vcf"
        vcf_add_sample.input_bam = "input.bam"
        vcf_add_sample.samplename = "dummy_name"
        vcf_add_sample.output = "output.vcf.gz"
        cmd = vcf_add_sample.command()
        self.assertIn('bgzip > output.vcf.gz', cmd)
        self.assertIn('tabix', cmd)

    def test_vcf_filter(self):
        vcf_filter = VcfFilter()
        vcf_filter.input = "input.vcf"
        vcf_filter.filter = "test_filter"
        vcf_filter.output = "output.vcf"
        cmd = vcf_filter.command()
        self.assertIn('input.vcf', cmd)
        self.assertIn('test_filter', cmd)
        self.assertIn('output.vcf', cmd)

    def test_curl_split_and_left_align(self):
        curl_split_and_left_align = CurlSplitAndLeftAlign()
        curl_split_and_left_align.remote = "dummy_remote"
        curl_split_and_left_align.input_reference_sequence = "dummy_reference.fasta"
        curl_split_and_left_align.input_reference_sequence_fai = "dummy_reference.fasta.fai"
        curl_split_and_left_align.output = "output.vcf"
        cmd = curl_split_and_left_align.command()
        self.assertIn('dummy_remote', cmd)
        self.assertIn('dummy_reference.fasta', cmd)
        self.assertIn('output.vcf', cmd)

    def test_install_vep(self):
        install_vep = InstallVep()
        install_vep.output_dir = "dummy_output_dir"
        cmd = install_vep.command()
        self.assertIn('dummy_output_dir', cmd)

    def test_call_somatic_variants(self):
        num_jobs_before_call = len(self.test_clinseq_pipeline.graph.nodes())
        output_dict = call_somatic_variants(self.test_clinseq_pipeline, "test_cancer.bam", "test_normal.bam",
                                            self.test_cancer_capture, self.test_normal_capture, "test-regions",
                                            "test_outdir")
        self.assertEquals(output_dict.keys(), ["vardict", "freebayes"])
        num_jobs_after_call = len(self.test_clinseq_pipeline.graph.nodes())
        self.assertEquals(num_jobs_after_call, num_jobs_before_call + 2)
