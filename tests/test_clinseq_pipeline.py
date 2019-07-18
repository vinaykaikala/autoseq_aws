import unittest
import itertools
from mock import patch
from autoseq.pipeline.clinseq import *
from autoseq.util.clinseq_barcode import UniqueCapture


class TestClinseq(unittest.TestCase):
    def setUp(self):
        self.test_single_panel_results = SinglePanelResults()
        self.test_cancer_vs_normal_results = CancerVsNormalPanelResults()
        self.sample_data = {
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
                }
            },
            "contest_vcfs": {
                "test-regions": "test_contest.vcf"
            },
            "vep_dir": None
        }
        self.test_cancer_capture = UniqueCapture("AL", "P-NA12877", "CFDNA", "03098850", "TD", "TT")
        self.test_normal_capture = UniqueCapture("AL", "P-NA12877", "N", "03098121", "TD", "TT")
        self.test_wg_capture = UniqueCapture("AL", "P-NA12877", "N", "03098121", "TD", "WG")
        self.test_clinseq_pipeline = ClinseqPipeline(self.sample_data, self.ref_data, {"cov-low-thresh-fraction": 0.8}, "/tmp", "/nfs/LIQBIO/INBOX/exomes")

    def test_single_panel_results_output(self):
        self.assertEquals(self.test_single_panel_results.merged_bamfile, None)

    def test_cancer_against_normal_results(self):
        self.assertEquals(self.test_cancer_vs_normal_results.somatic_vcf, None)

    def test_pipeline_constructor(self):
        self.assertEquals(self.test_clinseq_pipeline.job_params, {'cov-low-thresh-fraction': 0.8})
        self.assertEquals(self.test_clinseq_pipeline.sampledata["sdid"], "P-NA12877")
        self.assertEquals(self.test_clinseq_pipeline.refdata["icgc"], "variants/icgc_release_20_simple_somatic_mutation.aggregated.vcf.gz")

    def test_get_job_param_set(self):
        self.assertEquals(self.test_clinseq_pipeline.get_job_param("cov-low-thresh-fraction"), 0.8)

    def test_get_job_param_default(self):
        self.assertEquals(self.test_clinseq_pipeline.get_job_param("cov-high-thresh-fold-cov"), 100)

    def test_set_germline_vcf(self):
        self.test_clinseq_pipeline.set_germline_vcf(self.test_cancer_capture, ("test1.vcf", "test2.vcf"))
        self.assertEquals(self.test_clinseq_pipeline.normal_capture_to_vcf[self.test_cancer_capture],
                          ("test1.vcf", "test2.vcf"))

    def test_get_germline_vcf_exists(self):
        self.test_clinseq_pipeline.set_germline_vcf(self.test_cancer_capture, ("test1.vcf", "test2.vcf"))
        self.assertEquals(self.test_clinseq_pipeline.get_germline_vcf(self.test_cancer_capture), "test1.vcf")

    def test_get_germline_vcf_none(self):
        self.assertEquals(self.test_clinseq_pipeline.get_germline_vcf(self.test_cancer_capture), None)

    def test_set_capture_bam(self):
        self.test_clinseq_pipeline.set_capture_bam(self.test_cancer_capture, "test.bam")
        self.assertEquals(self.test_clinseq_pipeline.capture_to_results[self.test_cancer_capture].merged_bamfile,
                          "test.bam")

    def test_set_capture_cnr(self):
        self.test_clinseq_pipeline.set_capture_cnr(self.test_cancer_capture, "test.cnr")
        self.assertEquals(self.test_clinseq_pipeline.capture_to_results[self.test_cancer_capture].cnr,
                          "test.cnr")

    def test_set_capture_cns(self):
        self.test_clinseq_pipeline.set_capture_cns(self.test_cancer_capture, "test.cns")
        self.assertEquals(self.test_clinseq_pipeline.capture_to_results[self.test_cancer_capture].cns,
                          "test.cns")

    def test_get_capture_bam_exists(self):
        self.test_clinseq_pipeline.set_capture_bam(self.test_cancer_capture, "test.bam")
        self.assertEquals(self.test_clinseq_pipeline.get_capture_bam(self.test_cancer_capture),
                          "test.bam")

    def test_get_capture_bam_none(self):
        self.assertEquals(self.test_clinseq_pipeline.get_capture_bam(self.test_cancer_capture),
                          None)


    @patch('autoseq.pipeline.clinseq.data_available_for_clinseq_barcode')
    def test_check_sampledata_all_available(self, mock_data_available_for_clinseq_barcode):
        mock_data_available_for_clinseq_barcode.return_value = True
        # Test that no changes occur to the sample data if the data is all available:
        self.test_clinseq_pipeline.check_sampledata()
        self.assertEquals(self.test_clinseq_pipeline.sampledata,
                          {
                              "sdid": "P-NA12877",
                              "T": ["AL-P-NA12877-T-03098849-TD1-TT1", "AL-P-NA12877-T-03098849-TD1-WGS"],
                              "N": ["AL-P-NA12877-N-03098121-TD1-TT1", "AL-P-NA12877-N-03098121-TD1-WGS"],
                              "CFDNA": ["LB-P-NA12877-CFDNA-03098850-TD1-TT1", "LB-P-NA12877-CFDNA-03098850-TD1-TT2",
                                        "LB-P-NA12877-CFDNA-03098850-TD1-WGS"]
                          })

    @patch('autoseq.pipeline.clinseq.data_available_for_clinseq_barcode')
    def test_check_sampledata_none_available(self, mock_data_available_for_clinseq_barcode):
        mock_data_available_for_clinseq_barcode.return_value = False
        # Test that no changes occur to the sample data if the data is all available:
        self.test_clinseq_pipeline.check_sampledata()
        self.assertEquals(self.test_clinseq_pipeline.sampledata,
                          {
                              "sdid": "P-NA12877",
                              "T": [],
                              "N": [],
                              "CFDNA": []
                          })

    def test_vep_is_set(self):
        self.assertEquals(self.test_clinseq_pipeline.vep_data_is_available(), False)

    def test_get_all_unique_capture(self):
        self.assertEquals(self.test_clinseq_pipeline.get_mapped_captures_all(), [])

    def test_get_unique_captures_no_wgs(self):
        self.test_clinseq_pipeline.capture_to_results = {self.test_cancer_capture: 1}
        self.assertEquals(self.test_clinseq_pipeline.get_mapped_captures_no_wgs(), [self.test_cancer_capture])

    def test_get_unique_captures_only_wgs(self):
        self.test_clinseq_pipeline.capture_to_results = {self.test_cancer_capture: 1}
        self.assertEquals(self.test_clinseq_pipeline.get_mapped_captures_only_wgs(), [])

    def test_get_unique_normal_captures(self):
        self.test_clinseq_pipeline.capture_to_results = {self.test_cancer_capture: 1}
        self.assertEquals(self.test_clinseq_pipeline.get_mapped_captures_normal(), [])

    def test_get_unique_cancer_captures(self):
        self.test_clinseq_pipeline.capture_to_results = {self.test_cancer_capture: 1}
        self.assertEquals(self.test_clinseq_pipeline.get_mapped_captures_cancer(), [self.test_cancer_capture])

    def test_get_prep_kit_name(self):
        self.assertEquals(self.test_clinseq_pipeline.get_prep_kit_name("BN"), "BIOO_NEXTFLEX")

    def test_get_capture_name(self):
        self.assertEquals(self.test_clinseq_pipeline.get_capture_name("CM"), "monitor")

    def test_get_capture_name_wg(self):
        self.assertEquals(self.test_clinseq_pipeline.get_capture_name("WG"), "lowpass_wgs")

    def test_get_all_clinseq_barcodes(self):
        self.assertEquals(self.test_clinseq_pipeline.get_all_clinseq_barcodes(),
                          ["AL-P-NA12877-T-03098849-TD1-TT1", "AL-P-NA12877-T-03098849-TD1-WGS",
                           "AL-P-NA12877-N-03098121-TD1-TT1", "AL-P-NA12877-N-03098121-TD1-WGS",
                           "LB-P-NA12877-CFDNA-03098850-TD1-TT1", "LB-P-NA12877-CFDNA-03098850-TD1-TT2",
                           "LB-P-NA12877-CFDNA-03098850-TD1-WGS"])

    def test_get_unique_capture_to_clinseq_barcodes(self):
        l1 = list(itertools.chain.from_iterable(self.test_clinseq_pipeline.get_unique_capture_to_clinseq_barcodes().values()))
        l2 = ["AL-P-NA12877-T-03098849-TD1-TT1", "AL-P-NA12877-T-03098849-TD1-WGS",
              "AL-P-NA12877-N-03098121-TD1-TT1", "AL-P-NA12877-N-03098121-TD1-WGS",
              "LB-P-NA12877-CFDNA-03098850-TD1-TT1", "LB-P-NA12877-CFDNA-03098850-TD1-TT2",
              "LB-P-NA12877-CFDNA-03098850-TD1-WGS"]
        self.assertEquals(set(l1),
                          set(l2))

    def test_merge_and_rm_dup(self):
        self.test_clinseq_pipeline.merge_and_rm_dup(self.test_cancer_capture, ["test.bam"])
        self.assertNotEqual(
            self.test_clinseq_pipeline.capture_to_results[self.test_cancer_capture].merged_bamfile,
            None)
        self.assertEquals(\
            len(self.test_clinseq_pipeline.graph.nodes()), 2)
        self.assertEquals(\
            len(self.test_clinseq_pipeline.qc_files), 1)

#    @patch('autoseq.pipeline.clinseq.find_fastqs')
#    def test_configure_fastq_qcs(self, mock_find_fastqs):
#        mock_find_fastqs.return_value = ["dummy.fastq.gz"]
#        import pdb; pdb.set_trace()
#        self.assertEquals(len(self.test_clinseq_pipeline.configure_fastq_qcs()),
#                          len(self.test_clinseq_pipeline.get_all_clinseq_barcodes()))

    @patch('autoseq.pipeline.clinseq.align_library')
    @patch('autoseq.pipeline.clinseq.find_fastqs')
    def test_configure_align_and_merge(self, mock_find_fastqs, mock_align_library):
        mock_align_library.return_value = "dummy_merged.bam"
        mock_find_fastqs.return_value = "dummy.fastq.gz"
        self.test_clinseq_pipeline.configure_align_and_merge()
        self.assertTrue(mock_align_library.called)
        self.assertEquals(len(self.test_clinseq_pipeline.qc_files),
                          len(self.test_clinseq_pipeline.get_unique_capture_to_clinseq_barcodes()))

    def test_call_germline_variants(self):
        self.test_clinseq_pipeline.call_germline_variants(self.test_normal_capture, "test.bam")
        self.assertEquals(len(self.test_clinseq_pipeline.graph.nodes()), 1)

    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.call_germline_variants')
    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.get_mapped_captures_cancer')
    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.configure_panel_analysis_cancer_vs_normal')
    def test_configure_panel_analysis_with_normal(self, mock_configure_panel_analysis_cancer_vs_normal,
                                                  mock_get_mapped_captures_cancer,
                                                  mock_call_germline_variants):
        mock_get_mapped_captures_cancer.return_value = [self.test_cancer_capture]
        self.test_clinseq_pipeline.configure_panel_analysis_with_normal(self.test_normal_capture)
        self.assertTrue(mock_call_germline_variants.called)
        self.assertTrue(mock_configure_panel_analysis_cancer_vs_normal.called)

    def test_configure_panel_analysis_with_normal_invalid(self):
        self.assertRaises(ValueError,
                          lambda: self.test_clinseq_pipeline.configure_panel_analysis_with_normal(
                              self.test_cancer_capture))

    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.get_capture_name')
    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.get_capture_bam')
    def test_configure_single_capture_analysis(self, mock_get_capture_bam,
                                               mock_get_capture_name):
        mock_get_capture_name.return_value = "test-regions"
        mock_get_capture_bam.return_value = "test.bam"
        self.test_clinseq_pipeline.configure_single_capture_analysis(self.test_cancer_capture)
        self.assertEquals(len(self.test_clinseq_pipeline.graph.nodes()), 2)

    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.get_capture_name')
    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.get_capture_bam')
    def test_configure_single_capture_analysis_no_ref(self, mock_get_capture_bam,
                                                      mock_get_capture_name):
        mock_get_capture_name.return_value = "test-regions"
        mock_get_capture_bam.return_value = "test.bam"
        self.test_clinseq_pipeline.configure_single_capture_analysis(self.test_cancer_capture)
        self.assertEquals(len(self.test_clinseq_pipeline.graph.nodes()), 2)

    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.configure_single_wgs_analyses')
    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.get_mapped_captures_only_wgs')
    def test_configure_lowpass_analyses(self, mock_get_mapped_captures_only_wgs,
                                        mock_configure_single_wgs_analyses):
        mock_get_mapped_captures_only_wgs.return_value = [self.test_wg_capture]
        self.test_clinseq_pipeline.configure_lowpass_analyses()
        self.assertTrue(mock_configure_single_wgs_analyses.called)

    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.get_capture_bam')
    def test_configure_single_wgs_analyses(self, mock_get_capture_bam):
        mock_get_capture_bam.return_value = "test.bam"
        self.test_clinseq_pipeline.configure_single_wgs_analyses(self.test_wg_capture)
        self.assertEquals(len(self.test_clinseq_pipeline.graph.nodes()), 2)

    def test_run_wgs_bam_qc(self):
        self.assertEquals(len(self.test_clinseq_pipeline.run_wgs_bam_qc(["test1.bam", "test2.bam"])),
                          4)

    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.configure_single_capture_analysis')
    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.configure_panel_analysis_with_normal')
    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.get_mapped_captures_no_wgs')
    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.get_mapped_captures_normal')
    def test_configure_panel_analyses(self, mock_get_mapped_captures_normal,
                                      mock_get_mapped_captures_no_wgs,
                                      mock_configure_panel_analysis_with_normal,
                                      mock_configure_single_capture_analysis):
        mock_get_mapped_captures_normal.return_value = [self.test_normal_capture]
        mock_get_mapped_captures_no_wgs.return_value = [self.test_cancer_capture]
        self.test_clinseq_pipeline.configure_panel_analyses()
        self.assertTrue(mock_configure_single_capture_analysis.called)
        self.assertTrue(mock_configure_panel_analysis_with_normal.called)

    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.configure_single_capture_analysis')
    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.configure_panel_analysis_with_normal')
    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.get_mapped_captures_no_wgs')
    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.get_mapped_captures_normal')
    def test_configure_panel_analyses_cancer_only(self, mock_get_mapped_captures_normal,
                                                  mock_get_mapped_captures_no_wgs,
                                                  mock_configure_panel_analysis_with_normal,
                                                  mock_configure_single_capture_analysis):
        mock_get_mapped_captures_normal.return_value = []
        mock_get_mapped_captures_no_wgs.return_value = [self.test_cancer_capture]
        self.test_clinseq_pipeline.configure_panel_analyses()
        self.assertTrue(mock_configure_single_capture_analysis.called)
        self.assertFalse(mock_configure_panel_analysis_with_normal.called)

    @patch('autoseq.pipeline.clinseq.call_somatic_variants')
    def test_configure_somatic_calling(self, mock_call_somatic_variants):
        mock_call_somatic_variants.return_value = {"vardict": "test.vcf"}
        self.test_clinseq_pipeline.configure_somatic_calling(self.test_normal_capture,
                                                             self.test_cancer_capture)
        self.assertTrue(mock_call_somatic_variants.called)
        stored_vcf = self.test_clinseq_pipeline.normal_cancer_pair_to_results[
            (self.test_normal_capture, self.test_cancer_capture)].somatic_vcf
        self.assertEquals(stored_vcf, "test.vcf")

    def test_configure_vep_no_vep_data(self):
        self.assertRaises(ValueError, lambda: self.test_clinseq_pipeline.configure_vep(
            self.test_normal_capture, self.test_cancer_capture))

    def test_configure_vep_data_available(self):
        normal_capture_str = compose_lib_capture_str(self.test_normal_capture)
        cancer_capture_str = compose_lib_capture_str(self.test_cancer_capture)
        self.test_clinseq_pipeline.refdata['vep_dir'] = "dummy_vep_dir"
        self.test_clinseq_pipeline.configure_vep(self.test_normal_capture, self.test_cancer_capture)
        vepped_file = self.test_clinseq_pipeline.normal_cancer_pair_to_results[
            (self.test_normal_capture, self.test_cancer_capture)].vepped_vcf
        self.assertIn(cancer_capture_str, vepped_file)
        self.assertIn(normal_capture_str, vepped_file)

    def test_configure_vcf_add_sample(self):
        self.test_clinseq_pipeline.configure_vcf_add_sample(self.test_normal_capture,
                                                            self.test_cancer_capture)
        self.assertTrue(
            self.test_clinseq_pipeline.normal_cancer_pair_to_results[(
                self.test_normal_capture, self.test_cancer_capture)].vcf_addsample_output is not None)
        self.assertEquals(len(self.test_clinseq_pipeline.graph.nodes()), 1)

    def test_configure_msi_sensor(self):
        self.test_clinseq_pipeline.configure_msi_sensor(self.test_normal_capture,
                                                        self.test_cancer_capture)
        self.assertTrue(
            self.test_clinseq_pipeline.normal_cancer_pair_to_results[(
                self.test_normal_capture, self.test_cancer_capture)].msi_output is not None)
        self.assertEquals(len(self.test_clinseq_pipeline.graph.nodes()), 1)

    def test_configure_msings_invalid_refdata(self):
        self.assertRaises(InvalidRefDataException,
                          lambda: self.test_clinseq_pipeline.configure_msings(self.test_cancer_capture))


    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.get_capture_bam')
    def test_configure_msings(self, mock_get_capture_bam):
        capture_name = 'test-regions'
        self.ref_data['targets'][capture_name]['msings-baseline'] = "intervals/targets/test-regions.msings.baseline"
        self.ref_data['targets'][capture_name]['msings-bed'] = "intervals/targets/test-regions.msings.bed"
        self.ref_data['targets'][capture_name]['msings-msi_intervals'] = "intervals/targets/test-regions.msings.msi_intervals"
        test_clinseq_pipeline_2 = ClinseqPipeline(self.sample_data, self.ref_data, {"cov-low-thresh-fraction": 0.8}, "/tmp", "/nfs/LIQBIO/INBOX/exomes")
        mock_get_capture_bam.return_value = 'dummy.bam'
        test_clinseq_pipeline_2.configure_msings(self.test_cancer_capture)
        msings_result = test_clinseq_pipeline_2.capture_to_results[self.test_cancer_capture].msings_output
        self.assertTrue(msings_result is not None)
        self.assertEquals(len(test_clinseq_pipeline_2.graph.nodes()), 1)

    def test_configure_hz_conc(self):
        self.test_clinseq_pipeline.configure_hz_conc(self.test_normal_capture,
                                                     self.test_cancer_capture)
        self.assertTrue(
            self.test_clinseq_pipeline.normal_cancer_pair_to_results[(
                self.test_normal_capture, self.test_cancer_capture)].hzconcordance_output is not None)
        self.assertEquals(len(self.test_clinseq_pipeline.graph.nodes()), 1)

    def test_configure_contest_vcf_generation(self):
        self.test_clinseq_pipeline.configure_contest_vcf_generation(self.test_normal_capture,
                                                                    self.test_cancer_capture)
        self.assertEquals(len(self.test_clinseq_pipeline.graph.nodes()), 1)

    def test_configure_contest(self):
        self.test_clinseq_pipeline.configure_contest(self.test_normal_capture,
                                                     self.test_cancer_capture,
                                                     "test.vcf")
        self.assertEquals(len(self.test_clinseq_pipeline.graph.nodes()), 1)

    def test_configure_contam_qc_call(self):
        self.test_clinseq_pipeline.configure_contam_qc_call("dummy.txt",
                                                            self.test_cancer_capture)
        self.assertEquals(len(self.test_clinseq_pipeline.graph.nodes()), 1)

    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.configure_contest')
    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.configure_contam_qc_call')
    def test_configure_contamination_estimate(self, mock_configure_contam_qc_call, mock_configure_contest):
        mock_configure_contest.side_effect = ["cancer_vs_normal_output.txt", "normal_vs_cancer_output.txt"]
        mock_configure_contam_qc_call.return_value = "dummy_call.json"
        self.test_clinseq_pipeline.configure_contamination_estimate(self.test_normal_capture,
                                                                    self.test_cancer_capture)

        stored_normal_contest_output = self.test_clinseq_pipeline.normal_cancer_pair_to_results[
            (self.test_normal_capture, self.test_cancer_capture)].normal_contest_output
        self.assertEquals(stored_normal_contest_output, "normal_vs_cancer_output.txt")

        stored_cancer_contest_output = self.test_clinseq_pipeline.normal_cancer_pair_to_results[
            (self.test_normal_capture, self.test_cancer_capture)].cancer_contest_output
        self.assertEquals(stored_cancer_contest_output, "cancer_vs_normal_output.txt")

        stored_cancer_contam_call = self.test_clinseq_pipeline.normal_cancer_pair_to_results[
            (self.test_normal_capture, self.test_cancer_capture)].cancer_contam_call
        self.assertEquals(stored_cancer_contam_call, "dummy_call.json")

    def test_configure_panel_analysis_cancer_vs_normal(self):
        self.test_clinseq_pipeline.refdata['vep_dir'] = "dummy_vep_dir"
        self.test_clinseq_pipeline.configure_panel_analysis_cancer_vs_normal(self.test_normal_capture,
                                                                             self.test_cancer_capture)
        self.assertEquals(len(self.test_clinseq_pipeline.graph.nodes()), 10)

    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.get_mapped_captures_only_wgs')
    def test_configure_all_lowpass_qcs(self, mock_get_mapped_captures_only_wgs):
        mock_get_mapped_captures_only_wgs.return_value = [self.test_wg_capture]
        self.test_clinseq_pipeline.configure_all_lowpass_qcs()
        self.assertEquals(len(self.test_clinseq_pipeline.qc_files), 2)

    def test_configure_wgs_qc(self):
        qc_files = self.test_clinseq_pipeline.configure_wgs_qc(self.test_wg_capture)
        self.assertEquals(len(qc_files), 2)

    @patch('autoseq.pipeline.clinseq.ClinseqPipeline.get_mapped_captures_no_wgs')
    def test_configure_all_panel_qcs(self, mock_get_mapped_captures_no_wgs):
        mock_get_mapped_captures_no_wgs.return_value = [self.test_normal_capture, self.test_cancer_capture]
        self.test_clinseq_pipeline.configure_all_panel_qcs()
        self.assertEquals(len(self.test_clinseq_pipeline.qc_files), 12)

    def test_configure_multi_qc(self):
        self.test_clinseq_pipeline.configure_multi_qc()
        self.assertEquals(len(self.test_clinseq_pipeline.graph.nodes()), 1)

    def test_configure_panel_qc(self):
        qc_files = self.test_clinseq_pipeline.configure_panel_qc(self.test_cancer_capture)
        self.assertEquals(len(self.test_clinseq_pipeline.graph.nodes()), 6)
        self.assertEquals(len(qc_files), 6)
