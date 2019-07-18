import pdb, unittest

# from mock import patch, Mock
# from autoseq.pipeline.alascca import AlasccaPipeline
# 
# 
# class TestAlascca(unittest.TestCase):
#     @patch('autoseq.pipeline.alascca.MultiQC')
#     @patch('autoseq.pipeline.alascca.AlasccaPipeline.add')
#     @patch('autoseq.pipeline.alascca.AlasccaPipeline.run_panel_bam_qc')
#     @patch('autoseq.pipeline.alascca.AlasccaPipeline.run_fastq_qc')
#     @patch('autoseq.pipeline.alascca.AlasccaPipeline.analyze_panel')
#     @patch('autoseq.pipeline.alascca.AlasccaPipeline.get_all_fastqs')
#     @patch('autoseq.pipeline.alascca.get_libdict')
#     @patch('autoseq.pipeline.alascca.find_fastqs')
#     def setUp(self, mock_find_fastqs, mock_get_libdict, mock_get_all_fastqs, mock_analyze_panel, mock_run_fastq_qc, mock_run_panel_bam_qc, mock_add, mock_multiQC):
#         mock_find_fastqs.return_value = ["dummy_fastq_filename"]
#         mock_get_libdict.return_value = {"capture_kit_name": "big_design"}
#         mock_get_all_fastqs.return_value = None
#         mock_analyze_panel.return_value = {"tbam": "dummy", "nbam": "dummy"}
#         mock_run_fastq_qc.return_value = []
#         mock_run_panel_bam_qc.return_value = []
#         mock_add.return_value = None
# 
#         # Construct a AlasccaPipeline that is set up well enough for unit
#         # testing to be run on it's instance methods:
#         self.dummy_path_to_big_vcf = "dummy_path_to_big_vcf"
#         self.dummy_path_to_V3V4_vcf = "dummy_path_to_V3V4_vcf"
#         self.dummy_path_to_V3V4big_vcf = "dummy_path_to_V3V4big_vcf"
#         self.dummy_path_to_test_regions_vcf = "dummy_path_to_test_regions_vcf"
#         dummy_refdata = {"contest_vcfs": {"big": self.dummy_path_to_big_vcf,
#                                           "clinseqV3V4": self.dummy_path_to_V3V4_vcf,
#                                           "clinseqV3V4big": self.dummy_path_to_V3V4big_vcf,
#                                           "test-regions": self.dummy_path_to_test_regions_vcf}}
# 
#         self.dummy_alascca_pipeline = AlasccaPipeline({"panel":{"T":0, "N":0}},
#                                                       dummy_refdata,
#                                                       '/tmp/autoseq-test',
#                                                       '/tmp')
# 
#     @patch('autoseq.pipeline.alascca.get_libdict')
#     def test_get_pop_af_vcf_both_cb(self, mock_get_libdict):
#         mock_get_libdict.side_effect = [{"capture_kit_name": "big_design"}, {"capture_kit_name": "big_design"}]
# 
#         self.assertEquals(self.dummy_alascca_pipeline.get_pop_af_vcf(), self.dummy_path_to_big_vcf)
# 
#     @patch('autoseq.pipeline.alascca.get_libdict')
#     def test_get_pop_af_vcf_tcs_ncs(self, mock_get_libdict):
#         mock_get_libdict.side_effect = [{"capture_kit_name": "clinseq_v3_targets"}, {"capture_kit_name": "clinseq_v3_targets"}]
# 
#         self.assertEquals(self.dummy_alascca_pipeline.get_pop_af_vcf(), self.dummy_path_to_V3V4_vcf)
# 
#     @patch('autoseq.pipeline.alascca.get_libdict')
#     def test_get_pop_af_vcf_tcs_ncz(self, mock_get_libdict):
#         mock_get_libdict.side_effect = [{"capture_kit_name": "clinseq_v3_targets"}, {"capture_kit_name": "clinseq_v4"}]
# 
#         self.assertEquals(self.dummy_alascca_pipeline.get_pop_af_vcf(), self.dummy_path_to_V3V4_vcf)
# 
#     @patch('autoseq.pipeline.alascca.get_libdict')
#     def test_get_pop_af_vcf_tcb_ncz(self, mock_get_libdict):
#         mock_get_libdict.side_effect = [{"capture_kit_name": "big_design"}, {"capture_kit_name": "clinseq_v4"}]
# 
#         self.assertEquals(self.dummy_alascca_pipeline.get_pop_af_vcf(), self.dummy_path_to_V3V4big_vcf)
# 
#     @patch('autoseq.pipeline.alascca.get_libdict')
#     def test_get_pop_af_vcf_tcs_ncb(self, mock_get_libdict):
#         mock_get_libdict.side_effect = [{"capture_kit_name": "clinseq_v3_targets"}, {"capture_kit_name": "big_design"}]
# 
#         self.assertEquals(self.dummy_alascca_pipeline.get_pop_af_vcf(), self.dummy_path_to_V3V4big_vcf)
# 
#     @patch('autoseq.pipeline.alascca.get_libdict')
#     def test_get_pop_af_vcf_invalid_combo(self, mock_get_libdict):
#         mock_get_libdict.side_effect = [{"capture_kit_name": "core_design"}, {"capture_kit_name": "big_design"}]
# 
#         self.assertRaises(ValueError, lambda: self.dummy_alascca_pipeline.get_pop_af_vcf())

#    def test_get_pop_af_vcf_tpanel_cs_npanel_cs(self):
