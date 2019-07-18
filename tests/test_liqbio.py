# 
# 
# def test_make_sample_dict(self):
#     libs = parse_orderform(self.orderform)
#     sample_dicts = make_sample_dicts(libs)
# 
#     self.assertIn('NA12877', sample_dicts)
# 
#     sample_dict_NA12877 = {
#         "sdid": "NA12877",
#         "panel": {
#             "T": "NA12877-T-03098849-TD1-TT1",
#             "N": "NA12877-N-03098121-TD1-TT1",
#             "CFDNA": ["NA12877-CFDNA-03098850-TD1-TT1", "NA12877-CFDNA-03098850-TD1-TT2"]
#         },
#         "wgs": {
#             "T": "NA12877-T-03098849-TD1-WGS",
#             "N": "NA12877-N-03098121-TD1-WGS",
#             "CFDNA": ["NA12877-CFDNA-03098850-TD1-WGS"]
#         }
#     }
# 
#     self.assertEqual(sample_dict_NA12877, sample_dicts['NA12877'])
