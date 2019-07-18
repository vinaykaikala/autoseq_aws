import unittest
from mock import MagicMock, patch

from autoseq.util.orderform import *


class TestOrderform(unittest.TestCase):
    def test_parse_orderform_block_valid(self):
        fields_to_parse = ["<SAMPLE ENTRIES>",
                           "LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000",
                           "</SAMPLE ENTRIES>"]
        parsed_clinseq_barcodes = parse_orderform_block(fields_to_parse)
        self.assertEquals(parsed_clinseq_barcodes,
                          ["LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000"])

    def test_parse_orderform_block_empty(self):
        fields_to_parse = ["<SAMPLE ENTRIES>",
                           "</SAMPLE ENTRIES>"]
        parsed_clinseq_barcodes = parse_orderform_block(fields_to_parse)
        self.assertEquals(parsed_clinseq_barcodes, [])

    def test_parse_orderform_block_no_start(self):
        fields_to_parse = ["LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000",
                           "</SAMPLE ENTRIES>"]
        parsed_clinseq_barcodes = parse_orderform_block(fields_to_parse)
        self.assertEquals(parsed_clinseq_barcodes, [])

    def test_parse_orderform_block_no_end(self):
        fields_to_parse = ["<SAMPLE ENTRIES>",
                           "LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000"]
        parsed_clinseq_barcodes = parse_orderform_block(fields_to_parse)
        self.assertEquals(parsed_clinseq_barcodes,
                          ["LB-P-00000001-CFDNA-01234567-TP201701011540-CM2017001022000"])

    def test_parse_orderform_worksheet(self):
        # FIXME: Fiddly to mock behaviour of order_form_worksheet => Skipping proper testing presently.
        dummy_worksheet = MagicMock()
        dummy_worksheet.iter_rows = MagicMock()
        dummy_worksheet.iter_rows.return_value = []
        extracted_barcodes = parse_orderform_worksheet(dummy_worksheet)
        self.assertEquals(extracted_barcodes, [])

    @patch('autoseq.util.orderform.load_workbook')
    def test_parse_orderform(self, mock_load_workbook):
        # FIXME: Fiddly to mock behaviour of order_form_worksheet => Skipping proper testing presently.
        dummy_worksheet = MagicMock()
        dummy_worksheet.iter_rows = MagicMock()
        dummy_worksheet.iter_rows.return_value = []
        mock_load_workbook.return_value = dummy_worksheet

        extracted_barcodes = parse_orderform(dummy_worksheet)
        self.assertEquals(extracted_barcodes, [])
