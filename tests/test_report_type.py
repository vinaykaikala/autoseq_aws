import unittest
from mock import Mock, MagicMock, patch

from autoseq.util.report_type import *


class TestReportType(unittest.TestCase):
    def setUp(self):
        self.hospitalhundreds2country = {300: "country1",
                                         500: "country2"}
        
        self.country2reporttype = {"country1": "full",
                                   "country2": "alascca_only"}
    
    @patch('autoseq.util.report_type.create_sql_session')
    @patch('autoseq.util.report_type.query_database')
    def test_get_hospital_code_correct(self, mock_query_database, mock_create_sql_session):
        mock_referral_type = MagicMock()
        mock_db_config_file = MagicMock()
        mock_query_database.return_value.hospital_code = 353
        
        output = get_hospital_code("03098121", mock_referral_type, mock_db_config_file)
        self.assertEquals(output, 353)
    
    
    @patch('autoseq.util.report_type.create_sql_session')
    @patch('autoseq.util.report_type.query_database')
    def test_get_hospital_code_no_hospital_code(self, mock_query_database, mock_create_sql_session):
        mock_referral_type = MagicMock()
        mock_db_config_file = MagicMock()
        mock_query_database.return_value.hospital_code = None

        self.assertRaises(TypeError, get_hospital_code, "03098121", mock_referral_type, mock_db_config_file)


    def test_get_hospital_country_correct(self):
        hospital_code = 353
        output = get_hospital_country(hospital_code, self.hospitalhundreds2country)
        
        self.assertEquals(output, "country1")
    
    
    def test_get_hospital_country_no_int(self):
        hospital_code = "353"
        
        self.assertRaises(TypeError, get_hospital_country, hospital_code, self.hospitalhundreds2country)


    def test_get_hospital_country_unknown_code(self):
        hospital_code = 653
    
        self.assertRaises(KeyError, get_hospital_country, hospital_code, self.hospitalhundreds2country)


    def test_get_report_type_full(self):
        country = "country1"
        output = get_report_type(country, self.country2reporttype)
    
        self.assertEquals(output, "full")

    def test_get_report_type_alascca_only(self):
        country = "country2"
        output = get_report_type(country, self.country2reporttype)
    
        self.assertEquals(output, "alascca_only")


    def test_get_report_type_unknown_country(self):
        country = "country3"
    
        self.assertRaises(KeyError, get_report_type, country, self.country2reporttype)

    @patch('autoseq.util.report_type.get_hospital_code')
    def test_only_alascca_class_report_two_full(self, mock_get_hospital_code):
        mock_db_config_file = MagicMock()
        mock_get_hospital_code.side_effect = [353, 353]
        output = only_alascca_class_report("03098121", "03098849", mock_db_config_file)
    
        self.assertEquals(output, False)


    @patch('autoseq.util.report_type.get_hospital_code')
    def test_only_alascca_class_report_one_full_one_alascca_only(self, mock_get_hospital_code):
        mock_db_config_file = MagicMock()
        mock_get_hospital_code.side_effect = [353, 553]
        output = only_alascca_class_report("03098121", "03098849", mock_db_config_file)
    
        self.assertEquals(output, False)


    @patch('autoseq.util.report_type.get_hospital_code')
    def test_only_alascca_class_report_two_alascca_only(self, mock_get_hospital_code):
        mock_db_config_file = MagicMock()
        mock_get_hospital_code.side_effect = [553, 553]
        output = only_alascca_class_report("03098121", "03098849", mock_db_config_file)
    
        self.assertEquals(output, True)
