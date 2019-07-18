import unittest
from mock import patch, mock_open
from autoseq.cli.cli import *


class TestCLI(unittest.TestCase):
    def test_convert_to_absolute_path_non_string1(self):
        self.assertEquals(
            convert_to_absolute_path(None, "/dummy/base/dir"),
            None
        )

    def test_convert_to_absolute_path_non_string2(self):
        self.assertEquals(
            convert_to_absolute_path(1, "/dummy/base/dir"),
            1
        )

    def test_convert_to_absolute_path_absolute_path(self):
        self.assertEquals(
            convert_to_absolute_path("/an/absolute/path", "/dummy/base/dir"),
            "/an/absolute/path"
        )

    @patch('autoseq.cli.cli.os.path.isfile')
    def test_convert_to_absolute_path_relative_path_file(self, mock_isfile):
        mock_isfile.return_value = True
        self.assertEquals(
            convert_to_absolute_path("a_terminal_filename", "/dummy/base/dir"),
            "/dummy/base/dir/a_terminal_filename"
        )

    @patch('autoseq.cli.cli.os.path.isdir')
    def test_convert_to_absolute_path_relative_path_dir(self, mock_isdir):
        mock_isdir.return_value = True
        self.assertEquals(
            convert_to_absolute_path("a_terminal_dirname", "/dummy/base/dir"),
            "/dummy/base/dir/a_terminal_dirname"
        )

    @patch('autoseq.cli.cli.os.path.isfile')
    @patch('autoseq.cli.cli.os.path.isdir')
    def test_convert_to_absolute_path_relative_not_there(self, mock_isdir, mock_isfile):
        mock_isdir.return_value = False
        mock_isfile.return_value = False
        self.assertEquals(
            convert_to_absolute_path("a_relative_filename_not_existing", "/dummy/base/dir"),
            "a_relative_filename_not_existing"
        )

    @patch('autoseq.cli.cli.os.path.isfile')
    def test_make_paths_absolute(self, mock_isfile):
        mock_isfile.return_value = True
        input_dict = {
            "some_key1": None,
            "some_key2": 1,
            "some_key3": "a_terminal_filename",
            "nested_dict": {"some_key4": "filename2"}
        }
        output_dict = make_paths_absolute(input_dict, "/dummy/base/dir")
        self.assertEquals(output_dict["some_key1"], None)
        self.assertEquals(output_dict["some_key2"], 1)
        self.assertEquals(output_dict["some_key3"], "/dummy/base/dir/a_terminal_filename")
        self.assertEquals(output_dict["nested_dict"]["some_key4"], "/dummy/base/dir/filename2")

    @patch('autoseq.cli.cli.os.path.isfile')
    def test_make_paths_absolute(self, mock_isfile):
        mock_isfile.return_value = False
        input_dict = {
            "some_key1": None,
            "some_key2": 1,
            "some_key3": "a_terminal_filename",
            "nested_dict": {"some_key4": "filename2"}
        }
        output_dict = make_paths_absolute(input_dict, "/dummy/base/dir")
        self.assertEquals(output_dict["some_key1"], None)
        self.assertEquals(output_dict["some_key2"], 1)
        self.assertEquals(output_dict["some_key3"], "a_terminal_filename")
        self.assertEquals(output_dict["nested_dict"]["some_key4"], "filename2")

    @patch('autoseq.cli.cli.os.path.isfile')
    def test_load_ref(self, mock_isfile):
        mock_isfile.return_value = True
        mocked_open = mock_open(read_data='{"some_key": "a_terminal_filename"}')

        with patch('autoseq.cli.cli.open', mocked_open, create=True):
            loaded_ref = load_ref("/dummy/base/dir/dummy_file.json")
            self.assertEquals(loaded_ref["some_key"], "/dummy/base/dir/a_terminal_filename")
