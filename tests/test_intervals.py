import unittest
from autoseq.tools.intervals import *


class TestIntervals(unittest.TestCase):
    def test_slop_interval_list(self):
        slop_interval_list = SlopIntervalList()
        slop_interval_list.input = "test_input"
        slop_interval_list.output = "test_output"
        cmd = slop_interval_list.command()
        self.assertIn('test_input', cmd)
        self.assertIn('test_output', cmd)

    def test_interval_list_to_bed(self):
        interval_list_to_bed = IntervalListToBed()
        interval_list_to_bed.input = "test_input"
        interval_list_to_bed.output = "test_output"
        cmd = interval_list_to_bed.command()
        self.assertIn('test_input', cmd)
        self.assertIn('test_output', cmd)