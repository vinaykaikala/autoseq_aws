import unittest
from autoseq.tools.msi import *


class TestMsi(unittest.TestCase):
    def test_msi_sensor_scan(self):
        msi_sensor_scan = MsiSensorScan()
        msi_sensor_scan.input_fasta = "input.fa"
        msi_sensor_scan.output = "test_output"
        cmd = msi_sensor_scan.command()
        self.assertIn('input.fa', cmd)
        self.assertIn('test_output', cmd)
    
    
    def test_intersect_msi_sites(self):
        msi_sensor_scan = IntersectMsiSites()
        msi_sensor_scan.input_msi_sites = "test_input"
        msi_sensor_scan.target_bed = "test.bed"
        msi_sensor_scan.output_msi_sites = "test_output"
        cmd = msi_sensor_scan.command()
        self.assertIn('test_input', cmd)
        self.assertIn('test.bed', cmd)
        self.assertIn('test_output', cmd)
    
    
    def test_msi_sensor(self):
        msi_sensor = MsiSensor()
        msi_sensor.msi_sites = "dummy_msi_sites"
        msi_sensor.input_normal_bam = "normal.bam"
        msi_sensor.input_tumor_bam = "tumor.bam"
        msi_sensor.output = "test_output"
        cmd = msi_sensor.command()
        self.assertIn('dummy_msi_sites', cmd)
        self.assertIn('normal.bam', cmd)
        self.assertIn('tumor.bam', cmd)
        self.assertIn('test_output', cmd)
    
    
    def test_msings(self):
        msings = Msings()
        msings.input_fasta = "dummy.fasta"
        msings.outdir = "test_output"
        msings.msings_baseline = "dummy.baseline"
        msings.msings_bed = "dummy.bed"
        msings.msings_intervals = "dummy.intervals"
        msings.input_bam = "input.bam"
        cmd = msings.command()
        self.assertIn("dummy.fasta", cmd)
        self.assertIn("test_output", cmd)
        self.assertIn("dummy.baseline", cmd)
        self.assertIn("dummy.bed", cmd)
        self.assertIn("dummy.intervals", cmd)
        self.assertIn("input.bam", cmd)