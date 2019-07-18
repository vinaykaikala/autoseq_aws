import os
import unittest

from autoseq.util.library import find_fastqs


class TestLibrary(unittest.TestCase):
    library = 'NA12877-N-03098121-TD1-TT1'
    libdir = 'tests/libraries'

    def test_find_fastqs_for_no_library(self):
        """
        test that find_fastqs return (None,None) if called with library=None
        """
        files = find_fastqs(library=None, libdir=self.libdir)
        self.assertEqual(files, (None, None))

    def test_find_fastq_gz(self):
        """
        test that files on the format *_1.fastq.gz / *_2.fastq.gz are found
        """
        files = find_fastqs(library=self.library, libdir=self.libdir)
        files_basenames = [os.path.basename(f) for f in files[0]] + [os.path.basename(f) for f in files[1]]
        self.assertIn('foo_1.fastq.gz', files_basenames)
        self.assertIn('foo_2.fastq.gz', files_basenames)

    def test_find_fq_gz(self):
        """
        test that files on the format *_1.fq.gz / *_2.fq.gz are found
        """
        files = find_fastqs(library=self.library, libdir=self.libdir)
        files_basenames = [os.path.basename(f) for f in files[0]] + [os.path.basename(f) for f in files[1]]
        self.assertIn('bar_1.fq.gz', files_basenames)
        self.assertIn('bar_2.fq.gz', files_basenames)

    def test_find_RN_DDD(self):
        """
        test that files on the format *R1_nnn.fastq.gz/*R2_nnn.fastq.gz are found
        """
        files = find_fastqs(library=self.library, libdir=self.libdir)
        files_basenames = [os.path.basename(f) for f in files[0]] + [os.path.basename(f) for f in files[1]]
        self.assertIn('baz_R1_001.fastq.gz', files_basenames)
        self.assertIn('baz_R2_001.fastq.gz', files_basenames)
        self.assertIn('baz_R1_999.fastq.gz', files_basenames)
        self.assertIn('baz_R2_999.fastq.gz', files_basenames)

    def test_fastq_path(self):
        """
        test that the components of the path to the files are in place
        """
        files = find_fastqs(library='NA12877-N-03098121-TD1-TT1', libdir='tests/libraries')
        self.assertIn(self.libdir, files[0][0])
        self.assertIn(self.library, files[0][0])
