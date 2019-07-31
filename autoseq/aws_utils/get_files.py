import json
import logging
import os
import sys
import time
import subprocess
import shlex
import pathlib

class Awscli():
    """Get required files from s3 to run the pipeline"""
    def __init__(self, sampledata, refdata, outdir, libdir, s3bucket='probio-genome'):

        """
        :param sampledata:
        :param refdata:
        :param outdir:
        :param libdir:
        :param s3bucket:

        the base directory path should be
        base: /nfs/PROBIO
        Fastq : /nfs/PROBIO/INBOX/
        outputs: /nfs/PROBIO/autoseq-output/<sdid>
        confs: /nfs/PROBIO/config/<sdid>.json
        """

        self.base_dir = '/nfs/PROBIO'
        self.s3bucket = s3bucket
        self.outdir = outdir
        #self.sampledata = json.load(sampledata)
        self.sampledata_dir = os.path.dirname(self.sampledata)
        self.refdata = refdata
        self.refdata_dir = os.path.dirname(refdata)
        self.libdir = libdir

        #check and create directories for autoseq pipeline
        self.check_and_create_dir(self, self.base_dir)
        self.check_and_create_dir(self, self.outdir)
        self.check_and_create_dir(self, self.libdir)
        self.check_and_create_dir(self, self.refdata_dir)
        self.check_and_create_dir(self, self.sampledata_dir)
        # check if output directory exist  if not create it



    def get_common(self):
        """Get common files required for all steps"""
        cmd = "aws s3 ls s3://" + self.s3bucket
        self.run_awscmd(cmd)

    def get(self, step):
        """Get the files from s3 for given step"""
        pass

    def put(self):
        pass

    def download_folder(s3_path, directory_to_download):
        """
        Downloads a folder from s3
        :param s3_path: s3 folder path
        :param directory_to_download: path to download the directory to
        :return: directory that was downloaded
        """
        cmd = 'aws s3 cp --recursive %s %s' % (s3_path, directory_to_download)

    def run_awscmd(self, cmd):
        #add conda env
        cmd = 'source activate awscli && ' + cmd
        subprocess.check_call(shlex.split(cmd))
        return True

    def get_all_clinseq_barcodes(self):
        """
        :return: All clinseq barcodes included in this clinseq analysis pipeline's panel data.
        """
        all_clinseq_barcodes = \
            self.sampledata['T'] + \
            self.sampledata['N'] + \
            self.sampledata['CFDNA']
        return filter(lambda bc: bc != None, all_clinseq_barcodes)

    def check_and_create_dir(self, dirname):
        """check if base path is same for current step"""
        if not os.path.join(*pathlib.Path(dirname).parts[0:3]) == self.base_dir:
            raise Exception('base path should be /nfs/PROBIO')

        try:
            os.makedirs(dirname)
        except Exception as e:
            return dirname
        return dirname

