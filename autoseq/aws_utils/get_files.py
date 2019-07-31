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
    def __init__(self,  refdata, outdir, libdir, s3bucket='probio-genome'):

        """
        :param refdata:
        :param outdir:
        :param libdir:
        :param s3bucket:

        the base directory path should be
        base: /nfs/PROBIO
        Fastq : /nfs/PROBIO/INBOX/                    --libdir
        outputs: /nfs/PROBIO/autoseq-output/<sdid>    --outdir
        confs: /nfs/PROBIO/config/<sdid>.json         --sample
        refdata: /nfs/PROBIO/autoseq-genome           --ref
        """

        self.base_dir = '/nfs/PROBIO'
        self.s3bucket = s3bucket
        self.outdir = outdir
        self.refdata = refdata
        self.refdata_dir = os.path.dirname(refdata)
        self.libdir = libdir

        #check and create directories for autoseq pipeline
        self.check_and_create_dir(self.base_dir)
        self.check_and_create_dir(self.outdir)
        self.check_and_create_dir(self.libdir)
        self.check_and_create_dir(self.refdata_dir)
        self.check_and_create_dir(self.sampledata_dir)
        # check if output directory exist  if not create it

        #get ref file from s3
        self.get_s3files(refdata)




    def get_s3files(self, *args):
        """Get common files required for all steps"""

        for each_file in args:
            if not os.path.exists(each_file):
                logging.info("Coping file from s3 to {}".format(each_file))
                cmd ='aws s3 cp s3://{bucket}/{file_path}  {file_path}'.format(bucket=self.s3bucket, file_path=each_file)
                self.run_awscmd(cmd)

        return True

    def get_s3directories(self, step):
        """Get the files from s3 for given step"""
        pass

    def put_file_to_s3(self):
        pass

    def put_directories_s3(self):
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

