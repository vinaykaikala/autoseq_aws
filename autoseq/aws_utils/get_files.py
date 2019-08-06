import json
import logging
import os
import sys
import time
import subprocess
import shlex
import pathlib
from autoseq.aws_utils.s3_files_config import files_for_each_step

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
        self.sample_data = {}
        self.files_for_each_step = files_for_each_step
        self.sdid = "" # this value is set when  ran set_fastq_files()

        #check and create directories for autoseq pipeline
        self.check_and_create_dir(self.base_dir)
        self.check_and_create_dir(self.outdir)
        self.check_and_create_dir(self.libdir)
        self.check_and_create_dir(self.refdata_dir)
        # check if output directory exist  if not create it

        #get ref file from s3
        self.get_s3files(refdata)

    def get_files_for_current_step(self, current_step):
        """
        Get required files for current step from s3
        :param current_step :
        """
        print(self.files_for_each_step)
        if current_step not in self.files_for_each_step:
            raise Exception('Invalid step name , configure required files in s3_files_config.py')

        base_dir = self.files_for_each_step[current_step]['base_dir']
        for each_file in self.files_for_each_step[current_step]['files']:
            file_complete_path = os.path.join(base_dir ,each_file['name'])
            if each_file['type'] == 'dir':
                self.get_s3directories(file_complete_path)
            else:
                self.get_s3files(file_complete_path)
        return True

    def get_s3files(self, *args):
        """Get common files required for all steps"""
        for each_file in args:
            if not os.path.exists(each_file):
                logging.info("Coping file from s3://{bucket}{filepath} to {filepath}".format(bucket=self.s3bucket, filepath=each_file))
                cmd ='aws s3 cp s3://{bucket}{file_path}  /{file_path}'.format(bucket=self.s3bucket, file_path=each_file)
                self.run_awscmd(cmd)
        return True

    def get_s3directories(self, dirname):
        """Get the files from s3 for given step"""
        cmd = 'aws s3 cp --recursive s3://{bucket}{filepath}  /{filepath}'.format(bucket=self.s3bucket, filepath=dirname)
        if not os.path.exists(dirname):
            logging.info(cmd)
            self.run_awscmd(cmd)
        return True

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
        #add conda env aws cli to run the commands
        cmd = '/usr/local/conda3/envs/awscli/bin/' + cmd
        subprocess.check_call(shlex.split(cmd))
        return True

    def set_fastq_files(self, sample_file):
        """
        :return: All clinseq barcodes included in this clinseq analysis pipeline's panel data.
        """
        self.sample_data = json.load(open(sample_file, 'r'))
        self.sdid = self.sample_data['sdid']
        all_clinseq_barcodes = \
            self.sample_data['T'] + \
            self.sample_data['N'] + \
            self.sample_data['CFDNA']
        filter(lambda bc: bc != None, all_clinseq_barcodes)

        temp_dirnames = []
        for each_fastq_dir in all_clinseq_barcodes:
            temp_dirnames.append({'name': each_fastq_dir, 'type': 'dir'})

        self.files_for_each_step['qc']['files'] = temp_dirnames

    def check_and_create_dir(self, dirname):
        """check if base path is same for current step"""
        if not os.path.join(*pathlib.Path(dirname).parts[0:3]) == self.base_dir:
            raise Exception('base path should be /nfs/PROBIO')
        if os.path.isfile(dirname):
            dirname = os.path.dirname(dirname)
        try:
            os.makedirs(dirname)
        except Exception as e:
            return dirname
        return dirname

