import logging
import os

__author__ = 'dankle'

class Readpair(object):
    """
    Container for a pair of fastq files
    """
    def __init__(self, library, fq1, fq2):
        self.LIBRARY = library
        self.FQ1 = fq1
        self.FQ2 = fq2
        logging.debug("Instantiated Readpair from library {l}".format(l=self.LIBRARY))


    def __str__(self):
        return ",".join([self.LIBRARY,
                          self.FQ1,
                          self.FQ2])

    @staticmethod
    def fromDir(basedir):
        """Get a list of read pairs from a dir"""
        dirs = [d for d in os.listdir(basedir) if os.path.isdir(basedir+"/"+d)]
        logging.debug("dirs to consider = {}".format(dirs))
        readpairs = []
        for d in dirs:
            files = os.listdir(basedir+"/"+d)
            fq1 = [basedir+"/"+d+"/"+f for f in files if f.endswith("_1.fastq.gz")]
            fq2 = [basedir+"/"+d+"/"+f for f in files if f.endswith("_2.fastq.gz")]
            if fq1 is not []:
                if fq2 == []:
                    fq2 = None
                readpairs.append(Readpair(library=d, fq1=fq1, fq2=fq2))
 
        return readpairs
