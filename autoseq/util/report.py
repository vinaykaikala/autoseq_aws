import logging

__author__ = 'dankle'

def collapse(arr):
    """
    arr is an array of array like so
    arr = [[1,2], [3], [4,5,6]]
    collapse(arr) gives
    [1,2,3,4,5,6]
    """
    cc = []
    for elm in arr:
        if elm is not None:
            cc = cc + elm
    if cc == []:
        return None
    return cc


class Report(object):
    """
    Container for a single report
    """

    def __init__(self, line):
        elements = line.strip("\n").split("\t")
        elements = map(lambda x: None if x == "NA" or x == "" else x, elements)

        self.REPORTID = elements[0]
        self.PATIENTID = elements[1]
        self.TARGETS = elements[2]
        self.SUBSTUDY = elements[3]
        self.PANEL_TUMOR_LIB = elements[4]
        self.PANEL_NORMAL_LIB = elements[5]
        self.WGS_TUMOR_LIB = elements[6]
        self.WGS_NORMAL_LIB = elements[7]
        self.RNASEQ_LIB = elements[8]
        self.RNASEQCAP_LIB = elements[9]
        logging.debug("Instantiated Report with id {r}".format(r=self.REPORTID))

    def __str__(self):
        return "\t".join([self.REPORTID,
                          self.PATIENTID,
                          self.TARGETS,
                          self.SUBSTUDY,
                          self.PANEL_TUMOR_LIB,
                          self.PANEL_NORMAL_LIB,
                          self.WGS_TUMOR_LIB,
                          self.WGS_NORMAL_LIB,
                          self.RNASEQ_LIB,
                          self.RNASEQCAP_LIB])

    def to_dict(self, readpairs):
        d = {'REPORTID': self.REPORTID,
             'PATIENTID': self.PATIENTID,
             'TARGETS': self.TARGETS,
             'SUBSTUDY': self.SUBSTUDY,
             'PANEL_TUMOR_LIB': self.PANEL_TUMOR_LIB,
             'PANEL_NORMAL_LIB': self.PANEL_NORMAL_LIB,
             'WGS_TUMOR_LIB': self.WGS_TUMOR_LIB,
             'WGS_NORMAL_LIB': self.WGS_NORMAL_LIB,
             'RNASEQ_LIB': self.RNASEQ_LIB,
             'RNASEQCAP_LIB': self.RNASEQCAP_LIB}

        libraries = [self.PANEL_NORMAL_LIB, self.PANEL_TUMOR_LIB,
                     self.WGS_NORMAL_LIB, self.WGS_TUMOR_LIB,
                     self.RNASEQ_LIB, self.RNASEQCAP_LIB]

        parts = ['PANEL_TUMOR', 'PANEL_NORMAL', 'WGS_TUMOR', 'WGS_NORMAL', 'RNASEQ', 'RNASEQCAP']

        d['PANEL_TUMOR_FQ1'] = collapse([r.FQ1 for r in readpairs if r.LIBRARY == self.PANEL_TUMOR_LIB])
        d['PANEL_TUMOR_FQ2'] = collapse([r.FQ2 for r in readpairs if r.LIBRARY == self.PANEL_TUMOR_LIB])
        d['PANEL_NORMAL_FQ1'] = collapse([r.FQ1 for r in readpairs if r.LIBRARY == self.PANEL_NORMAL_LIB])
        d['PANEL_NORMAL_FQ2'] = collapse([r.FQ2 for r in readpairs if r.LIBRARY == self.PANEL_NORMAL_LIB])
        d['WGS_TUMOR_FQ1'] = collapse([r.FQ1 for r in readpairs if r.LIBRARY == self.WGS_TUMOR_LIB])
        d['WGS_TUMOR_FQ2'] = collapse([r.FQ2 for r in readpairs if r.LIBRARY == self.WGS_TUMOR_LIB])
        d['WGS_NORMAL_FQ1'] = collapse([r.FQ1 for r in readpairs if r.LIBRARY == self.WGS_NORMAL_LIB])
        d['WGS_NORMAL_FQ2'] = collapse([r.FQ2 for r in readpairs if r.LIBRARY == self.WGS_NORMAL_LIB])
        d['RNASEQ_FQ1'] = collapse([r.FQ1 for r in readpairs if r.LIBRARY == self.RNASEQ_LIB])
        d['RNASEQ_FQ2'] = collapse([r.FQ2 for r in readpairs if r.LIBRARY == self.RNASEQ_LIB])
        d['RNASEQCAP_FQ1'] = collapse([r.FQ1 for r in readpairs if r.LIBRARY == self.RNASEQCAP_LIB])
        d['RNASEQCAP_FQ2'] = collapse([r.FQ2 for r in readpairs if r.LIBRARY == self.RNASEQCAP_LIB])

        return d

    @staticmethod
    def fromFile(f):
        """
        :param f: get a list of reports read from a reports file
        :return: list of Reports
        """
        reportsf = open(f, 'r')
        reports = [Report(line) for line in reportsf if line.strip() != "" and not line.strip().startswith("#")]
        reportsf.close()
        return reports

