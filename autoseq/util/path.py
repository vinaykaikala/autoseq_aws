import logging
import os
import subprocess


def normpath(path):
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))


def stripsuffix(thestring, suffix):
    """
    Strip the given suffix from the given string
    :param thestring: string from which to strip
    :param suffix: suffix to strip
    :return: string without the given suffix. If the given string
    doesn't end in suffix, the original string is returned
    """
    if thestring.endswith(suffix):
        return thestring[:-len(suffix)]
    return thestring


def mkdir(dir):
    """ Create a directory if it doesn't exist
    :param dir: dir to create
    """
    if not os.path.isdir(dir):
        try:
            os.makedirs(dir)
        except OSError:
            logging.error("Couldn't create directory {}".format(dir))
    else:  # if dir already exists, do nothing
        pass


def transfer_data_to_tmpdir(sampledata, tmpdir, final_outdir):
    rawdata_dir = "{}/raw/".format(tmpdir)
    new_outdir = "{}/{}".format(tmpdir, sampledata['REPORTID'])
    mkdir(rawdata_dir)
    mkdir(new_outdir)

    logging.debug("Syncing final output dir")
    rsync_dirs(final_outdir, new_outdir)

    logging.debug("Fetching raw data")
    new_sampledata = fetch_raw_data(sampledata, rawdata_dir)
    return new_sampledata


def rsync_dirs(src, target):
    """
    Sync two a source dir to target
    :param src:
    :param target:
    :return:
    """
    src = os.path.expandvars(os.path.expanduser(src))
    target = os.path.expandvars(os.path.expanduser(target))
    rsync_results_cmd = ['rsync', '--stats', '--size-only', '--delete', '-avP', src + "/", target + "/"]
    logging.debug("Running rsync with command")
    logging.debug("{}".format(" ".join(rsync_results_cmd)))
    subprocess.check_call(rsync_results_cmd, stderr=open('/dev/null', 'w'), stdout=open('/dev/null', 'w'))


def rsync_file(src, target):
    """
    Sync src file to target
    :param src:
    :param target:
    """
    src = os.path.expandvars(os.path.expanduser(src))
    target = os.path.expandvars(os.path.expanduser(target))
    logging.info("Copying {} to local work dir".format(src))
    rsync_command = ["rsync", "--stats", "--size-only", src, target]
    mkdir(os.path.dirname(target))
    logging.debug("Executing {}".format(rsync_command))
    subprocess.check_call(rsync_command, stderr=open('/dev/null', 'w'), stdout=open('/dev/null', 'w'))


def fetch_raw_data(sampledata, rawdata_dir):
    """
    Copy data for a single report to a directory
    """
    items_to_copy = ["PANEL_TUMOR_FQ1", "PANEL_TUMOR_FQ2", "PANEL_NORMAL_FQ1", "PANEL_NORMAL_FQ2",
                     "WGS_TUMOR_FQ1", "WGS_TUMOR_FQ2", "WGS_NORMAL_FQ1", "WGS_NORMAL_FQ2",
                     "RNASEQ_FQ1", "RNASEQ_FQ2", "RNASEQCAP_FQ1", "RNASEQCAP_FQ2"]

    for item in items_to_copy:
        if sampledata[item] is not None and sampledata[item] is not "NA":
            new_fqs = []
            for f in sampledata[item]:
                local_file = os.path.abspath("{}/{}".format(rawdata_dir, os.path.expanduser(f)))
                new_fqs.append(local_file)
                rsync_file(f, local_file)
            sampledata[item] = new_fqs
    return sampledata