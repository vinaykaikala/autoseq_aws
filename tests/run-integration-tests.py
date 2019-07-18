import logging
import os
import shutil
import subprocess
import sys

import click

from autoseq.tests import alascca_test_outdir, liqbio_test_outdir

test_files = ['integration/test_alascca.py', 'integration/test_liqbio.py']


@click.command()
@click.option('--force-download', is_flag=True, help="Download test libs and ref even if in place")
def run_tests(force_download):
    """Run the autoseq integration tests"""

    logging.basicConfig(level=logging.INFO)

    download_test_libraries(force_download)
    download_test_genome(force_download)

    if os.path.exists(alascca_test_outdir) or os.path.exists(liqbio_test_outdir):
        if query_yes_no("Output path exists, delete old files? (Use 'N' to continue with the present state)"):
            for p in [alascca_test_outdir, liqbio_test_outdir]:
                if os.path.exists(p):
                    shutil.rmtree(p)

    for f in test_files:
        subprocess.check_call(['py.test', '--capture=no', '--ignore=tests/integration', 'tests/{}'.format(f)])


def download_test_libraries(force_download):
    """
    Download test libraries from uppmax
    """
    url = "https://export.uppmax.uu.se/b2010040/test-libraries.tar.gz"
    local_f = "/tmp/test-libraries.tar.gz"

    if force_download and os.path.exists(local_f):
        os.remove(local_f)

    if not os.path.exists(local_f):
        logging.info("Downloading and unpacking test libraries")
        subprocess.check_call(["wget", "--continue", "--no-clobber", "-O", local_f, url, "-q"])
        subprocess.check_call(["tar", "xvfz", local_f, "-C", "/tmp"])


def download_test_genome(force_download):
    """
    Download test genome from uppmax
    """
    url = "https://export.uppmax.uu.se/b2010040/test-genome.tar.gz"
    local_f = "/tmp/test-genome.tar.gz"

    if force_download and os.path.exists(local_f):
        os.remove(local_f)

    if not os.path.exists(local_f):
        logging.info("Downloading and unpacking test genome")
        subprocess.check_call(["wget", "--continue", "--no-clobber", "-O", local_f, url, "-q"])
        subprocess.check_call(["tar", "xvfz", local_f, "-C", "/tmp"])


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".

    from http://stackoverflow.com/questions/3041986/apt-command-line-interface-like-yes-no-input
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")


if __name__ == '__main__':
    run_tests()
