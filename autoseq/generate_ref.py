import logging
import sys
import click

from autoseq.cli.cli import setup_logging, get_runner
from autoseq.util.path import mkdir
from autoseq.pipeline.generate_ref_files_pipeline import GenerateRefFilesPipeline

__author__ = 'dankle'


@click.command()
@click.option('--genome-resources', help='reference sequence', type=str)
@click.option('--outdir', default='/tmp/autoseq-test', help='output directory', type=click.Path())
@click.option('--runner_name', default='shellrunner', help='Runner to use.')
@click.option('--loglevel', default='INFO', help='level of logging')
@click.option('--cores', default=1, help="write graph to dot file with this name")
@click.option('--debug', default=False, is_flag=True)
def main(genome_resources, outdir, runner_name, loglevel, cores, debug):
    setup_logging(loglevel)

    mkdir(outdir)
    logging.info("Writing to {}".format(outdir))

    runner = get_runner(runner_name, cores)
    p = GenerateRefFilesPipeline(genome_resources, outdir,
                                 maxcores=cores, runner=runner)

    # start main analysis
    p.start()

    # block thread until done
    p.join()

    # return with exitcode
    sys.exit(p.exitcode)

if __name__ == "__main__":
    sys.exit(main())
