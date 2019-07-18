import json
import logging
import datetime

import click

from autoseq.cli.cli import setup_logging


@click.command()
@click.option('--jobdb', required=True, help='json with reference files to use', type=str)
@click.option('--svg', required=True, help="output svg file", type=str)
@click.option('--loglevel', default='INFO', help='level of logging')
def cli(jobdb, svg, loglevel):
    setup_logging(loglevel)

    jdb = json.load(open(jobdb, 'r'))
    jobs = jdb['jobs']

    import matplotlib.pyplot as plt
    plt.switch_backend('agg')
    import matplotlib.patches as patches

    x1 = [datetime_to_timestamp(deserialize_date(j['starttime'])) for j in jobs if j['starttime'] and j['endtime']]
    x2 = [datetime_to_timestamp(deserialize_date(j['endtime'])) for j in jobs if j['starttime'] and j['endtime']]
    y = xrange(len([1 for j in jobs if j['starttime'] and j['endtime']]))
    status = [j['status'] for j in jobs if j['starttime'] and j['endtime']]
    jobnames = [j['jobname'] for j in jobs if j['starttime'] and j['endtime']]
    color_mapper = {'FAILED': 'red',
                    'CANCELLED': 'pink',
                    'RUNNING': 'blue',
                    'COMPLETED': '#00ff00',
                    'PENDING': 'lightblue'}

    fig6 = plt.figure()
    plt.axis([0, 1.01*(max(x2)-min(x1)), 0, len(x1)])
    ax6 = fig6.add_subplot(111)


    for idx, j in enumerate([jb for jb in jobs if jb['starttime'] and jb['endtime']]):
        tstart = x1[idx] - min(x1)
        tend = x2[idx] - min(x1)
        p = patches.Rectangle(
            (tstart, idx),
            tend - tstart,
            0.9,
            facecolor=color_mapper[status[idx]]
        )
        ax6.add_patch(p)
        plt.text(tstart+20, idx+.25, jobnames[idx], fontdict={'size': '5'})

    # for p in [
    #     patches.Rectangle(
    #         (0.03, 0.1), 0.2, 0.6,
    #         facecolor=None  # Default
    #     ),
    #     patches.Rectangle(
    #         (0.26, 0.1), 0.2, 0.6,
    #         facecolor="none"  # No background
    #     ),
    #     patches.Rectangle(
    #         (0.49, 0.1), 0.2, 0.6,
    #         facecolor="red"
    #     ),
    #     patches.Rectangle(
    #         (0.72, 0.1), 0.2, 0.6,
    #         facecolor="#00ffff"
    #     ),
    # ]:
    #     ax6.add_patch(p)
    fig6.savefig(svg)



    # import numpy as np
    # import matplotlib.pyplot as plt
    #
    # # Read data from file into variables
    # y, c, x1, x2 = np.loadtxt('data.txt', unpack=True)
    #
    # # Map value to color
    # color_mapper = np.vectorize(lambda x: {0: 'red', 1: 'blue'}.get(x))
    #
    # # Plot a line for every line of data in your file
    # plt.hlines(y, x1, x2, colors=color_mapper(c))

    # import numpy as np
    # import matplotlib.pyplot as plt
    # plt.switch_backend('agg')
    #
    # x1 = [deserialize_date(j['starttime']) for j in jobs if j['starttime'] and j['endtime']]
    # x2 = [deserialize_date(j['endtime']) for j in jobs if j['starttime'] and j['endtime']]
    # y = xrange(len([1 for j in jobs if j['starttime'] and j['endtime']]))
    # status = [j['status'] for j in jobs if j['starttime'] and j['endtime']]
    #
    # color_mapper = np.vectorize(lambda x: {'FAILED': 'red',
    #                                        'CANCELLED': 'pink',
    #                                        'RUNNING': 'blue',
    #                                        'COMPLETED': 'green',
    #                                        'PENDING': 'lightblue'}.get(x))
    #
    # plt.hlines(y, x1, x2, colors=color_mapper(status))
    # plt.savefig(svg)
    #

    # Change font default
    # gantt.define_font_attributes(fill='black',
    #                              stroke='black',
    #                              stroke_width=0,
    #                              font_family="Verdana")
    #
    # res = gantt.Resource('Autoseq jobs')
    # t1 = gantt.Task(name='tache1',
    #                 start=datetime.datetime(2014, 12, 25, 12, 00, 00),
    #                 duration=4,
    #                 percent_done=44,
    #                 resources=[res],
    #                 color="#FF8080")
    #
    # # Create a project
    # p1 = gantt.Project(name='Projet 1')
    #
    # # Add tasks to this project
    # p1.add_task(t1)
    #
    # p = gantt.Project(name='Gantt')
    # # wich contains the first two projects
    # # and a single task
    # p.add_task(p1)
    #
    # p.make_svg_for_tasks(filename=svg,
    #                      today=datetime.datetime(2014, 12, 26),
    #                      start=datetime.datetime(2014, 12, 22),
    #                      end=datetime.datetime(2015, 01, 1))


def deserialize_date(d):
    obj = None
    try:
        obj = datetime.datetime.strptime(d, "%Y-%m-%dT%H:%M:%S.%f")
    except ValueError:
        obj = datetime.datetime.strptime(d, "%Y-%m-%dT%H:%M:%S")
    return obj


def datetime_to_timestamp(d):
    return (d - datetime.datetime(1970, 1, 1)).total_seconds()
