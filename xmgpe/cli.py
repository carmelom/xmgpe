#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Create: 08-2019 - Carmelo Mordini <carmelo> <carmelo.mordini@unitn.it>

"""Module docstring

"""
import click
import subprocess
import shutil
from pathlib import Path
from xmgpe.main import write_simulation
from xmgpe.movie import make_movie


@click.group()
def cli():
    """Command line interface
    """
    pass


@cli.command()
@click.argument('config')
@click.option('--compile', is_flag=True, help="Compile")
def write(config, compile=False):
    print(config)
    # scriptname = write_simulation(config)
    # if compile:
    #     subprocess.call(f"xmds2 {scriptname}", shell=True)


@cli.command()
@click.option('--filename', default='template.xmds')
def template(filename):
    source = Path(__file__).parent / 'templates' / 'template.xmds'
    dest = Path.cwd() / filename
    shutil.copy2(source, dest)


@cli.command()
@click.argument('Natoms', type=float)
@click.argument('trap-freqs', nargs=3, type=float)
@click.argument('ndims', type=int)
@click.option('--dx', default=0.33, type=float, help="lattice spacing (units of healing length, default: 0.33)")
@click.option('--box-size', default=2, type=float, help="simulation boundaries (units of TF radius, default: 2)")
@click.option('--name', default='groundstate', help="simulation filename (default: groundstate)")
@click.option('-n', '--dry-run', is_flag=True, default=False, help='Do not write the xmds script')
def setup(*args, **kwargs):
    write_simulation(*args, **kwargs)


@cli.command()
@click.argument('h5filename')
@click.option('--fps', default=20, type=float)
@click.option('--output', default='movie.mp4')
def movie(h5filename, fps, output):
    make_movie(h5filename, fps=fps, output=output)
