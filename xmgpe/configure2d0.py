#!/usr/bin/env python
import os
import sys
import subprocess
from pathlib import Path
from math import pi, sqrt
from lxml import etree

import argparse
import numpy as np
import matplotlib.pyplot as plt
import h5py

# Physics
from scipy.constants import hbar, Boltzmann as kB
# from lib_hf import hartree_fock as hf
from lib_hf.load_params import load_params
# Sodium constants
mass = 3.81924e-26
sclength = 2.75172e-09

parser = argparse.ArgumentParser()
parser.add_argument(
    '-p', '--params', default='params2d0.yaml', help='Params file')
parser.add_argument('-t', '--tag', help='Tag name')
parser.add_argument('--plot', action='store_true', help='Final plot')
group = parser.add_mutually_exclusive_group()
group.add_argument('-n', '--dry-run', action='store_true', help='Dry run')
group.add_argument('-f', '--force', action='store_true',
                   help='Overwrite existing file')
args = parser.parse_args()

print(f"Loading {args.params}")
_, pars = load_params(args.params, bunch=True)

tag = '_' + args.tag if args.tag is not None else ''
if args.tag:
    pars.my_dict['tag'] = tag


# -- Set globals
N = pars.natoms     # n atoms
nu_x = pars.geometry['x']['trap_freq']
nu_z = pars.geometry['z']['trap_freq']
AR = nu_z / nu_x


results_file = f"gpe2d_T0{tag}_results_{pars.natoms*1e-6:.1f}M_AR_{AR:.1f}.h5"
results_file = Path(results_file)

print(f"results file {results_file} exists? {results_file.exists()}")
if results_file.exists() and not args.dry_run:
    if args.force:
        print(f"Removing {results_file}")
        results_file.unlink()
    else:
        raise FileExistsError(f"Results file {results_file} already exists.")

omega_x = 2 * pi * nu_x
omega_z = 2 * pi * nu_z
a_z = sqrt(hbar / (mass * omega_z))
g1 = 4 * pi * sclength / a_z


# -- Set geometry
# dimension_name = 'x'
stepx = [dk['box_size'] / dk['npoints'] for dk in pars.geometry.values()]
# stepx = [print(dk['box_size'],  dk['npoints']) for dk in pars.geometry.keys()]
# print(pars.geometry)

# -- set optimal fixed timestep
runtime = pars.runtime
timestep = pars.timestep

sample_timestep = 2 * [pars.sample_timestep]

stepx = min(stepx)
ts = 0.25 * stepx**2  # Heuristic based on Von Neumann criterion
nsteps = int(runtime / timestep)
timestep = float(runtime) / nsteps  # the actual one

mu0_tilde = 0.5*(15*N*g1/4/pi/AR)**(2./5)
mu0 = mu0_tilde*hbar*omega_z

print(f" -- Initialize gpe with {N:.2e} atoms--")
q = '\n'.join(['simulation:', f'runtime: {runtime}',
               f'nsteps: {nsteps}', f'stepx: {stepx}',
               f'timestep : {timestep:.4e}',
               f'VN timest: {ts:.4e}'])
print(q)
print(f"HF mu0: {mu0*1e9/kB:.2f} nK ({mu0/hbar/omega_z:.2f} hbar omega)")

if args.dry_run:
    sys.exit()

# compile
const = [f'dt = {timestep:g}',
         f'Uint = {g1:g}',
         f'AR = {AR:g}',
         f'Nparticles = {N:e}',
         'mu0 = 0.5*pow(15*Nparticles*Uint/4/M_PI/AR, 2.0/5.0)'
         ]

sep = '\n\t'

cdata_str = sep.join([f'const real {s};' for s in const])

scriptname = pars.scriptname

# -- Prepare folder for run and storage
# -- Write xml
parser = etree.XMLParser(strip_cdata=False, remove_blank_text=True)
with open(scriptname, "r") as source:
    et = etree.parse(source, parser=parser)


def indent(elem, level=0):
    i = "\n" + level * "  "
    i2 = "\n\n" if level < 2 else "\n"
    i2 += level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i2
        for elem in elem:
            indent(elem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i2
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i2


name = et.find('name')
name.text = os.path.splitext(pars.scriptname)[0]

glob = et.find('features').find('globals')
glob.text = etree.CDATA(sep + cdata_str + sep)
# print(etree.tostring(glob, encoding='unicode', pretty_print=True))


def rewrite_geometry(et):
    geom = et.find('geometry')
    td = etree.Element('transverse_dimensions')
    for dimension_name, values in pars.geometry.items():
        dim = etree.Element('dimension')
        dim.attrib['name'] = dimension_name
        dim.attrib['lattice'] = str(values['npoints'])
        dim.attrib['domain'] = f"(0, {values['box_size']:.3f})"
        dim.attrib['transform'] = values['transform']
        dim.attrib['volume_prefactor'] = values['prefactor']
        td.append(dim)
    geom.replace(geom.find('transverse_dimensions'), td)
    # print(etree.tostring(geom, encoding='unicode', pretty_print=True))


rewrite_geometry(et)

sequence = et.find('sequence')
# int_real = list(sequence.iterfind('integrate'))[-1]
# samples = int_real.find('samples')
#
# nsamples = [str(int(args.runtime / timestep)) for timestep in args.sample_timestep]
# int_real.attrib['interval'] = str(args.runtime)
# samples.text = ' '.join(nsamples)
# print samples.text

integrate = sequence.find('integrate')
print('set runtime = %.1f, steps = %d for integrate element' % (runtime, nsteps))
samples = integrate.find('samples')
nsamples = [str(int(runtime / sts)) for sts in sample_timestep]
integrate.attrib['interval'] = str(runtime)
integrate.attrib['steps'] = str(nsteps)
samples.text = ' '.join(nsamples)
print(samples.text)

output = et.find('output')
output.attrib['filename'] = str(results_file.stem)

dest = scriptname
print("Rewrite on %s...\n" % dest)
indent(et.getroot(), level=0)
et.write(dest, )

print('recompiling simulation')
subprocess.run('rm {0:s} {0:s}.cc'.format(
    os.path.splitext(dest)[0]), shell=True)
subprocess.run('xmds2 {:s}'.format(dest), shell=True)
