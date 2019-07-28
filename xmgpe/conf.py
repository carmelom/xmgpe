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
from lib_hf import hartree_fock as hf
from lib_hf.load_params import load_params
# Sodium constants
mass = 3.81924e-26
sclength = 2.75172e-09

parser = argparse.ArgumentParser()
parser.add_argument(
    '-p', '--params', default='params2d.yaml', help='Params file')
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

valid_solver = ['popov', 'hartree_fock']
if pars.solver not in valid_solver:
    raise ValueError(
        f"Invalid value for 'solver' param: must be one of {', '.join(valid_solver)} ({pars.solver})")

# -- Trap freqs are nice in the filename
nu_x = pars.geometry['x']['trap_freq']
nu_z = pars.geometry['z']['trap_freq']

results_file = f"gpe2d_{pars.solver}{tag}_{pars.mode}_{nu_x:.1f}_{nu_z:.1f}_{pars.natoms*1e-6:.1f}M_{pars.temperature*1e9:.0f}nK.h5"
results_file = Path(results_file)

print(f"results file {results_file} exists? {results_file.exists()}")
if results_file.exists() and not args.dry_run:
    if args.force:
        print(f"Removing {results_file}")
        results_file.unlink()
    else:
        raise FileExistsError(f"Results file {results_file} already exists.")

# -- Set globals
N = pars.natoms     # n atoms
T = pars.temperature  # temperature [K]
AR = nu_z / nu_x

omega_x = 2 * pi * nu_x
omega_z = 2 * pi * nu_z
a_z = sqrt(hbar / (mass * omega_z))
g1 = 4 * pi * sclength / a_z


hf_kwargs = dict(omega_ho=omega_x, dr=1e-7, Rmax=5e-3)


def n_hf(z, x, mu0, T, split=False):
    z, x = np.broadcast_arrays(z.reshape(-1, 1), x.reshape(1, -1))
    Vr = 0.5 * hbar * omega_z * (x**2/AR**2 + z**2)
    if split:
        n0, nt = hf.n_hartree_fock(
            mu0, T, Vr, split=split, solver_kwargs=hf_kwargs)
        n0 = np.clip(n0, a_min=0, a_max=None)
        return n0 * a_z**3, nt * a_z**3
    else:
        return hf.n_hartree_fock(mu0, T, Vr, split=split, solver_kwargs=hf_kwargs) * a_z**3


def integrate_N2d(z, x, n):
    n1 = 2*pi*np.trapz(n*z.reshape(-1, 1), x=z, axis=0)
    N = 2*np.trapz(n1, x=x)
    return N


def plot_hf(x, z, n0, nt):
    fig, (ax_x, ax_z) = plt.subplots(2,1)
    ax_x.plot(x, n0[0], '--')
    ax_x.plot(x, nt[0])
    ax_x.plot(x, n0[0] + nt[0])
    ax_x.set_xlabel('x')

    ax_z.plot(z, n0[:,0], '--')
    ax_z.plot(z, nt[:,0])
    ax_z.plot(z, n0[:,0] + nt[:,0])
    ax_z.set_xlabel('z')
    plt.show()


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

print(f" -- Initialize gpe with {N:.2e} atoms, T = {T*1e9:.2f} nK --")
mu0 = hf.get_mu0(N, T, mu0_lims=(-50, 200), omega_rho=omega_z, AR=AR) * 1e-9 * kB
q = '\n'.join(['simulation:', f'runtime: {runtime}',
               f'nsteps: {nsteps}', f'stepx: {stepx}',
               f'timestep : {timestep:.4e}',
               f'VN timest: {ts:.4e}'])
print(q)
print(f"HF mu0: {mu0*1e9/kB:.2f} nK ({mu0/hbar/omega_z:.2f} hbar omega)")

if args.dry_run:
    if args.plot:
        _x = np.linspace(0, pars.geometry['x']['box_size'], 200)
        _z = np.linspace(0, pars.geometry['z']['box_size'], 200)
        _n0, _nt = n_hf(_z, _x, mu0, T, split=True)
        plot_hf(_x, _z, _n0, _nt)
    sys.exit()

# compile
const = ['dt = %g' % timestep,
         'Uint = %g' % g1,
         'AR = %g' % AR,
         ]

sep = '\n\t'

cdata_str = sep.join([f'const real {s};' for s in const])

scriptname = pars.scriptname
scriptname_grid = 'grid_specifier.xmds'

# -- Prepare folder for run and storage
# -- Write xml
parser = etree.XMLParser(strip_cdata=False, remove_blank_text=True)
with open(scriptname, "r") as source:
    et = etree.parse(source, parser=parser)

with open(scriptname_grid, "r") as source:
    et_grid = etree.parse(source, parser=parser)


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
rewrite_geometry(et_grid)

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

# dest = 'test_configure.xmds'
dest = scriptname
print("Rewrite on %s...\n" % dest)
indent(et.getroot(), level=0)
et.write(dest, )

print('recompiling simulation')
subprocess.run('rm {0:s} {0:s}.cc'.format(
    os.path.splitext(dest)[0]), shell=True)
subprocess.run('xmds2 {:s}'.format(dest), shell=True)

dest = scriptname_grid
print("Rewrite geometry on %s...\n" % dest)
indent(et_grid.getroot(), level=0)
et_grid.write(dest, )

print('recompiling grid')
subprocess.run('rm {0:s} {0:s}.cc'.format(
    os.path.splitext(dest)[0]), shell=True)
subprocess.run('xmds2 {:s}'.format(dest), shell=True)

# -- initialisation
script_grid = os.path.splitext(scriptname_grid)[0]
subprocess.run(f"./{script_grid}")

with h5py.File('grid.h5') as f:
    x = f['1/x'][:]
    z = f['1/z'][:]
    shape = f['1/dummy'].shape


print('Creating new results file')
with h5py.File(results_file, 'w') as f:
    f['x'] = x
    f['z'] = z
    for key in 'mu', 'mu_used', 'energy', 'N0', 'Nt':
        f.create_dataset(key, shape=(0,), maxshape=(None,), dtype=float)
    for key in 'psi', 'nt':
        f.create_dataset(key, shape=(0,) + shape,
                         maxshape=(None,) + shape, dtype=float)
    for k, v in pars.my_dict.items():
        v = repr(v) if isinstance(v, (dict, type(None))) else v
        f.attrs[k] = v
    f.attrs['params_file'] = str(args.params)

n0, nt = n_hf(z, x, mu0, T, split=True)

N0 = integrate_N2d(z, x, n0)
Nt = integrate_N2d(z, x, nt)
mu0 = mu0 / hbar / omega_z
e0 = 5. / 7 * mu0 * N0  # Thomas-fermi result

# check that the number of atoms is consistent with LDA starting point
dN = (N0 + Nt - N) / N
print(f"N0: {N0:.3f} + Nt: {Nt:.3f} = {N0 + Nt:.3f} / {N:.3f} ({100*dN:.2f} %)")

results = {'mu': mu0,
           'mu_used': mu0,
           'energy': e0,
           'N0': N0,
           'Nt': Nt,
           'psi': np.sqrt(n0),
           'nt': nt
           }
# save results
print('save results')
with h5py.File(results_file) as f:
    for name, data in results.items():
        dset = f[name]
        dset.resize(dset.shape[0] + 1, axis=0)
        print(dset)
        dset[-1] = data

# some optional plots
if args.plot:
    plot_hf(x, z, n0, nt)
