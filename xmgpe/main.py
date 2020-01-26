#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Create: 08-2019 - Carmelo Mordini <carmelo> <carmelo.mordini@unitn.it>

"""Module docstring

"""
from pathlib import Path
from .physics import gpe_physics

from jinja2 import Environment, FileSystemLoader

dir = Path(__file__).parent
file_loader = FileSystemLoader(str(dir / 'templates'))
env = Environment(loader=file_loader, trim_blocks=True, lstrip_blocks=True)


def write_simulation(natoms, trap_freqs, ndims,
                     dx=0.33, box_size=2,
                     name='groundstate', dry_run=False,):
    """Writes a xmds GPE simulation to compute the ground state of BEC

    Params:
        natoms
        trap_freqs
        ndims
    """
    template = env.get_template('template.xmds')
    conf = {}

    conf['name'] = f"gpe_{name}"

    mu, a_scatt, xi, radii, msg = gpe_physics(natoms, trap_freqs, ndims)
    g = f"4*M_PI*a_scatt"
    if dry_run:
        print("--- xmgpe dry-run ---")
        print(msg)
        return

    print(msg)

    dim_names = ['x', 'y', 'z'][:ndims]
    trap_freqs = trap_freqs[:ndims]
    conf['dimensions'] = " ".join(dim_names)

    # write globals
    conf['globals'] = dict(N=natoms, a_scatt=a_scatt, g=g)

    # write geometry
    geometry = []
    dx = xi * dx
    for d, r in zip(dim_names, radii):
        L = int(round(r * box_size))
        domain = (-L, L)
        lattice = int(round(2 * L / dx))
        dim = dict(name=d, domain=domain, lattice=lattice)
        geometry.append(dim)

    conf['geometry'] = geometry

    # wavefunction
    conf['gaussian_wavefunction'] = gaussian_wavefunction(dim_names, radii, amp=1)

    # potential
    conf['harmonic_potential'] = harmonic_potential(dim_names, trap_freqs)

    # set integration time
    interval = 1
    dt = 0.5 * dx**2
    steps = interval / dt
    steps = int(round(steps / 100) * 100)
    conf['integrate_interval'] = interval
    conf['integrate_steps'] = steps

    sign = ''
    tx_list = [f"T{dim}[psi]" for dim in dim_names]
    eqn = f"dpsi_dt = {' + '.join(tx_list)} - {sign}(V1 + g*mod2(psi))*psi;"
    conf['gpe'] = eqn
    print(eqn)

    # write to file
    scriptname = f"{name}.xmds"
    output = template.render(conf=conf)
    with open(scriptname, 'w') as f:
        f.write(output)

    print('Done')


def harmonic_potential(names, trap_freqs):
    f_ho = trap_freqs[0]
    v = []
    for name, freq in zip(names, trap_freqs):
        s = "0.5"
        if freq != f_ho:
            q = freq / f_ho
            s += f"*{q:.2f}*{q:.2f}"
        s += f"*{name}*{name}"
        v.append(s)
    v = " + ".join(v)
    return v


def gaussian_wavefunction(dim_names, radii, amp=1):
    g = []
    for name, r in zip(dim_names, radii):
        sigma = r / 3
        s = f"-{name}*{name} / {2.0*sigma**2:.2f}"
        g.append(s)
    g = " ".join(g)
    g = f"{amp:.2f} * exp({g})"
    return g
