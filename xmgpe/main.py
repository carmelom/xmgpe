#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Create: 08-2019 - Carmelo Mordini <carmelo> <carmelo.mordini@unitn.it>

"""Module docstring

"""
from pathlib import Path
from lxml import etree
from . import writer, physics

from ruamel.yaml import YAML
yaml = YAML(typ='safe')

dir = Path(__file__).parent

# with open(conf_filename) as f:
#     conf = yaml.load(f)


def write_simulation(natoms, trap_freqs, ndims,
                     dx=0.33, box_size=2,
                     name='groundstate', dry_run=False,):
    """Writes a xmds GPE simulation to compute the ground state of BEC

    Params:
        natoms
        trap_freqs
        ndims
    """
    template = dir / 'templates/template.xmds'
    main = writer.read_xmds(template)
    # for tag in ['name', 'author', 'description']:
    #     main.find(tag).text = conf[tag]
    main.find('name').text = f"gpe_{name}"

    mu, a_scatt, xi, radii, msg = physics.gpe_physics(natoms, trap_freqs, ndims)
    g = f"4*M_PI*a_scatt"
    if dry_run:
        print("--- xmgpe dry-run ---")
        print(msg)
        return

    print(msg)

    dim_names = ['x', 'y', 'z'][:ndims]
    trap_freqs = trap_freqs[:ndims]

    # write globals
    features = main.find('features')
    globals = features.find('globals')
    globs = dict(N=natoms, a_scatt=a_scatt, g=g)
    globals.text = writer.globals_cdata(**globs)

    # write geometry
    dimensions = []
    dx = xi * dx
    for d, r in zip(dim_names, radii):
        L = int(round(r * box_size))
        domain = (-L, L)
        lattice = int(round(2 * L / dx))
        dim = writer.dimension(d, domain, lattice)
        dimensions.append(dim)

    geom = writer.geometry(dimensions)
    main.getroot().replace(main.find("geometry"), geom)

    # wavefunction
    wavef = main.xpath("//vector[@name='wavefunction']")[0]
    wavef.attrib['dimensions'] = " ".join(dim_names)
    psi = writer.gaussian_wavefunction(dim_names, radii, amp=1)
    wavef.find('initialisation').text = etree.CDATA(f"psi = {psi};")

    # potential
    vpot = main.xpath("//vector[@name='potential']")[0]
    vpot.attrib['dimensions'] = " ".join(dim_names)
    v1 = writer.harmonic_potential(dim_names, trap_freqs)
    vpot.find('initialisation').text = etree.CDATA(f"V1 = {v1};")

    sequence = main.find('sequence')
    integrate = sequence.find('integrate')

    # set integration time
    interval = 1
    dt = 0.5 * dx**2
    steps = interval / dt
    steps = int(round(steps / 100) * 100)
    integrate.attrib['interval'] = str(interval)
    integrate.attrib['steps'] = str(steps)

    operators = integrate.find('operators')
    for child in operators:
        if child.tag == 'operator':
            operators.remove(child)
    tx_list = []
    for j, dim in enumerate(dim_names):
        op = writer.ip_operator(dim, evolution='imaginary')
        operators.insert(j, op)
        tx_list.append(f"T{dim}[psi]")
    # sign = 'i*' if conf['time_evolution'] == 'real' else ''
    sign = ''
    eqn = f"dpsi_dt = {' + '.join(tx_list)} - {sign}(V1 + g*mod2(psi))*psi;"
    print(eqn)
    integrate.replace(integrate.find("operators"), operators)
    sep = '\n\t\t'
    operators.find('dependencies').tail = etree.CDATA(sep + eqn + sep)

    # explicit output names
    breakpoint = sequence.find('breakpoint')
    breakpoint.attrib['filename'] = f"gpe_{name}_final"
    breakpoint.find('dependencies').attrib['basis'] = " ".join(dim_names)

    main.find('output').attrib['filename'] = f"gpe_{name}_results"

    # fix geometry in first sampling group
    samp0 = main.xpath("//output/sampling_group")[0]
    samp0.attrib['basis'] = " ".join(dim_names)

    # write to file
    scriptname = f"{name}.xmds"
    writer.write_xmds(main, scriptname)
    print('Done')
