#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Create: 08-2019 - Carmelo Mordini <carmelo> <carmelo.mordini@unitn.it>

"""Module docstring

"""
from lxml import etree
from math import sqrt
from . import writer
from .thomasfermi import pi, m_u, hbar, abohr
from ruamel.yaml import YAML
yaml = YAML(typ='safe')


solid_angle_d = {1: 2, 2: 2 * pi, 3: 4 * pi}


def chem_pot_d(N, Uint, ndim):
    return 0.5*(ndim*(ndim+2)/2/solid_angle_d[ndim] * Uint * N)**(2./(ndim+2))


def write_simulation(filename=None,):
    with open('xmgpe/configure.yaml') as f:
        conf = yaml.load(f)
    parser = etree.XMLParser(strip_cdata=False, remove_blank_text=True)
    with open('xmgpe/template.xmds', "r") as source:
        main = etree.parse(source, parser=parser)

    # read physical parameters
    Natoms = conf['Natoms']
    mass = conf['mass'] * m_u
    scatt_len = conf['scatt_len'] * abohr
    geometry_conf = conf['geometry']
    ndim = len(geometry_conf)
    dim_names = [dim['name'] for dim in geometry_conf]
    omega_ho = 2 * pi * geometry_conf[-1]['trap_freq']
    print(f"in main: omega ho = {omega_ho / 2 / pi}")

    a_ho = sqrt(hbar / mass / omega_ho)
    Uint = 4 * pi * scatt_len / a_ho
    if ndim < 3:
        omega_tr = 2 * pi * conf['transverse_trap_freq']
        a_tr = sqrt(hbar / mass / omega_tr)
        Uint *= (a_ho / sqrt(2 * pi) / a_tr)**(3 - ndim)
    mu0 = chem_pot_d(Natoms, Uint, ndim)

    # write globals
    features = main.find('features')
    globals = features.find('globals')
    globals.text = writer.globals_cdata(Natoms=Natoms, Uint=Uint, mu0=mu0)

    # write geometry
    geom = writer.geometry(geometry_conf, mu0=mu0)
    main.getroot().replace(main.find("geometry"), geom)

    # write vectors
    vectors = []
    vpot = writer.potential(geometry_conf)
    writer.pretty_print(vpot)

    wavef = writer.wavefunction(geometry_conf, initialisation='thomasfermi')
    writer.pretty_print(wavef)

    # write to file
    writer.indent(main.getroot(), level=0)
    # writer.pretty_print(main)
    main.write('newscript.xmds')
