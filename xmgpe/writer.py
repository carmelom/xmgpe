#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Create: 07-2019 - Carmelo Mordini <carmelo> <carmelo.mordini@unitn.it>

"""Module docstring

"""
from math import pi, sqrt
from lxml import etree
from .validation import transform_prefactors


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


def pretty_print(elem):
    print(etree.tostring(elem, encoding='unicode', pretty_print=True))


def globals_cdata(**kwargs):
    const = [f"{key} = {value:.6f}" for key, value in kwargs.items()]
    sep = '\n\t'
    cdata_str = sep.join([f'const real {s};' for s in const])
    return etree.CDATA(sep + cdata_str + sep)


def dimension(name, lattice, domain, transform=None, mu0=1, f_ho=1, **kwargs):
    dim = etree.Element('dimension')
    dim.attrib['name'] = name
    print(f"Building dimension {name}")
    if transform is not None and transform != "dft":
        dim.attrib['transform'] = transform
        dim.attrib['volume_prefactor'] = transform_prefactors[transform]
    if isinstance(domain, list):
        left, right = domain
    else:
        if isinstance(domain, (int, float)):
            right = domain
        elif domain.startswith("x"):
            mult = float(domain[1:])
            print(f"domain len: {mult} Rtf")
            Rtf = sqrt(2 * mu0) * f_ho / kwargs['trap_freq']
            right = mult * Rtf
        if transform is not None and transform != "dft":
            left = 0
        else:
            left = -right
    dim.attrib['domain'] = f"({left:.3f}, {right:.3f})"
    if isinstance(lattice, str) and lattice.startswith("x"):
        mult = float(lattice[1:])
        print(f"lattice spacing: {mult} healing length")
        healing_length = 1 / sqrt(2 * mu0)
        dx = mult * healing_length
        lattice = int((right - left) / dx)
    dim.attrib['lattice'] = str(lattice)
    return dim


def geometry(geometry_conf, mu0=1):
    # geom = et.find('geometry')
    geom = etree.Element('geometry')
    td = etree.SubElement(geom, 'propagation_dimension')
    td.text = 't'
    td = etree.SubElement(geom, 'transverse_dimensions')
    f_ho = geometry_conf[-1]['trap_freq']
    for dimension_conf in geometry_conf:
        dim = dimension(mu0=mu0, f_ho=f_ho, **dimension_conf)
        td.append(dim)
    return geom


def potential(geometry_conf):
    f_ho = geometry_conf[-1]['trap_freq']
    dimensions = " ".join([dim['name'] for dim in geometry_conf])
    cdata = "V1 = "
    v = []
    for dim in geometry_conf:
        name = dim['name']
        s = f"0.5*{name}*{name}"
        if dim['trap_freq'] != f_ho:
            q = dim['trap_freq'] / f_ho
            s += f"*{q:.3f}*{q:.3f}"
        v.append(s)
    cdata = cdata + " + ".join(v) + ";"
    vpot = etree.Element("vector", name='potential', dimensions=dimensions, type='real')
    comp = etree.SubElement(vpot, 'components')
    comp.text = 'V1'
    init = etree.SubElement(vpot, 'initialisation')
    init.text = etree.CDATA(cdata)
    return vpot


def wavefunction(geometry_conf, initialisation='thomasfermi'):
    dimensions = " ".join([dim['name'] for dim in geometry_conf])
    wf = etree.Element("vector", name='wavefunction', dimensions=dimensions, type='complex')
    comp = etree.SubElement(wf, 'components')
    comp.text = 'psi'
    init = etree.SubElement(wf, 'initialisation')
    if initialisation == 'thomasfermi':
        dep = etree.SubElement(init, 'dependencies')
        dep.text = 'potential'
        init.text = etree.CDATA("psi = (mu0 - V1)/Uint;\n\tif (psi >= 0) psi = sqrt(psi); else psi = 0.0;")
    elif initialisation.endswith('.h5'):
        init.attrib['kind'] = 'hdf5'
        filen = etree.SubElement(init, 'filename')
        filen.text = initialisation
    else:
        raise ValueError('Wrong initialisation argument')
    return wf
