#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Create: 07-2019 - Carmelo Mordini <carmelo> <carmelo.mordini@unitn.it>

"""Module docstring

"""
from lxml import etree
from math import pi, sqrt
from .physics import transform_prefactors

s2pi = sqrt(2 * pi)


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


def read_xmds(filename):
    parser = etree.XMLParser(strip_cdata=False, remove_blank_text=True)
    with open(filename, "r") as source:
        main = etree.parse(source, parser=parser)
    return main


def write_xmds(tree, scriptname):
    print(f"Write {scriptname}")
    indent(tree.getroot(), level=0)
    # writer.pretty_print(main)
    tree.write(scriptname)


def globals_cdata(**kwargs):
    globs = {}
    for k, v in kwargs.items():
        globs[k] = v if isinstance(v, str) else f"{v:.6f}"
    sep = '\n\t'
    globs_str = sep.join([f"const real {key} = {value};" for key, value in globs.items()])
    return etree.CDATA(sep + globs_str + sep)


def dimension(name, domain, lattice, transform=None,):
    dim = etree.Element('dimension')
    dim.attrib['name'] = name
    print(f"Building dimension {name}")
    if transform is not None and transform != "dft":
        dim.attrib['transform'] = transform
        dim.attrib['volume_prefactor'] = transform_prefactors[transform]
    left, right = domain
    dim.attrib['domain'] = f"({left:.3f}, {right:.3f})"
    dim.attrib['lattice'] = str(lattice)
    return dim


def geometry(dimensions):
    # geom = et.find('geometry')
    geom = etree.Element('geometry')
    td = etree.SubElement(geom, 'propagation_dimension')
    td.text = 't'
    td = etree.SubElement(geom, 'transverse_dimensions')
    for dim in dimensions:
        td.append(dim)
    return geom


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


def wavefunction(geometry_conf, initialisation='thomasfermi'):
    dimensions = " ".join([dim['name'] for dim in geometry_conf])
    wf = etree.Element("vector", name='wavefunction',
                       dimensions=dimensions, type='complex')
    comp = etree.SubElement(wf, 'components')
    comp.text = 'psi'
    init = etree.SubElement(wf, 'initialisation')
    if initialisation == 'thomasfermi':
        dep = etree.SubElement(init, 'dependencies')
        dep.text = 'potential'
        init.text = etree.CDATA(
            "psi = (mu0 - V1)/Uint;\n\tif (psi.Re() >= 0) psi = sqrt(psi); else psi = 0.0;")
    elif initialisation.endswith('.h5'):
        init.attrib['kind'] = 'hdf5'
        filen = etree.SubElement(init, 'filename')
        filen.text = initialisation
    else:
        raise ValueError('Wrong initialisation argument')
    return wf


def filter(dependencies, cdata):
    filt = etree.Element('filter')
    dep = etree.SubElement(filt, 'dependencies')
    dep.text = dependencies
    dep.tail = etree.CDATA(cdata)
    return filt


def ip_operator(dim, evolution='real'):
    op = etree.Element('operator', kind='ip')
    opn = etree.SubElement(op, 'operator_names')
    opn.text = f"T{dim}"
    if evolution == 'real':
        sign = '-i*'
    elif evolution == 'imaginary':
        sign = '-'
    else:
        raise TypeError("Wrong time evolution: real or imaginary")
    opn.tail = etree.CDATA(f"T{dim} = {sign}0.5*k{dim}*k{dim};")
    return op
