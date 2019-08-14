#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Create: 08-2019 - Carmelo Mordini <carmelo> <carmelo.mordini@unitn.it>

"""Module docstring

"""
from pprint import pprint
from xmgpe import writer
from xmgpe.main import write_simulation

from ruamel.yaml import YAML
yaml = YAML(typ='safe')


def test_geometry():
    with open('xmgpe/configure.yaml') as f:
        params = yaml.load(f)
        pprint(params)
        geom = writer.geometry(params['geometry'])
        writer.pretty_print(geom)


if __name__ == '__main__':
    # test_geometry()
    write_simulation()
