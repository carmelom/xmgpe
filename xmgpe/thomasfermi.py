#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Create: 07-2019 - Carmelo Mordini <carmelo> <carmelo.mordini@unitn.it>

"""Module docstring

"""

# import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import pi, hbar, m_u, physical_constants
from ruamel.yaml import YAML
yaml = YAML(typ='safe')

abohr = physical_constants['Bohr radius'][0]


def Uint():
    return 4*np.pi * scatt_len / a_ho


def mu_tf(N, Uint):
    return 0.5*(15*N*Uint / 4/np.pi)**2/5
