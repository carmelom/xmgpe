#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Create: 10-2019 - Carmelo Mordini <carmelo> <carmelo.mordini@unitn.it>

"""Module docstring

"""
import numpy as np
from math import sqrt
from scipy.constants import pi, hbar, atomic_mass, Boltzmann as kB, physical_constants

a_bohr = physical_constants['Bohr radius'][0]

solid_angle_d = {1: 2, 2: 2 * pi, 3: 4 * pi}

transform_prefactors = {'dct': '2.0',
                        'dst': '2.0',
                        'bessel': '2*M_PI',
                        'bessel-neumann': '2*M_PI',
                        'spherical-bessel': '4*M_PI', }


def chem_pot_d(N, g, ndims=3):
    return 0.5 * (ndims * (ndims + 2) / solid_angle_d[ndims] * g * N)**(2. / (ndims + 2))


def gpe_physics(natoms, trap_freqs, ndims, mass=23, a_scatt=52):
    assert len(trap_freqs) == 3
    assert ndims <= 3
    trap_freqs = np.asarray(trap_freqs)
    mass = mass * atomic_mass
    a_scatt = a_scatt * a_bohr
    omega = 2 * pi * trap_freqs
    omega_x = omega[0]
    omega_ho = np.prod(omega)**(1. / 3)
    a_x = sqrt(hbar / mass / omega_x)
    a_ho = sqrt(hbar / mass / omega_ho)
    mu = chem_pot_d(natoms, g=4 * pi * a_scatt / a_ho, ndims=3) * \
        omega_ho / omega_x  # norm. wrt omega_x
    a_scatt_1 = a_scatt / a_x
    g = 4 * pi * a_scatt_1
    xi = 1 / sqrt(2 * mu)  # norm. wrt a_x
    # compute TF radii along the simulated  dimensions
    _omega = omega[:ndims]
    Rtf = sqrt(2 * mu) * omega_x / _omega
    np.set_printoptions(precision=2)
    msg = f"""Setting up xmds for an imaginary time GPE with:
N: {natoms:.2e} atoms
trap freqs:
    {trap_freqs}
simulated dimensions: {ndims}
normalized wrt a_x = {a_x*1e6:.2f} um
TF radii:
    {Rtf*a_x*1e6} nm ({Rtf})
scattering length:
    {a_scatt_1:.6f} ({a_scatt*1e9:.2f} nm)
interaction constant:
    {g:.6f}
chemical potential:
    {mu*hbar*omega_x*1e9/kB:.2f} nK ({mu:.2f} hbar*omega_x)
healing length:
    {xi*a_x*1e6:.2f} nm ({xi:.2f})
"""
    return mu, a_scatt_1, xi, Rtf, msg
