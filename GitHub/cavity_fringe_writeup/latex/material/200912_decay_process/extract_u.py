#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2020-12-09 11:15:56 ajp"

#  file       extract_u.py
#  author     Annie Park
#  created    2020-12-09 11:15:56

import srlab.units as u
import numpy as np
import scipy.constants as const
import srlab.plot.pyplot as plt

# constants
c = const.c
e_0 = const.epsilon_0
h = const.h
hbar = const.hbar
kB = const.k
a_0 = const.physical_constants['atomic unit of length'][0]
amu = const.physical_constants['atomic mass constant'][0]
au = 4 * np.pi* e_0 * a_0**3
m = 88 * amu

########################################
# Function to calculate polarizability #
# Reference: Rauschenbeutel paper      #
########################################

# calculate recoil energy
def cal_E_rec(wavelength):
    wl = wavelength *1e-9
    E_rec = h**2 /(2 * m * wl**2)

    return E_rec

if __name__ == '__main__':
	
	ag=253.3
	ae=214.7
	Erec=cal_E_rec(914.322)
	vt=(106500+1*Erec/h)*np.sqrt(ag/ae)
	U=((vt*h)**2/(4*Erec))*2
	U=U*(1-ae/ag)
	print("potential in kHz",U/h/1000)

