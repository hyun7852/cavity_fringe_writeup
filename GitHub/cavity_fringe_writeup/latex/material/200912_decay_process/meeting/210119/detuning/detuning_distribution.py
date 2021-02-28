#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2021-01-20 11:15:56 ajp"

#  file       detuning_distribution.py
#  author     Annie Park
#  created    2021-01-20 11:15:56

import numpy as np
import srlab.plot.pyplot as plt
import random

# constants
c = const.c
e_0 = const.epsilon_0
h = const.h
hbar = const.hbar
kB = const.k
a_0 = const.physical_constants['atomic unit of length'][0]
amu = const.physical_constants['atomic mass constant'][0]
m = 88 * amu

CAVITY_WAIST=455.3e-6*1.0

class Sample(object):
	def __init__(self,waist,v_lattice_trap,v_z_trap,wavelength=914e-9):
		self.waist=waist
		self.v_lattice_trap=v_lattice_trap
		self.v_z_trap=v_z_trap
		self.wavelength=wavelength

	def cal_E_rec(self):
		return h**2 /(2 * m * self.wavelength**2)

	def get_acstark(self):
		# 3P0: 223
		# 1S0: 253
	    au = 4 * np.pi* e_0 * a_0**3
	    dp=223-253
	    dp=dp*au	

	    Erec=self.cal_E_rec()
	    




if __name__ == '__main__':
	main()