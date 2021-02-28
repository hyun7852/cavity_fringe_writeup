#!/usr/bin/env python
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2020-12-09 11:15:56 ajp"

#  file       plot_acstark.py
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
m = 88 * amu

CAVITY_WAIST=455.3e-6*1.0

# calculate recoil energy
def cal_E_rec(wavelength):
    wl = wavelength *1e-9
    E_rec = h**2 /(2 * m * wl**2)

    return E_rec

def intensity_profile(r,I0):
	
	# # at z=0
	# I0=2*P0/np.pi/CAVITY_WAIST

	return I0*np.exp(-2*((r/CAVITY_WAIST)**2))

def plot_gaussian_envelope(r,v_trap,norm=False,figname='ac_stark_914nm.pdf'):
	
	#Sr property at 914 lattice
    SS0_pol = 253
    TP0_pol = 215
    
    au = 4 * np.pi* e_0 * a_0**3
    dp=223-253
    dp=dp*au
 
    Erec=cal_E_rec(914)
    I0=((v_trap*h/4)**2)*(2*e_0*c/Erec/(SS0_pol*au))

    potential=-1*dp*intensity_profile(r,I0)/(2*e_0*c)
    ac_stark=potential/hbar

    # multiply by 2 for 2D lattice
    plt.figure(figsize=(12*0.7,9*0.7))
    if norm==True:
        plt.plot(r*1e3,(ac_stark-np.max(ac_stark))*2)
        plt.ylabel('AC stark shift (Hz)')
    else:
        plt.plot(r*1e3,(ac_stark)*2/1e3)
        plt.ylabel('AC stark shift (kHz)')
    plt.xlabel('x (mm)')
    # plt.ylim([385.1,385.6])
    plt.grid()
    plt.savefig(figname)
    # plt.show()

if __name__ == '__main__':

	# r=np.linspace(-1e-3,1e-3,500)
	# plot_gaussian_envelope(r,106000)

	#in the center for four pixels 5.45*4
	numb_pix=1
	r=np.linspace(-0,5.45e-6*numb_pix/2,200)
	plot_gaussian_envelope(r,106800,norm=True,figname='ac_stark_914nm_center_assym'+str(numb_pix)+'pix_detuning.pdf')

