import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

# constants
c = const.c
e_0 = const.epsilon_0
h = const.h
hbar = const.hbar
kB = const.k
a_0 = const.physical_constants['atomic unit of length'][0]
amu = const.physical_constants['atomic mass constant'][0]

def boltzmann(vtrap,T):
	return np.exp(-h*vtrap/kB/T)


if __name__ == '__main__':
	
	angle=50e-3
	nx=angle*3.2
	nz=0.24

	vz=80000 # in Hz
	vx=450
	zx=boltzmann(vx,3e-6)

	t=np.linspace(0,0.1,300)
	Omega=55 # in Hz
	phi=2*np.pi*t*Omega*np.exp(-nx**2/2)*np.exp(-nz**2/2)
	exc_fraction=0.5+((1-zx)/2)*(zx*np.cos(phi*(1-nx**2))-np.cos(phi))/(1+zx**2-2*zx*np.cos(phi*nx*nx))

	plt.plot(t,exc_fraction)
	plt.show()




