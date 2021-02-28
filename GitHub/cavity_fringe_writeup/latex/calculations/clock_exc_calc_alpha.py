import srlab.units as u
import numpy as np

if __name__ == '__main__':

	q=1.60217662e-19*u.Real("C")
	D=0.151*u.bohrradius*q#*np.sqrt(3)

	mu=np.sqrt(2/3)*u.bohrmagneton

	Delta_32=5.6*u.Real("THz")*2*np.pi

	# print(Test)
	prefactor=mu*D/((u.hbar**2)*Delta_32)#*((2/u.c.value/u.vacuumpermittivity.value)**0.5)
	postfactor=(1/2/u.c/u.vacuumpermittivity)#**0.5

	alpha=mu.value*D.value/((u.hbar.value**2)*Delta_32.value)/((2*u.c.value*u.vacuumpermittivity.value)**0.5)

	print(alpha/2/np.pi*np.sqrt(10),prefactor.value*(postfactor.value**0.5)/2/np.pi*np.sqrt(10))
