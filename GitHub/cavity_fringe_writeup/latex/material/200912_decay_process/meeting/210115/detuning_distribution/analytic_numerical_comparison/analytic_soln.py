import numpy as np
import matplotlib.pyplot as plt

def Pe_time_avg_2d(a,b,d):
	arg=2*a*b*(d**2)/(4+(a**2+b**2)*(d**2))
	print("arg is:",np.arctanh(arg))
	t=np.arctanh(arg)
	prefactor=0.5*((d)**(-2))
	return prefactor/a/b*((b-a)*d*np.arctan((a-b)*d/2)+(a+b)*d*np.arctan((a+b)*d/2)-2*t)

def Pe_time_avg_1d(a,d):
	prefactor=0.5*((d)**(-1))
	return prefactor*2*np.arctan(a*d/2)/a

if __name__ == '__main__':
	
	tpi=1.5e-3
	Omega=np.pi/tpi
	d=5.5e-6/2*10
	J=25/(5.5e-6/2)

	a=J/Omega
	b=J/Omega

	print("Pe",Pe_time_avg_1d(a,d))

