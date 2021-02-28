import numpy as np
import matplotlib.pyplot as plt
from palettable.cmocean.sequential import Matter_7
import palettable

def Pe_time_avg_1d(a,d):
	prefactor=0.5*((d)**(-1))
	return prefactor*2*np.arctan(a*d/2)/a

def get_lattice_detuning(r,detuning_type="gaussian",v_trap=106800,norm=True):

	if detuning_type=="gaussian":

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

		ac_stark_2d=2*ac_stark
		if norm==True:
			ac_stark_2d=ac_stark_2d-np.max(ac_stark_2d)

	else:

		half_pix=5.5e-6/2
		slope=(2*np.pi*25)/half_pix
		ac_stark_2d=r*slope

	return slope,ac_stark_2d

def excited_state_fraction(t_array,Omega,Delta):
	excited_fraction=np.zeros([int(np.size(t_array))])

	# print(t_array)

	for n,t in enumerate(t_array):
		A=Omega**2/(Omega**2+Delta**2)
		argument=(Omega**2+Delta**2)**0.5
		excited_fraction[n]=A*(np.sin(argument*t/2)**2)

	return excited_fraction

if __name__ == '__main__':

	# print(s)
	N=1000

	SMALL_SIZE = 12
	MEDIUM_SIZE = 14
	BIGGER_SIZE = 16

	plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
	plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
	plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
	plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
	plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
	plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

	numb_pix_array=[1,4,10,20,50,100]
	fig,ax=plt.subplots(figsize=[10,5])

	for j,numb_pix in enumerate(numb_pix_array):

		tpi=1.5e-3
		Omega=np.pi/tpi
		t=np.linspace(0,20e-3,400)

		a=5.45e-6
		r=np.linspace(-a*numb_pix/2,a*numb_pix/2,N)
		slope,Detuning_array=get_lattice_detuning(r,detuning_type="linear")
		time_avg=Pe_time_avg_1d(slope/Omega,np.abs(r[0])*2)
		print(time_avg)
		# detuning_array_lin=get_lattice_detuning(r,detuning_type="linear")
		# plt.plot(r,detuning_array,label="potential")
		# plt.plot(r,-1*np.abs(detuning_array_lin),label="lin")
		# plt.show()

		# detuning_max=200 #in Hz
		# detuning=np.linspace(0,detuning_max*2*np.pi,N)
		
		for i, Delta in enumerate(Detuning_array):

			if i ==0:
				dynamics=excited_state_fraction(t,Omega,Delta)
			else:
				temp=excited_state_fraction(t,Omega,Delta)
				dynamics+=temp

		# on_res=excited_state_fraction(t,Omega,0)
		# detuned=excited_state_fraction(t,Omega,2*np.pi*400)
		dynamics=dynamics/N
		ax.plot(t*1e3,dynamics,":o",label=r"numb pix: "+str(numb_pix),color=palettable.colorbrewer.qualitative.Set1_6.hex_colors[j],"linewdith")
		ax.plot(t*1e3,time_avg*np.ones(len(t)),color=palettable.colorbrewer.qualitative.Set1_6.hex_colors[j])

	ax.set_xlabel('time (ms)',fontsize=15)
	ax.set_ylabel('excitation fraction',fontsize=15)
	ax.legend(loc="upper left")

	fig.savefig('excfraction_VS_time_lin_1D.pdf')
	# plt.plot(t,detuned,":o",label="detuned")	
	plt.show()

		