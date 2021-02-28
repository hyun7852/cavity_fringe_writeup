import numpy as np
import srlab.plot.pyplot as plt
import scipy.constants as const

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

def sample_gaussian(mu=0,sigma=200,N=1000,plot=False):
	# mu, sigma = 0, 200 # mean and standard deviation
	s = np.random.normal(mu, sigma, N)
	if plot==True:
		fig2,ax2=plt.subplots(figsize=[7,5])
		count, bins, ignored = plt.hist(s, 30, density=True)
		ax2.plot(bins, 1/(sigma*np.sqrt(2*np.pi))*np.exp(-(bins-mu)**2 / (2*sigma**2)), linewidth=2, color='r')	
		ax2.set_xlabel(r'Detuning ($\Delta$)')
		ax2.set_ylabel('Counts')
		fig2.savefig('Sigma_'+str(sigma)+'.pdf')
	# plt.show()
	return s
    # return 

def get_lattice_detuning(r,v_trap=106800,norm=True):

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

    return ac_stark_2d

def excited_state_fraction(t_array,Omega,Delta):
	excited_fraction=np.zeros([int(np.size(t_array))])

	# print(t_array)

	for n,t in enumerate(t_array):
		A=Omega**2/(Omega**2+Delta**2)
		argument=(Omega**2+Delta**2)**0.5
		excited_fraction[n]=A*(np.sin(argument*t/2)**2)

	return excited_fraction

# check whether 1D makes sense and extend to 2D

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

	numb_pix_array=[1,4,10,20,50]
	fig,ax=plt.subplots(figsize=[10,5])

	for j,numb_pix in enumerate(numb_pix_array):

		r=np.linspace(-5.45e-6*numb_pix/2,5.45e-6*numb_pix/2,N)
		detuning_array=get_lattice_detuning(r)


		tpi=1.5e-3
		Omega=np.pi/tpi
		t=np.linspace(0,20e-3,400)

		# detuning_max=200 #in Hz
		# detuning=np.linspace(0,detuning_max*2*np.pi,N)
		
		for i, delta in enumerate(detuning_array):

			if i ==0:
				dynamics=excited_state_fraction(t,Omega,delta)
			else:
				temp=excited_state_fraction(t,Omega,delta)
				dynamics+=temp

		# on_res=excited_state_fraction(t,Omega,0)
		# detuned=excited_state_fraction(t,Omega,2*np.pi*400)
		dynamics=dynamics/N
		ax.plot(t*1e3,dynamics,":o",label=r"numb pix: "+str(numb_pix))

	ax.set_xlabel('time (ms)',fontsize=15)
	ax.set_ylabel('excitation fraction',fontsize=15)
	ax.legend()

	fig.savefig('excfraction_VS_time.pdf')
	# plt.plot(t,detuned,":o",label="detuned")	
	plt.show()


