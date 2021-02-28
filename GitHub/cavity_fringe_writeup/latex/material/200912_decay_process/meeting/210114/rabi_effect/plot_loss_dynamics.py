import numpy as np
import srlab.plot.pyplot as plt

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
	N=5000

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

	fig,ax=plt.subplots(figsize=[10,5])

	tpi=1.5e-3
	Omega=np.pi/tpi
	max_omega_array=[Omega*0.01,Omega*0.1,Omega*0.3]


	for j, max_omega in enumerate(max_omega_array):
		if j==0 or j==4:
			Omega_array=sample_gaussian(mu=Omega,sigma=max_omega,N=N,plot=True)
		else:
			Omega_array=sample_gaussian(mu=Omega,sigma=max_omega,N=N)


		t=np.linspace(0,20e-3,400)

		# detuning_max=200 #in Hz
		# detuning=np.linspace(0,detuning_max*2*np.pi,N)
		
		for i, Omegai in enumerate(Omega_array):

			if i ==0:
				dynamics=excited_state_fraction(t,Omegai,0)
			else:
				temp=excited_state_fraction(t,Omegai,0)
				dynamics+=temp

		# on_res=excited_state_fraction(t,Omega,0)
		# detuned=excited_state_fraction(t,Omega,2*np.pi*400)
		dynamics=dynamics/N
		ax.plot(t*1e3,dynamics,":o",label=r"$\Omega_{\sigma}$= 2*pi*"+str(max_omega/2/np.pi) + " Hz")
	ax.set_xlabel('time (ms)',fontsize=15)
	ax.set_ylabel('excitation fraction',fontsize=15)
	ax.legend()

	fig.savefig('excfraction_VS_time.pdf')
	# plt.plot(t,detuned,":o",label="detuned")	
	plt.show()


