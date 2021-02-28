import numpy as np
import srlab.plot.pyplot as plt
import random

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

def excited_state_fraction(t_array,Omega,Delta):
	excited_fraction=np.zeros([int(np.size(t_array))])

	# print(t_array)

	if np.size(t_array)==1:
		A=Omega**2/(Omega**2+Delta**2)
		argument=(Omega**2+Delta**2)**0.5
		excited_fraction=A*(np.sin(argument*t_array/2)**2)

	else:
		for n,t in enumerate(t_array):
			A=Omega**2/(Omega**2+Delta**2)
			argument=(Omega**2+Delta**2)**0.5
			excited_fraction[n]=A*(np.sin(argument*t/2)**2)

	return excited_fraction

if __name__ == '__main__':
	
	# print(s)
	N=100
	# Nt=len(np.array([0.1,0.3,0.5,0.75,1,1.15,2,2.5,3.5,5,6,7,10,15])*1e-3)
	Nt=100

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

	max_detuning_array=[400]
	fig,ax=plt.subplots(figsize=[10,5])

	avg_dynamics=np.zeros([N,Nt])
	# print(np.size(avg_dynamics))
	for k in range(N):
		dynamics=[]
		for j, delta_max in enumerate(max_detuning_array):


			tpi=1.7e-3
			Omega=np.pi/tpi
			# t=np.array([0.1,0.3,0.5,0.75,1,1.15,2,2.5,3.5,5,6,7,10,15])*1e-3
			t=np.linspace(0,20e-3,Nt)
			# random.shuffle(t)

			# detuning_max=200 #in Hz
			# detuning=np.linspace(0,detuning_max*2*np.pi,N)
			
			
			for i, ti in enumerate(t):

				Delta=sample_gaussian(mu=0,sigma=delta_max*2*np.pi,N=1)
				dynamics.append(excited_state_fraction(ti,Omega,Delta)[0])

		# on_res=excited_state_fraction(t,Omega,0)
		# detuned=excited_state_fraction(t,Omega,2*np.pi*400)
	
		SORT=True
		# print(dynamics)

		dynamics=np.array(dynamics)
		if SORT:
			sort_order=np.argsort(t)
			t_sorted=t[sort_order]
			dynamics_sorted=dynamics[sort_order]
		else:
			t_sorted=t
			dynamics_sorted=dynamics
		# print("dynamcis size",np.shape(dynamics))

		# print("dynamcis sorted size",np.shape(dynamics_sorted))
		# print("avg dynamics size",np.shape(avg_dynamics))

		avg_dynamics[k,:]=dynamics_sorted
	# print(np.shape(avg_dynamics[0]))
		# avg_dynamics=np.array(avg_dynamics) + np.array(dynamics_sorted)

		# print(avg_dynamics)
	# avg_dynamics=avg_dynamics/N

	# print(np.size(avg_dynamics))

		
	# ax.plot(t_sorted*1e3,avg_dynamics[k],":o",label=r"$\delta_{\sigma}$= "+str(delta_max) + " Hz")
	y=np.mean(avg_dynamics,0)
	ax.plot(t_sorted*1e3,y,":o",color='b',alpha=0.6,label=r"$\delta_{\sigma}$= "+str(delta_max) + " Hz")
	ci=1.96*np.std(avg_dynamics,0)/np.sqrt(N)
	ax.fill_between(t_sorted*1e3,y-ci,y+ci,color='b',alpha=0.1)
	ax.set_xlabel('time (ms)',fontsize=15)
	ax.set_ylabel('excitation fraction',fontsize=15)
	ax.legend()

	fig.savefig('excfraction_VS_time.pdf')
	# plt.plot(t,detuned,":o",label="detuned")	
	plt.show()


