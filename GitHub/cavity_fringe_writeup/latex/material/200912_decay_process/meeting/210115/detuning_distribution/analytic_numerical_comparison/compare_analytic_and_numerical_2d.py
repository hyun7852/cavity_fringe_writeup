import numpy as np
import matplotlib.pyplot as plt
from palettable.cmocean.sequential import Matter_7
import palettable

def Pe_time_avg_2d(a,b,d):
    arg=2*a*b*(d**2)/(4+(a**2+b**2)*(d**2))
    print("arg is:",np.arctanh(arg))
    t=np.arctanh(arg)
    prefactor=0.5*((d)**(-2))
    return prefactor/a/b*((b-a)*d*np.arctan((a-b)*d/2)+(a+b)*d*np.arctan((a+b)*d/2)-2*t)

def get_lattice_detuning(xx,yy,v_trap=106800,norm=True):

        half_pix=5.5e-6/2
        slope=(2*np.pi*25)/half_pix
        ac_stark_2d=slope*(xx+yy)
        # ac_stark_2d=r*slope

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
    N=100

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
        x=np.linspace(-a*numb_pix/2,a*numb_pix/2,N)
        xx, yy = np.meshgrid(x,x)
        

        slope,Detuning_array=get_lattice_detuning(xx,yy)
        time_avg=Pe_time_avg_2d(slope/Omega,slope/Omega,np.abs(x[0])*2)

        print(time_avg)

        # detuning_array_lin=get_lattice_detuning(r,detuning_type="linear")
        # plt.plot(r,detuning_array,label="potential")
        # plt.plot(r,-1*np.abs(detuning_array_lin),label="lin")
        # plt.show()

        # detuning_max=200 #in Hz
        # detuning=np.linspace(0,detuning_max*2*np.pi,N)
        
        dynamics=0
        for i, Delta_i in enumerate(Detuning_array):
            for k, Delta_ij in enumerate(Delta_i):
                temp=excited_state_fraction(t,Omega,Delta_ij)
                dynamics+=temp

        # on_res=excited_state_fraction(t,Omega,0)
        # detuned=excited_state_fraction(t,Omega,2*np.pi*400)
        dynamics=dynamics/N/N
        ax.plot(t*1e3,dynamics,":o",label=r"numb pix: "+str(numb_pix),color=palettable.colorbrewer.qualitative.Set1_6.hex_colors[j])
        ax.plot(t*1e3,time_avg*np.ones(len(t)),color=palettable.colorbrewer.qualitative.Set1_6.hex_colors[j])

    ax.set_xlabel('time (ms)',fontsize=15)
    ax.set_ylabel('excitation fraction',fontsize=15)
    ax.legend(loc="upper left")

    fig.savefig('excfraction_VS_time_lin_2D.pdf')
    # plt.plot(t,detuned,":o",label="detuned")  
    plt.show()

        