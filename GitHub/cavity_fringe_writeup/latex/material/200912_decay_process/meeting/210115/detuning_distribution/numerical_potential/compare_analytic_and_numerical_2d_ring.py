import numpy as np
import matplotlib.pyplot as plt
from palettable.cmocean.sequential import Matter_7
import palettable
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


def Pe_time_avg_2d(a,b,d):
    arg=2*a*b*(d**2)/(4+(a**2+b**2)*(d**2))
    print("arg is:",np.arctanh(arg))
    t=np.arctanh(arg)
    prefactor=0.5*((d)**(-2))
    return prefactor/a/b*((b-a)*d*np.arctan((a-b)*d/2)+(a+b)*d*np.arctan((a+b)*d/2)-2*t)

def intensity_profile(xx,yy,I0):
    
    # # at z=0
    # I0=2*P0/np.pi/CAVITY_WAIST
    r2=(xx**2+yy**2)
    return I0*np.exp(-2*((r2/(CAVITY_WAIST**2))))

def cal_E_rec(wavelength):
    wl = wavelength *1e-9
    E_rec = h**2 /(2 * m * wl**2)

    return E_rec

def get_lattice_detuning(xx,yy,v_trap=106800,norm=True):

    #Sr property at 914 lattice
    SS0_pol = 253
    TP0_pol = 215
    
    au = 4 * np.pi* e_0 * a_0**3
    dp=223-253
    dp=dp*au
 
    Erec=cal_E_rec(914)
    I0=((v_trap*h/4)**2)*(2*e_0*c/Erec/(SS0_pol*au))

    potential=-1*dp*intensity_profile(xx,yy,I0)/(2*e_0*c)
    ac_stark=potential/hbar

    ac_stark_2d=2*ac_stark
    if norm==True:
        print(np.max(ac_stark_2d))
        ac_stark_2d=ac_stark_2d-np.max(ac_stark_2d)

    return ac_stark_2d*2*np.pi

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

    numb_pix_array=[1,2,3,4]
    fig,ax=plt.subplots(figsize=[10,5])

    for j,numb_pix in enumerate(numb_pix_array):

        tpi=1.5e-3
        Omega=np.pi/tpi
        t=np.linspace(0,20e-3,400)

        a=5.45e-6
        x=np.linspace(19*a,19*a+a*numb_pix,N)
        xx, yy = np.meshgrid(x,x)
        

        Detuning_array=get_lattice_detuning(xx,yy)
        # time_avg=Pe_time_avg_2d(slope/Omega,slope/Omega,np.abs(x[0])*2)

        # print(time_avg)

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
        # ax.plot(t*1e3,time_avg*np.ones(len(t)),color=palettable.colorbrewer.qualitative.Set1_6.hex_colors[j])

    ax.set_xlabel('time (ms)',fontsize=15)
    ax.set_ylabel('excitation fraction',fontsize=15)
    ax.legend(loc="lower right")
    ax.set_ylim([0,1])

    fig.savefig('excfraction_VS_time_potent_ring.pdf')
    # plt.plot(t,detuned,":o",label="detuned")  
    plt.show()

        