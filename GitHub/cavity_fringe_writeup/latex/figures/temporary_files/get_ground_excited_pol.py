import numpy as np
import srlab.plot.pyplot as plt
import matplotlib

# matplotlib.rcParams.update({'font.size': 30})
import matplotlib.gridspec as gridspec

# gs = gridspec.GridSpec(1,7,hspace=0.05,wspace=0.5, bottom=0.3,
#                        left=0.02, right=0.95, width_ratios=[1,1,1,1,1,1,0.1])

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# matplotlib.rcParams.update({' font.size': 30})
# 

def linear_interpl(x,y,dx,dy,a):

    new_value=x+(y-x)*a
    new_value_std=np.sqrt(((1-a)*dx)**2+(a*dy)**2)
    return new_value,new_value_std


def clock_detuning_with_excited_trap(trap_freq,recoil_freq,ratio_ge,dratio_ge):
    detuning=(trap_freq**2)*(ratio_ge-1)/(2*recoil_freq)+trap_freq*(1-(ratio_ge)**0.5)
    ddetuning=dratio_ge*((trap_freq**2)/(2*recoil_freq)-trap_freq/2/np.sqrt(ratio))
    return detuning,ddetuning

if __name__ == '__main__':
    v,dv=linear_interpl(261.09,260.91,1.16,1.16,0.332)
    print(v,dv)

    #polarizability at 914.332 nm

    alpha_g = v
    dalpha_g = dv
    alpha_e = 220.8
    dalpha_e = 2.3

    ratio=alpha_g/alpha_e
    dratio = ((dalpha_g/alpha_e)**2+(((alpha_g/(alpha_e**2))*dalpha_e)**2))**0.5

    ratio_measured=1.179
    dratio_measured=0.002

    # print(ratio_measured, dratio_measu)

    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=[(297/ 25.4), 210/2 / 25.4], sharex=False, sharey=False)
    fig.tight_layout(pad=3, w_pad=4.0, h_pad=10)

    axs[0].errorbar(['M.'],ratio,dratio,fmt='o',markersize=8,capsize=20)
    axs[0].errorbar(['measured'],ratio_measured,dratio_measured,fmt='o',markersize=8,capsize=20)
    axs[0].set_xlim([-1,2])
    axs[0].set_ylabel(r'ratio $\alpha_g/\alpha_e$')

    recoil_freq=2.712
    trap_freq_high=106.56+recoil_freq
    trap_freq_low=74.74+recoil_freq

    detuning_high,ddetuning_high=clock_detuning_with_excited_trap(trap_freq_high,recoil_freq,ratio_measured,dratio_measured)
    detuning_low,ddetuning_low=clock_detuning_with_excited_trap(trap_freq_low,recoil_freq,ratio_measured,dratio_measured)

    Mdetuning_high,Mddetuning_high=clock_detuning_with_excited_trap(trap_freq_high,recoil_freq,ratio,dratio)
    Mdetuning_low,Mddetuning_low=clock_detuning_with_excited_trap(trap_freq_low,recoil_freq,ratio,dratio)

    print(detuning_high,detuning_low)

    axs[1].errorbar(['comb'],399,1,fmt='o',markersize=8,capsize=20)
    axs[1].errorbar(['measured'],detuning_high,ddetuning_high,fmt='o',markersize=8,capsize=20)
    axs[1].errorbar(['M.'],Mdetuning_high,Mddetuning_high,fmt='o',markersize=8,capsize=20)

    # axs[1].set_xlim([-0.5,1.6])
    axs[1].set_ylabel('detuning (kHz)')

    axs[2].errorbar(['comb'],194,1,fmt='o',markersize=8,capsize=20)
    axs[2].errorbar(['measured'],detuning_low,ddetuning_low,fmt='o',markersize=8,capsize=20)
    axs[2].errorbar(['M.'],Mdetuning_low,Mddetuning_low,fmt='o',markersize=8,capsize=20)

    # axs[2].set_xlim([-0.5,1.6])
    axs[2].set_ylabel('detuning (kHz)')
    plt.savefig('ratio_comparison_plot.pdf')


    # plt.show() 