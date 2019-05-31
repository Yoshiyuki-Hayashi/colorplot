import os
import sys
import numpy as np
from matplotlib import pyplot as pl
from scipy.interpolate import UnivariateSpline, InterpolatedUnivariateSpline
from matplotlib.colors import LogNorm
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from mpms2 import *

path = os.getcwd()


'''
Usage: Find index of element of array closest to value.
'''
def find_idx(array,value):
#    print(array)
#    print(value)
#    print(array - float(value))
    idx = np.abs(array - float(value)).argmin() #I'm not sure why I had to convert value to a float but it was throwing an error for like 4 hours and im just glad I fixed it and can go to bed
    return idx



"""
[filename, P, erP, Tmag, erTmag, Td, erTd, a]
"""
filenames = [['/AP_straw/AP_straw_.txt', 0.0, 0.0, 38.0, 2.0, 0.0, 0.0, 0.0, 'indianred', 'o'],\
             ['/setup2/re-analyze/run2_.txt', 3.068e-01, 1.376e-01, 38.0, 2.0, 0.0, 0.0, 0.001168, 'coral','s'],\
             ['/setup1/re-analyze/run2_.txt', 8.520e-01, 1.888e-02, 38.0, 2.0, 0.0, 0.0, 0.0, 'goldenrod','D'],\
             ['/setup1/re-analyze/run3_.txt', 1.002, 1.797e-02, 38.0, 2.0, 98.0, 4.0, 0.0, 'darkgoldenrod', '^'],\
             ['/setup1/re-analyze/run4_.txt', 1.386, 2.424e-02, 36.0, 2.0, 154.0, 4.0, -0.000842, 'olivedrab','v'],\
             ['/setup6-1903/run1_ZFC_.txt', 1.608, 2.469e-02, 0.0, 0.0, 182.0, 4.0, 0.0, 'green', '>'],\
             #['/setup6-1903/run1_FC_.txt',1.608, 2.469e-02, 0.0, 0.0, 170.0, 4.0],\
             #['/setup6-1903/run2_ZFC_.txt', 2.241, 2.876e-02, 0.0, 0.0, 218.0, 4.0, 0.0, 'lightseagreen', '<'],\
             ['/setup6-1903/run2_FC_.txt', 2.241, 2.876e-02, 0.0, 0.0, 217.0, 4.0, 0.0, 'lightseagreen', '<'],\
             ['/setup6-1903/run3_ZFC_.txt', 3.032, 3.250e-02, 0.0, 0.0, 292.0, 4.0, 0.0, 'royalblue', 'd'],\
             #['/setup6-1903/run3_FC_.txt', 3.032, 3.250e-02, 0.0, 0.0, 300.0, 4.0, 0.0]
             ]
NMR_Tds = [[1.3, 0.0, 90.0, 0.0],\
           [2.6, 0.0, 230.0, 0.0],\
           [3.5, 0.0, 292.0, 0.0],\
           [3.9, 0.0, 315.0, 0.0]]

NMR_Tds = np.array(NMR_Tds)

#synthesis data plot3
fig = pl.figure(figsize=(6,8))

synthesis_impurity = [['synthesis/171212_Li2IrO3_S42_crystal1/171212_MT_betaLi2IrO3_s42_crystal1.dc.dat', 0.50, 2.26],\
             ['synthesis/171213_Li2IrO3_S42_hishi2/171213_MT_1T_5-1.6K_betaLi2IrO3_S42_hishi2_2.dc.dat', 0.43, 1.12],\
             ['synthesis/171213_Li2IrO3_S44_hishi1/171213_MT_1T_5-1.6K_betaLi2IrO3_S44_hishi1.dc.dat', 0.44, 2.17],\
             ['synthesis/171214_Li2IrO3_S44_hishi2/171214_MT_1T_5-1.6K_betaLi2IrO3_S44_hishi2_2.dc.dat', 0.4, 1.055],\
             ['synthesis/171219_Li2IrO3_S47_hishi1/171218_MT_1T_5-1.6K_betaLi2IrO3_S47_hishi1_2.dc.dat', 0.69, 1.008],\
             ['synthesis/180518_Li2IrO3_S47_hishi2/180518_MT_betaLi2IrO3_s47hishi2.dc.dat', 0.54, 1.01]]




#synthesis data plot2
fig = pl.figure(figsize=(6,8))
path = os.getcwd()
synthesis_Tmag = [['synthesis/171116_MT_1T_5-300K_betaLi2IrO3_S14.dc.dat', 'indianred', 'o', 'No.14'],\
             ['synthesis/171123_MT_1T_5-300K_betaLi2IrO3_S24.dc.dat', 'coral', 's', 'No.24'],\
             ['synthesis/171124_MT_1T_5-300K_betaLi2IrO3_S27_hishi.dc.dat', 'goldenred', 'D', 'No.27']]

for syn in synthesis_Tmag:
    dataset, name = load_MPMS_dc(syn[0], 1.0, 1000.0)#emu = erg/Oe at 1 T
    dataset = np.array(dataset)
    temps = dataset[:,1]
    Ms = dataset[:,2]
    idx=find_idx(temps, 60.0)
    Ms_calib = Ms/Ms[idx]
    pl.plot(temps, Ms_calib, color=syn[1], linestyle='-', marker=syn[2], label=syn[3])
pl.set_xlim(0.0, 60.0)
pl.set_xlabel('$\it{T}$ (K)')
pl.set_ylabel('$\it{M}$(arb.units)')
pl.legend(loc='upper right')

pl.tight_layout()

pl.show()

"""
#synthesis data plot1
fig = pl.figure(figsize=(6,8))

batches_Li2CO3 = ['No.42', 'No.47']
batch_LiOHH2O =['No.44']
ax_1 = fig.add_subplot(2,1,1)
ratios_Li2CO3 = [1.12, 1.008]
ratio_LiOHH2O =[2.17]
ax_1.plot(batches_Li2CO3, ratios_Li2CO3, 'o-', ms=10, color='indianred', label='Li$_2$CO$_3$')
ax_1.plot(batch_LiOHH2O, ratio_LiOHH2O, '^', ms=10, color='royalblue', label='LiOH$\cdot$H$_2$O')
ax_1.set_ylabel('$\it{M}(1.8\mathrm{K})$/$\it{M}(10\mathrm{K})$', fontsize=25)
ax_1.tick_params(labelbottom=False, direction='inout', top=True, right=True, length=8, labelsize=25)
ax_1.legend(loc='upper left', frameon=False, fontsize=20)
arrow_dict=dict(width=1, shrink=0.1, color='black')#arrowstyle='->',
ax_1.annotate(s='', xy=(batches_Li2CO3[1], ratios_Li2CO3[1]),xycoords='data', xytext=(batches_Li2CO3[1], ratios_Li2CO3[1]+0.25), textcoords='data', arrowprops = arrow_dict, va='center', ha='left')
ax_1.text(batches_Li2CO3[1], ratios_Li2CO3[1]+0.26, '{}'.format(ratios_Li2CO3[1]), fontsize=25)


ax_2 = fig.add_subplot(2,1,2)
sizes_Li2CO3 = [0.43, 0.69]
size_LiOHH2O = [0.44]
ax_2.plot(batches_Li2CO3, sizes_Li2CO3, 'o-', ms=10, color='indianred')
ax_2.plot(batch_LiOHH2O, size_LiOHH2O, '^', ms=10, color='royalblue')
ax_2.set_ylabel('length along b-axis (mm)', fontsize=25)
ax_2.tick_params(direction='inout', top=True, right=True, length=8, labelsize=25)

pl.tight_layout()

pl.show()


#chi plot
fig = pl.figure(figsize=(7,9))
G = GridSpec(nrows=4, ncols=1)
ax_1 = pl.subplot(G[0:3])
ax_2 = pl.subplot(G[3])
left, bottom, width, height = [0.5, 0.65, 0.4, 0.3]
ax_inset = fig.add_axes([left, bottom, width, height])

#for f in filenames:
#    data = np.loadtxt(path + f[0])    
#    temps = data[:,0]
#    chis = data[:,1]
#    ax_1.plot(temps, chis,'o-')
#    ax_1.set_ylim(-0.003,0.01)

Ps = []
erPs = []
temp_cut = 6.0
chis_cut = []
arrow_dict=dict(arrowstyle='->', color='black')
for f in filenames:
    Ps.append(f[1])
    erPs.append(f[2])
    data = np.loadtxt(path + f[0])    
    temps = data[:,0]
    chis = data[:,1]-f[7]
    for i, temp in enumerate(temps):
        if np.abs(temp - temp_cut) < 0.1:
            chis_cut.append(chis[i])
        if np.abs(temp - f[5]) < 2.0:
            ax_2.annotate(s='', xy=(f[5], chis[i]),xycoords='data', xytext=(f[5], chis[i]+0.003), textcoords='data', arrowprops = arrow_dict, va='center', ha='left')
            ax_2.text(f[5]-8.0, chis[i]+0.0035, '$\it{T}_\mathrm{d}$', fontsize=14)
    ax_1.plot(temps, chis, ms=4, linestyle='-', color=f[8], marker=f[9], label='{:.2f} GPa'.format(f[1]))
    ax_2.plot(temps, chis, ms=4, linestyle='-', color=f[8], marker=f[9])

if int(len(chis_cut)) != int(len(Ps)):
    sys.exit('Inappropriate choice of temp_cut={}K!'.format(temp_cut))
ax_1.set_xlim(0.0, 300.0)
ax_1.set_ylim(0.0,0.09)
ax_1.set_ylabel('$\it{M}$/$\it{B}$ (emu/mol)', fontsize=20)
ax_1.tick_params(labelbottom=False, direction='inout', top=True, right=True, length=8, labelsize=20)
ax_1.legend(bbox_to_anchor=(1,0.4), loc='upper right', borderaxespad=0.5, fontsize=12, frameon=False)
ax_1.text(120.0, 0.025, '$\it{B}$ = 1 T // b-axis', fontsize=12)
ax_1.axvline(38.0, linestyle='dashed', linewidth=1, color='black')
ax_1.text(39.0, 0.085, '$\it{T}_\mathrm{mag}$', fontsize=14)

ax_2.set_xlim(0.0, 300.0)
ax_2.set_ylim(0.0,0.01)
ax_2.set_xlabel('$\it{T}$ (K)', fontsize=20)
ax_2.set_ylabel('$\it{M}$/$\it{B}$ (emu/mol)', fontsize=20)
ax_2.tick_params(direction='inout', top=True, right=True, length=8, labelsize=20)
ax_2.axvline(38.0, linestyle='dashed', linewidth=1, color='black')

ax_inset.errorbar(Ps, chis_cut, xerr=erPs, capsize=3, color='indianred', marker='o', label='$\it{T} = 6 K$')
ax_inset.set_xlim(0.0, np.max(Ps + erPs))
ax_inset.set_ylim(0.0, 0.09)
ax_inset.set_xlabel('$\it{P}$ (GPa)', fontsize=15)
ax_inset.set_ylabel('$\it{M}$/$\it{B}$ (emu/mol)', fontsize=15)
ax_inset.tick_params(direction='inout', top=True, right=True, length=6, labelsize=15)
ax_inset.legend(loc='upper right', fontsize=12, frameon=False)


pl.tight_layout()

pl.show()

fig = pl.figure(figsize=(11,9))
Tmags = []
erTmags = []
Tds = []
erTds = []
chiss_T_interp = []
dT = 0.1
temp_min = 2.0
temp_max = 300.0
for f in filenames:
    Tmags.append(f[3])
    erTmags.append(f[4])
    Tds.append(f[5])
    erTds.append(f[6])
    data = np.loadtxt(path + f[0])    
    temps = data[:,0]
    #Interpolate chi about temperature at each pressure
    temps_interp = np.arange(temp_min, temp_max, dT)
    chis = data[:,1]-f[7]
    spl_T = InterpolatedUnivariateSpline(temps[np.argsort(temps)], chis[np.argsort(temps)], k=1)
    chiss_T_interp.append(spl_T(temps_interp))
chiss_T_interp = np.array(chiss_T_interp)

#Interpolate pressure
dP = 0.01
P_interp = np.arange(Ps[0], Ps[-1], dP)

#Interpolate chi about pressure at each temperature
chiss_PT_interp = []
for j in range(len(temps_interp)):
    spl_P = InterpolatedUnivariateSpline(Ps, chiss_T_interp[:,j], k=1)
    chiss_PT_interp.append(spl_P(P_interp))
chiss_PT_interp = np.array(chiss_PT_interp)

#colorplot
#ax = fig.add_subplot(1,2,2)
ax = fig.add_subplot(1,1,1)


X, Y = np.meshgrid(P_interp, temps_interp)
Z = chiss_PT_interp

pl.pcolormesh(X, Y, Z, cmap='RdYlBu_r',  norm=LogNorm()) 
pp = pl.colorbar (orientation="vertical") 
#pp.set_clim(vmin=0.0005, vmax=0.1)
pp.set_label("$\it{M}/\it{B}$ (emu/mol)", fontname="Arial", fontsize=20) 
ax.set_xlabel('$\it{P}$ (GPa)', fontsize=20)
ax.set_ylabel('$\it{T}$ (K)', fontsize=20)
ax.tick_params(labelsize=20)


Ps_for_Tmag = []
erPs_for_Tmag = []
Tmags_nonzero = []
erTmags_nonzero = []
for i, Tmag in enumerate(Tmags):
    if Tmag != 0.0:
        Tmags_nonzero.append(Tmag)
        erTmags_nonzero.append(erTmags[i])
        Ps_for_Tmag.append(Ps[i])
        erPs_for_Tmag.append(erPs[i])
ax.errorbar(Ps_for_Tmag, Tmags_nonzero, xerr=erPs_for_Tmag, yerr=erTmags_nonzero, capsize=3, fmt='s', color='black', label='Tmag from M/B')

Ps_for_Td = []
erPs_for_Td = []
Tds_nonzero = []
erTds_nonzero = []
for i, Td in enumerate(Tds):
    if Td != 0.0:
        Tds_nonzero.append(Td)
        erTds_nonzero.append(erTds[i])
        Ps_for_Td.append(Ps[i])
        erPs_for_Td.append(erPs[i])
ax.errorbar(Ps_for_Td, Tds_nonzero, xerr=erPs_for_Td, yerr=erTds_nonzero, capsize=3, fmt='o', color='black', label='Td from M/B')
ax.plot(NMR_Tds[:,0], NMR_Tds[:,2], '^', color='black', label='Td from NMR')


#ax.legend(loc='upper left')
pl.show()
"""