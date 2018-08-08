import numpy as np
from scipy.signal import argrelmax
from scipy.optimize import curve_fit
from scipy import interpolate
import matplotlib.pylab as plt

### factors of MPMS ####

LRF = 1.825 ##longitudinal regression factor
SCF = 7476.481 ##SQUID cal. factor

def sens_factor (r,g):
    if (r==0.0 and g==3.0):
        return 10.00
    if (r==0.0 and g==2.0):
        return 5.00
    if (r==0.0 and g==1.0):
        return 2.00
    if (r==0.0 and g==0.0):
        return 1.00
    if (r==1.0 and g==2.0):
        return 0.50
    if (r==1.0 and g==1.0):
        return 0.20
    if (r==1.0 and g==0.0):
        return 0.10
    if (r==2.0 and g==2.0):
        return 0.05
    if (r==2.0 and g==1.0):
        return 0.02
    if (r==2.0 and g==0.0):
        return 0.01
    if (r==3.0 and g==2.0):
        return 0.005
    if (r==3.0 and g==1.0):
        return 0.002
    if (r==3.0 and g==0.0):
        return 0.001

CF = 0.9125 #correlation factor


#### parameters of measurement #####
dp = 64 #data points per scan
sc = 3 #the number of scans per temperature
scanlength = 5.0 #cm
B = 10000#Oe
g = 0.000587#g
A = 6.941*2 + 192.217 + 15.9994*3#g/mol

#### parameters for program #####
ext = 0#4 #the number of temperature to extrpolate

delta = 0.001
#avg_delta_peak = 7.2-6.66-(0.55-0.5) #cm #peak hosei at T = 8 K

ini = [0.0,0.0,0.0,-scanlength*0.5]

weight1 = 1.0
weight2 = 10.0 # cms[i][j] < LL, cms[i][j] > UL -> 1, LL < cms[i][j] < UL -> weight
LL =   2.0#cm, lower limit of weight.   
UL =  4.6#cm, upper limit of weight. 
fit_weight_cm =  np.linspace(0.0, scanlength, num=dp) # = cms[i]
fit_weight = np.zeros(dp)
for i in range(dp):
    if fit_weight_cm[i] <= LL:
        fit_weight[i] = weight1
    if LL < fit_weight_cm[i] < UL:
        fit_weight[i] = weight2
    if fit_weight_cm[i] >= UL:
        fit_weight[i] = weight1



LLcut = 1.0#cm
ULcut = 4.8#cm

### file loading ###
diag = "180228_MT_1T_6-300K_betaLi2IrO3_b-axis_run1_AP.dc.diag"
raw = "180228_MT_1T_6-300K_betaLi2IrO3_b-axis_run1_AP.dc.raw"

diagBG = "1802011_MT_1T_6-300K_betaLi2IrO3_highP_bkgr.dc.diag_modified.csv"
skip_diagBG = 0#8*sc

rawBG = "1802011_MT_1T_6-300K_betaLi2IrO3_highP_bkgr.dc.raw_modified.csv"
skip_rawBG = 0#8*sc*dp


###Sensitivity factor of signal###
filename = diag
data = np.loadtxt(filename, delimiter=',', skiprows=27, usecols=(7,8))
#[0]-[26]:header
ranges = data[:,0]
gains = data[:,1]

SFs = np.zeros(ranges.size)
for i in range(int(ranges.size)):
    SFs[i] = sens_factor (ranges[i],gains[i]) #sensitivity factor


###Signal###
filename = raw
data = np.loadtxt(filename, delimiter=',', skiprows=31, usecols=(7,9))
#lines[0]-lines[30] : header 
#data = [[cm,v],
#         ...,
#          [cm,v]]
cms = data[:,0]# = 0th retsu of all gyou
vs = data[:,1]# = 1st retsu of all gyou
sharp_scans = cms.size /dp
sharp_temps = sharp_scans/sc
cms = cms.reshape(int(sharp_scans),dp)#cms = (scans, 48) matrix = [[48pt],
                            #                                ...,
                            #                               [48pt]]. cms[i] = i th list

vs = vs.reshape(int(sharp_scans), dp)
svs = np.zeros (cms.size)
svs = svs.reshape(int(sharp_scans), dp)
for i in range(int(sharp_scans)):
    svs[i] = vs[i]*LRF/SCF/SFs[i] #scaled voltage

 
###check temp stability of signal###
data = np.loadtxt(filename, delimiter=',', skiprows=31, usecols=(0,3,4))
time = data[:,0]
start_temps = data[:,1]
end_temps = data[:,2]
temp_stability = start_temps-end_temps

plt.figure()
plt.plot (time, temp_stability)
plt.savefig('temp_stability.eps')


start_temps = start_temps.reshape(int(sharp_temps),dp*sc)
temps = start_temps[:,0]
np.savetxt ('temps.csv',temps)


###Sensitivity factor of background###
filename = diagBG 
data = np.loadtxt(filename, delimiter=',', skiprows=27+skip_diagBG, usecols=(7,8))
ranges = data[:,0]
gains = data[:,1]

SFs = np.zeros(ranges.size)
for i in range(int(ranges.size)):
    SFs[i] = sens_factor (ranges[i],gains[i])


###background###
filename = rawBG 
data = np.loadtxt(filename, delimiter=',', skiprows=31+skip_rawBG, usecols=(7,9))
#lines[0]-lines[30] : header 
#data = [[cm,v],
#         ...,
#          [cm,v]]
cmsBG = data[:,0]# = 0th retsu of all gyou
vsBG = data[:,1]# = 1st retsu of all gyou
sharp_scansBG = cmsBG.size /dp
sharp_tempsBG = sharp_scansBG/sc
cmsBG = cmsBG.reshape(int(sharp_scansBG),dp)#cms = (scans, 48) matrix = [[48pt],
                            #                                ...,
                            #                               [48pt]]. cms[i] = i th list

vsBG = vsBG.reshape(int(sharp_scansBG), dp)
svsBG = np.zeros (cmsBG.size)
svsBG = svsBG.reshape(int(sharp_scansBG), dp)
for i in range(int(sharp_scansBG)):
    svsBG[i] = vsBG[i]*LRF/SCF/SFs[i] #scaled voltage


###check temp stability of backgorund###
data = np.loadtxt(filename, delimiter=',', skiprows=31+skip_rawBG, usecols=(0,3,4))
time = data[:,0]
start_temps = data[:,1]
end_temps = data[:,2]
temp_stability = start_temps-end_temps

plt.plot (time, temp_stability)
plt.savefig('temp_stabilityBG.eps')

start_temps = start_temps.reshape(int(sharp_tempsBG),dp*sc)
tempsBG = start_temps[:,0]
np.savetxt ('tempsBG.csv',tempsBG)

### remove unnecessray BG ###
while True:
    cmsBG = cmsBG.reshape(int(sharp_tempsBG),dp*sc)
    svsBG = svsBG.reshape(int(sharp_tempsBG), dp*sc)
    print(sharp_scansBG)
    print(sharp_tempsBG)
    for j in range(int(sharp_tempsBG)):
        if j == int(sharp_tempsBG)-1:
            break
        if abs(tempsBG[j]-tempsBG[j+1]) < 0.05:
            print(tempsBG[j])
            cmsBG = np.delete(cmsBG,j,0)
            svsBG = np.delete(svsBG,j,0)
            tempsBG = np.delete(tempsBG,j,0)
            break
    if j == int(sharp_tempsBG)-1:
            break
    cmsBG = cmsBG.reshape(int(sharp_scansBG)-sc,dp)
    svsBG = svsBG.reshape(int(sharp_scansBG)-sc, dp)
    sharp_scansBG = cmsBG.size /dp
    print(sharp_scansBG)
    sharp_tempsBG = sharp_scansBG/sc
    print(sharp_tempsBG)


##### extrapolate BG ####
cmsBG = cmsBG.reshape(int(sharp_tempsBG),dp*sc)
svsBG = svsBG.reshape(int(sharp_tempsBG), dp*sc)

##cmsBG_ext = np.zeros((int(sharp_tempsBG)+ext,dp*sc))
##svsBG_ext = np.zeros((int(sharp_tempsBG)+ext, dp*sc))
##tempsBG_ext = np.zeros(int(sharp_tempsBG)+ext)
##
##for j in range(int(sharp_tempsBG)+ext):
##    if j < int(sharp_tempsBG):
##        cmsBG_ext[j] = cmsBG[j]
##        svsBG_ext[j] = svsBG[j]
##        tempsBG_ext[j] = tempsBG[j]
##    if j > int(sharp_tempsBG)-1:
##        cmsBG_ext[j] = cmsBG[int(sharp_tempsBG)-1]
##        svsBG_ext[j] = svsBG[int(sharp_tempsBG)-1]
##    
##tempsBG_ext[int(sharp_tempsBG)] = 305.0
##tempsBG_ext[int(sharp_tempsBG)+1] = 310.0
##tempsBG_ext[int(sharp_tempsBG)+2] = 315.0
##tempsBG_ext[int(sharp_tempsBG)+3] = 320.0
##    
##cmsBG = cmsBG_ext
##svsBG = svsBG_ext
##tempsBG =tempsBG_ext

sharp_scansBG = cmsBG.size /dp
print(sharp_scansBG)
sharp_tempsBG = sharp_scansBG/sc
print(sharp_tempsBG)
print(tempsBG.size)

print(sharp_scans-sharp_scansBG)
print(sharp_temps-sharp_tempsBG)

cmsBG = cmsBG.reshape(int(sharp_scansBG),dp)
svsBG = svsBG.reshape(int(sharp_scansBG),dp)

### take average of signal ###
cms = cms.reshape(int(sharp_scans),dp)
svs = svs.reshape(int(sharp_scans),dp)
avg_cm = np.zeros((int(sharp_temps),dp))
avg_sv = np.zeros((int(sharp_temps),dp))
for j in range(int(sharp_temps)):
    for i in range(sc):
        avg_cm[j] += cms[i+sc*j]
        avg_sv[j] += svs[i+sc*j]
    avg_cm[j] = avg_cm[j]/sc
    avg_sv[j] = avg_sv[j]/sc


### take average of BG ###
cmsBG = cmsBG.reshape(int(sharp_scansBG),dp)
svsBG = svsBG.reshape(int(sharp_scansBG),dp)
avg_cmBG = np.zeros((int(sharp_tempsBG),dp))
avg_svBG = np.zeros((int(sharp_tempsBG),dp))
for j in range(int(sharp_tempsBG)):
    for i in range(sc):
        avg_cmBG[j] += cmsBG[i+sc*j]
        avg_svBG[j] += svsBG[i+sc*j]
    avg_cmBG[j] = avg_cmBG[j]/sc
    avg_svBG[j] = avg_svBG[j]/sc

np.savetxt('avg_cm.csv', avg_cm[0])
## avg_svs plot ###
plt.figure()
for j in range(int(sharp_temps)):
    plt.plot(avg_cm[j], avg_sv[j], 'o')
    np.savetxt('avg_sv_'+str(temps[j])+'K.csv', avg_sv[j])
filename = 'avg_svs.eps'
plt.savefig(filename)


## avg_svsBG plot ###
plt.figure()
for j in range(int(sharp_temps)):
    plt.plot(avg_cmBG[j], avg_svBG[j], '.')
    np.savetxt('avg_svBG_'+str(temps[j])+'K.csv', avg_svBG[j])
filename = 'avg_svsBG.eps'
plt.savefig(filename)

## netsvs_delta0 plot ###
netsv = np.zeros((int(sharp_tempsBG),dp))
plt.figure()
for j in range(int(sharp_temps)):
    netsv[j] = avg_sv[j]-avg_svBG[j]
    plt.plot(avg_cm[j], netsv[j], '+')
filename = 'netsvs_delta0.eps'
plt.savefig(filename)


#### fitting ####

def svfit(Z, X1, X2, X3, X4):
    R = 0.97 ##cm
    lmbda = 1.519 ##cm
    return X1 + X2*Z + X3*(2*((R**2 +(Z + X4)**2)**(-1.5))\
                           -(R**2 + (lmbda + Z + X4)**2)**(-1.5)\
                           -(R**2 + (-lmbda + Z + X4)**2)**(-1.5))



emu = np.zeros(int(sharp_temps))
emu_weight = np.zeros(int(sharp_temps))
emu_err = np.zeros(int(sharp_temps))
emu_err_percent = np.zeros(int(sharp_temps))

residual = np.zeros(int(sharp_temps))
residual_perc = np.zeros(int(sharp_temps))
sum_netsv = np.zeros(int(sharp_temps))

fit_weights = np.zeros((int(sharp_temps), dp))
for j in range(int(sharp_temps)):
    fit_weights[j] = fit_weight



for j in range(int(sharp_temps)):
    plt.figure()
    plt.plot(avg_cm[j], avg_sv[j], 'o')
    #plt.plot(avg_cmBG[j], avg_svBG[j], '.')
    
    avg_cmBG[j] = avg_cmBG[j]+delta
    plt.plot(avg_cmBG[j], avg_svBG[j], '.')
    if delta > 0: #truncate signal and BG to the data from delta to scanlength
        for k in range(avg_cm[j].size):
            if avg_cm[j][k] > delta:
                delpts = np.arange(0, k, 1)#0,1,...,k-1
                cm = np.delete(avg_cm[j], delpts)
                sv = np.delete(avg_sv[j], delpts)
                fit_weight0 = np.delete(fit_weights[j], delpts)
                break
        for k in range(avg_cmBG[j].size):
            if avg_cmBG[j][k] > scanlength:
                delpts = np.arange(k+1, avg_cmBG[j].size, 1)
                cmBG = np.delete(avg_cmBG[j], delpts)
                svBG = np.delete(avg_svBG[j], delpts)
                break
##    if delta < 0:
##        for k in range(avg_cm[j].size):
##            if avg_cm[j][k] >= delta:
##                delpts = np.arange(0, k, 1)#0,1,...,k-1
##                cm = np.delete(avg_cm[j], delpts)
##                netsv = np.delete(netsv[j], delpts)
##                fit_weight = np.delete(fit_weight, delpts)
##                break
##        for k in range(avg_cmBG[j].size):
##            if avg_cmBG[j][k] > scanlength:
##                delpts = np.arange(k, avg_cmBG[j].size, 1)
##                cmBG = np.delete(avg_cmBG[j], delpts)
##                netsvBG = np.delete(netsvBG[j], delpts)
##                break
   
    f_bg = interpolate.interp1d(cmBG, svBG, kind="cubic")
    netsv = sv - f_bg(cm)
    plt.plot(cm, netsv, '+')
        
    ## limit fitting region ##
    for k in range(cm.size):
        if cm[k] > LLcut:
            delpts = np.arange(0, k, 1)#0,1,...,k-1
            cm_lmt1 = np.delete(cm, delpts)
            netsv_lmt1 = np.delete(netsv, delpts)
            fit_weight_lmt1 = np.delete(fit_weight0, delpts)
            break
    for k in range(cm_lmt1.size):
        if cm_lmt1[k] > ULcut:
            delpts = np.arange(k, cm_lmt1.size, 1)
            cm_lmt2 = np.delete(cm_lmt1, delpts)
            netsv_lmt2 = np.delete(netsv_lmt1, delpts)
            fit_weight_lmt2 = np.delete(fit_weight_lmt1, delpts)
            break
        
    param, cov = curve_fit(svfit, cm_lmt2, netsv_lmt2, p0=ini, sigma=1/fit_weight_lmt2)
    ini=param
    error = np.sqrt(np.diag(cov))
    plt.plot (cm_lmt2, svfit(cm_lmt2, *param))
    filename = 'temp'+str(temps[j])+'K.eps'
    plt.savefig(filename)

    diff = netsv_lmt2 - svfit(cm_lmt2, *param)
    for k in range(diff.size):
        residual[j] += diff[k]**2
        sum_netsv[j] += netsv_lmt2[k]**2 
    residual_perc[j] = residual[j]/sum_netsv[j]*100
   
    
    emu[j] = param[2]/CF    
    emu_err[j] = error[2]/CF
    emu_err_percent[j] = emu_err[j]/emu[j]*100

np.savetxt('emu.csv', emu)
np.savetxt('emu_err.csv', emu_err)
np.savetxt('emu_err_percent.csv', emu_err_percent)
np.savetxt('residual.csv', residual)
np.savetxt('residual_perc.csv', residual_perc)

chi = emu/g*A/B
chi_err = emu_err/g*A/B
chi_inv = 1/chi
chi_inv_err = chi_err/chi/chi

np.savetxt('chi.csv', chi)
np.savetxt('chi_err.csv', chi_err)
np.savetxt('chi_inv.csv', chi_inv)
np.savetxt('chi_inv_err.csv', chi_inv_err)

plt.figure()
plt.plot(temps, emu, '+')
filename = 'emu.eps'
plt.savefig(filename)


plt.figure()
plt.plot(temps, chi, '+')
filename = 'chi.eps'
plt.savefig(filename)


plt.figure()
plt.plot(temps, residual_perc, '+')
filename = 'residual_perc.eps'
plt.savefig(filename)







