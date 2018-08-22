import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

P = [0, 0.9, 1.1, 1.4]

for name in glob.glob('temps_setup1_run1.csv'):
    temp_file = name

filename = temp_file
temp = np.loadtxt(filename)
#print(temp)

chi = np.zeros([len(P), len(temp)])

for run_num in range(4):
    for name in glob.glob('chi_*'+str(run_num+1)+'.csv'):
        chi_file = name
    filename = chi_file
    chi[run_num]= np.loadtxt(filename)

PT = []
for i in range(len(P)):
    for j in range(len(temp)):
        PT.append([P[i],temp[j]])

#print (chi)
x = np.arange(0, 1.4, 0.01)

chi2d = []
plt.figure()
for j in range(len(temp)):
    f = interpolate.interp1d(P,chi[:,j],kind="linear")
    chi2d.append(f(x))
    plt.plot(x, f(x))
plt.show()




y = temp#np.arange(0, 300, 0.5)
#
##plt.plot(x, f[0](x))
##plt.plot(x, f[100](x))
##plt.show()
#
X, Y = np.meshgrid(x, y)
##print(X)
##print(Y)
Z = chi2d
#print (chi_2Dinterp(0.0,270))
plt.pcolormesh(X, Y, Z, cmap='hsv') 
#
pp = plt.colorbar (orientation="vertical") 
pp.set_label("chi", fontname="Arial", fontsize=24) 
#plt.clim(0.06, 0.2)
plt.xlabel('Pressure', fontsize=24)
plt.ylabel('Temperature', fontsize=24)
#
plt.show()

#for j in range(len(temp)):
#    chi_T_fix = []
#    for i in range(len(P)):
#        chi_T_fix.append(chi[i][j])
#    chi_interp1d[i] = interpolate.interp1d(P, chi_T_fix)
#print(chi_interp1d[20](0))
#T = [6, 8, 10, 296, 300]
#chi=np.zeros([len(P),len(T)])
#for i in range(len(P)):
#    for j in range(len(T)):
#        chi[i,j] = i+j
#
#plt.imshow(chi_interp1d)
#plt.colorbar () # カラーバーの表示 
#plt.xlabel('Pressure')
#plt.ylabel('Temperature')
#plt.show()
