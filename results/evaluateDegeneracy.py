import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import scipy as sy
import scipy.optimize as so

filelist=os.listdir(os.getcwd())

taskname=sys.argv[1]
J=sys.argv[2]
g=sys.argv[3]

def f(x,a,n):
    return a/(x**n)


deltaE=[]
errors=[]
positions=[]
for filename in filelist:
    if filename[0:len(taskname)]==taskname and filename[(len(filename)-6):(len(filename)-4)]!='ES' and filename[(len(filename)-11):len(filename)-4]!='state_2':
            refName=filename.split('_J_')[1]
            getPars=refName.split('_g_')
            JVal=getPars[0]
            gVal=getPars[1].split('.')[0]
            if J==JVal and g==gVal:
                secondName=filename.split('.')[0]+'_state_2.txt'
                with open(filename) as readEnergy:
                    readEnergy.readline()
                    readEnergy.readline()
                    eList=readEnergy.readline().split('\t')
                    energyGS=(float(eList[6]))
                    errGS=(float(eList[7]))
                with open(secondName) as readSEnergy:
                    readSEnergy.readline()
                    readSEnergy.readline()
                    eList=readSEnergy.readline().split('\t')
                    energyES=(float(eList[6]))
                    errES=(float(eList[7]))                    
                deltaE.append(energyES-energyGS)
                errors.append((errES+errGS)/2.0)    
                getPos=filename.split('_N_')[0]
                pos=float(getPos.split('L_')[1])
                positions.append(pos)

fs=32
ls=28                

p0=sy.array([0.01,1])
fpars, acc=so.curve_fit(f,positions,deltaE,p0)

plt.figure(figsize=(14,10))
plt.tick_params(labelsize=ls)
plt.semilogy(positions,deltaE,'o',ms=11.5)
#sort positions
plt.semilogy(positions,f(positions,fpars[0],fpars[1]),'k')
print fpars
plt.xlabel('L',fontsize=fs)
plt.ylabel('Energy gap',fontsize=fs)
#plt.savefig('thesis_plots/top_gap_scaling.pdf')
plt.show()
