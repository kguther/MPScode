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

def getkey(x):
    return x[0]

errors=[]
baseData=[]
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
            errors.append((errES+errGS)/2.0)    
            getPos=filename.split('_N_')[0]
            pos=float(getPos.split('L_')[1])
            print pos
            print energyES
            print energyGS
            baseData.append([pos,energyES-energyGS])
fs=32
ls=28         
marksize=14

x0=33
x1=110
p0=sy.array([0.01,1])
colors=[(52.0/255.0,137.0/255.0,197.0/255.0), (137.0/255.0,199.0/255.0,58.0/255.0)]
plt.figure(figsize=(12,10))
plt.tick_params(labelsize=ls)
for i in [0]:
    sorter=sorted(baseData,key=getkey)
    positions=[x[0] for x in sorter]
    deltaE=[x[1] for x in sorter]
    fpars, acc=so.curve_fit(f,positions,deltaE,p0)
    fplot,=plt.loglog(range(x0-1,x1+4),f(range(x0-1,x1+4),fpars[0],fpars[1]),color='k')
    cplot,=plt.plot(positions[::2],deltaE[::2],'o',ms=marksize,color=colors[0])
    #print fpars
plt.xlabel('$L$',fontsize=fs)
plt.ylabel('Energy gap $\Delta$',fontsize=fs)
plt.xlim(xmin=x0,xmax=x1)
plt.ylim(ymin=0.015,ymax=0.055)
yt=[0.02,0.03,0.04,0.05]
plt.yticks(yt,map(str,yt))
xt=[40,50,60,70,80,90,100]
plt.xticks(xt,map(str,xt))
cplot.set_label('numerical data')
fplot.set_label('$\Delta=$'+str(fpars[0])[:6]+'$\cdot\,L^{-'+str(fpars[1])[:6]+'}$') 
plt.legend(loc=1,fontsize=fs,numpoints=1)
plt.savefig('../../draft/plots/top_gap_scaling.pdf')
plt.show()
