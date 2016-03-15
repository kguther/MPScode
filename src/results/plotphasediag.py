import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import os

filelist=os.listdir(os.getcwd())

filename='phasediagram'

def colormap(x):
    rgbcol=(0,abs(x),abs(x/2))
    return rgbcol

readOneD=False
writeOneD=True

x=[]
y=[]
c=[]
fs=24

for file in filelist:
    if file[0:len(filename)]==filename and file[(len(file)-4):len(file)]=='.txt' and file[(len(file)-8):(len(file)-4)]!='oneD':
        print file
        diag=np.loadtxt(file,skiprows=1)
        diag=diag.transpose()
        x.extend(diag[0])
        y.extend(diag[1])
        c.extend(diag[2])
plt.figure()
plt.tick_params(labelsize=fs)
for i in range(0,len(x)):
    plt.plot(x[i],y[i],'o',color=colormap(c[i]))
plt.xlabel('J',color=(170/255.0,0,45/255.0),fontsize=fs)
plt.ylabel('g',color=(25/255.0,170/255.0,70/255.0),fontsize=fs)
plt.ylim(ymin=0,ymax=1.5)
plt.xlim(xmin=0,xmax=1.5)
#plt.title('$\\frac{\\min(n_i^a)}{\\max(n_i^a)}$')
plt.savefig(filename+'_large_raw.eps',bbox_inches='tight')

if writeOneD:
    alpha=[]
    p=[]
    for i in range(0,len(x)):
        J=x[i]
        g=y[i]
        a=math.atan2((g-1),(J-1))
        if alpha.count(a)==0:
            alpha.append(a)
            p.append(c[i])

    with open(filename[0:(len(filename)-4)]+'_oneD.txt','w') as pd:
        for i in range(0,len(alpha)):
            point=str(alpha[i])+'\t'+str(p[i])+'\n'
            pd.write(point)

if readOneD:
    diag=np.loadtxt(filename[0:(len(filename)-4)]+'_oneD.txt')
    diag=diag.transpose()
    a=diag[0]
    p=diag[1]
    plt.figure()
    plt.plot(a,p,'o')
    plt.ylim([-.05,1.05])
    plt.xlabel('$\\arctan\\left(\\frac{g}{J} \\right)$')
    plt.ylabel('$\\frac{\\min(n_i^a)}{\\max(n_i^a)}$')
    plt.savefig(filename[0:(len(filename)-4)]+'_oneD.txt')
