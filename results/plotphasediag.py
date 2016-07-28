import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import os

filelist=os.listdir(os.getcwd())

filename=sys.argv[1]

def theta(x):
    if x>0:
        return 1
    else:
        return 0

upperTreshold=0.75
grey=0.4

def colormap(x,y):
    if abs(x)<1e-4:
        rgbcol=(grey,grey,grey)
        if abs(x)<1e-12:
            rgbcol=(0,0,0)
    else:
        if abs(y)<0.1:
            rgbcol=(1,0,x)
        else:
            rgbcol=(1-x,0,1)
    return rgbcol

readOneD=False
writeOneD=False

x=[]
y=[]
c=[]
r=[]
fs=24
markersize=11.5

for file in filelist:
    if file[0:len(filename)]==filename and file[(len(file)-4):len(file)]=='.txt' and file[(len(file)-8):(len(file)-4)]!='oneD':
        print file
        diag=np.loadtxt(file,skiprows=1)
        diag=diag.transpose()
        x.extend(diag[0])
        y.extend(diag[1])
        c.extend(diag[2])
        #r.extend(diag[3])
plt.figure(figsize=(12,10))
plt.tick_params(labelsize=fs)
for i in range(0,len(x)):
    plt.plot(x[i],y[i],'o',color=(0,c[i],0),ms=markersize)
plt.xlabel('J',color=(170/255.0,0,45/255.0),fontsize=fs)
plt.ylabel('g',color=(25/255.0,170/255.0,70/255.0),fontsize=fs)
plt.ylim(ymin=.94,ymax=1.055)
plt.xlim(xmin=.94,xmax=1.055)
#plt.title('$\\frac{\\min(n_i^a)}{\\max(n_i^a)}$')
plt.savefig('thesis_plots/'+filename+'_zoom.pdf',bbox_inches='tight')
#plt.show()

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
