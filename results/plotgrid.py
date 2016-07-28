import matplotlib.pyplot as plt
import numpy as np
import sys
import os

filename=sys.argv[1]
filelist=os.listdir(os.getcwd())

fs=24

def revivalScaling(rev,cmin):
    t=0.069
    a=0.7
    if cmin>0.25:
        return 0.0
    elif rev<t:
        return a/t*rev
    else:
        return a+(1-a)*np.sin(np.pi/2*(rev-t)/(1-t))

def cdwScaling(cdw):
    if cdw<0:
        return 0
    elif cdw>0.5:
        return 1
    else:
        return 2*cdw

def colormap(data):
    t=1e-6
    b=0.7
    #return (0,(data[k]-minE)/(maxE-minE),0)
    if data[9]>0.25 and data[5]<1e-10:
        return (0,1,0)
    elif data[8]<0.1 and data[9]<0.3:
        return (1,0,1)
    else:
        return (np.exp(-data[4]),revivalScaling(data[6],data[7]),cdwScaling(data[9]))
    #return (np.exp(-data[4]),0,0)

if filename in filelist:
    gridData=np.loadtxt(filename,skiprows=1)
    plt.figure(figsize=(10,9))
    plt.tick_params(labelsize=fs)
    for i in range(0,len(gridData)):
        plt.plot(gridData[i][2],gridData[i][3],'o',color=colormap(gridData[i]))
    plt.xlabel('J',color=(170/255.0,0,45/255.0),fontsize=fs)
    plt.ylabel('g',color=(25/255.0,170/255.0,70/255.0),fontsize=fs)
    #plt.show()
    plt.savefig('thesis_plots/perturbativePD.pdf')
