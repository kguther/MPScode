import colorsys
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

filename=sys.argv[1]
filelist=os.listdir(os.getcwd())

fs=18

def revivalScaling(rev,cmin):
    t=0.1
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
    elif cdw>0.6:
        return 1
    else:
        return 1/(0.6)*cdw

def accEv(x):
    return np.exp(-30*np.sqrt(x))

def colormap(data):
    t=1e-6
    b=0.7
    #check for phase separation: density fluctiations without CDW
    if data[9]>0.25 and data[5]<0.02:
        colrgb=(0,1,0)
    #check for interwire phase separation: 
    elif data[8]<0.2 and data[9]<0.3:
        colrgb=(1,0,1)
    else:
        colrgb=(np.exp(-data[4]),revivalScaling(data[6],data[7]),cdwScaling(data[9]))
    h,s,v=colorsys.rgb_to_hsv(*colrgb)
    s=s*accEv(data[11])
    return colorsys.hsv_to_rgb(h,s,v)
    #return (np.exp(-data[4]),0,0)


markerwidth=15
if filename in filelist:
    gridData=np.loadtxt(filename,skiprows=1)
    plt.figure(figsize=(18,13))
    plt.tick_params(labelsize=fs)
    for i in range(0,len(gridData)):
        plt.plot(gridData[i][2],gridData[i][3],'o',ms=markerwidth,color=colormap(gridData[i]))
    #plt.xlim(xmin=-4.9,xmax=2.15)
    #plt.ylim(ymin=-6,ymax=6.3)
    plt.xlabel('J',fontsize=fs)
    plt.ylabel('g',fontsize=fs)
    plt.savefig('thesis_plots/perturbativePD_raw.pdf',bbox_inches='tight')
    #plt.show()

plt.close()
