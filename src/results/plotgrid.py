import matplotlib.pyplot as plt
import numpy as np
import sys
import os

filename=sys.argv[1]
filelist=os.listdir(os.getcwd())

fs=24
k=9

def colormap(data):
    return (0,(data[k]-minE)/(maxE-minE),0)
    #return (0,data[6],0)

if filename in filelist:
    gridData=np.loadtxt(filename,skiprows=1)
    maxE=max((gridData.transpose())[k])
    minE=min((gridData.transpose())[k])
    c=colormap(gridData)
    plt.figure()
    plt.tick_params(labelsize=fs)
    for i in range(0,len(gridData)):
        plt.plot(gridData[i][2],gridData[i][3],'o',color=colormap(gridData[i]))
    plt.xlabel('J',color=(170/255.0,0,45/255.0),fontsize=fs)
    plt.ylabel('g',color=(25/255.0,170/255.0,70/255.0),fontsize=fs)
    plt.show()
