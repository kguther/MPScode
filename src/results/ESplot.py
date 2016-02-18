import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import math

def getSpec(x):
    return -2*math.log(x)

filelist=os.listdir(os.getcwd())
if 'plots' not in filelist:
    os.mkdir('plots')

taskname=sys.argv[1]
site=int(sys.argv[2])

gather=[]
labels=[]

for filename in filelist:
    if filename[0:(len(taskname)+6)]==taskname+'_point':
        if filename[(len(filename)-6):(len(filename)-4)]=='ES':
            data=[]
            print filename
            with open(filename,'r') as readData:
                for line in readData:
                    buf=line.split('\t')
                    if len(buf)>(site+1):
                       data.append(float(buf[site]))
            gather.append(data[1:len(data)])
        else:
            with open(filename,'r') as lbf:
                print filename
                lbf.readline()
                qnsr=lbf.readline()
                hparsr=lbf.readline()
                qnsr=qnsr[0:(len(qnsr)-1)]
                hparsr=hparsr[0:(len(hparsr)-1)]
                qns=qnsr.split('\t')
                hpars=hparsr.split('\t')
                labels.append(qns+hpars)

plt.figure()
for i in range(0,len(gather)):
    xdata=range(0,len(gather[i]))
    ydata=map(getSpec,gather[i])
    ydata.sort()
    print ydata
    esplot,=plt.plot(xdata,ydata,'o')
    legendstring='L='+labels[i][0]+' N='+labels[i][1]+' $\\alpha$='+labels[i][2]+' J='+labels[i][3]+' g='+labels[i][4]
    esplot.set_label(legendstring)
plt.xlabel('j')
plt.ylabel('$E_j$')
plt.legend(loc=1)
#plt.savefig('plots/'+taskname+'_entanglement_spectrum.pdf')
plt.close()
