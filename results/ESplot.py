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

cfiles=[]
site=int(sys.argv[1])
if len(sys.argv)==3:
    taskname=sys.argv[2]
else:
    taskname=''
    for i in range(2,len(sys.argv)):
        cfiles.append(sys.argv[i])

gather=[]
labels=[]

for filename in filelist:
    if filename[0:(len(taskname))]==taskname or filename[0:(len(filename)-7)] in cfiles or filename[0:(len(filename)-4)] in cfiles:
        if filename[(len(filename)-6):(len(filename)-4)]=='ES':
            data=[]
            print filename
            with open(filename,'r') as readData:
                lineindex=0
                for line in readData:
                    if lineindex==site+1:
                        buf=(line.strip()).split('\t')
                        data=map(float,buf)
                    lineindex+=1
            gather.append(data)
        else:
            with open(filename,'r') as lbf:
                print filename
                lbf.readline()
                lbf.readline()
                hparsr=lbf.readline()
                hparsr=hparsr[0:(len(hparsr)-1)]
                hpars=hparsr.split('\t')
                labels.append(hpars)

plt.figure()
for i in range(0,len(gather)):
    xdata=range(0,len(gather[i]))
    ydata=map(getSpec,gather[i])
    ydata.sort()
    print ydata
    esplot,=plt.plot(xdata,ydata,'o')
    legendstring='L='+labels[i][0]+' N='+labels[i][1]+' $\\alpha$='+labels[i][2]+' J='+labels[i][3]+' g='+labels[i][4]+ 'W='+labels[i][5]
    esplot.set_label(legendstring)
plt.xlabel('j')
plt.ylabel('$E_j$')
plt.legend(loc=0,fontsize='xx-small',numpoints=1)
#plt.savefig('plots/'+taskname+'_entanglement_spectrum.pdf')
plt.show()
plt.close()
