#import matplotlib 
#from matplotlib import rc

#rc('text',usetex=True)
#rc('text.latex',preamble='\usepackage{color} \definecolor{jred}{RGB}{170,0,45} \definecolor{gviolet}{RGB}{125,0,125}')


import numpy as np
import matplotlib.pyplot as plt
import os
import sys
labellist=['$\\left|\\langle a^\dagger_i a_0^{} \\rangle \\right|$','$\\left|\\langle a^\dagger_i b^{}_i b_0^{\dagger} a_0 \\rangle \\right|$','$\\left|\\langle n^{a}_i n_0^{a} \\rangle \\right|$','$\\left|\\langle n^{a}_i n_0^{b} \\rangle \\right|$','$\\left|\\langle n^{a}_i \\rangle \\right|$','$\\left|\\langle a^\dagger_i b^\dagger_i a_0^{} b_0^{} \\rangle \\right|$','$\\left|\\langle \\right|\\rangle$','S','$\\left|\\langle n^{a}_i n^{b}_i \\rangle\\right|$','$\\left|\\langle n^{b}_i\\rangle\\right|$','$\\langle a_i^{\dagger} a_{i+1}^{\dagger} a_0 a_1 \\rangle$','other']

def tasknum(n):
    if n=="Intrachain correlation" or n=="Bulk correlation function":
        taskindex=0
    elif n=="Interchain hopping correlation" or n=="Bulk interchain hopping correlation":
        taskindex=1
    elif n=="Intrachain density correlation" or n=="Bulk intrachain density correlation":
        taskindex=2
    elif n=="Interchain density correlation" or n=="Bulk interchain density correlation":
        taskindex=3
    elif n=="Local density":
        taskindex=4
    elif n=="Interchain pairwise correlation" or n=="Bulk interchain pairwise correlation":
        taskindex=5
    elif n=="Second Order":
        taskindex=6
    elif n=="Entanglement Entropy":
        taskindex=7
    elif n=="Local density product":
        taskindex=8
    elif n=="Local density B":
        taskindex=9
    elif n=="Bulk superconducting corrleation":
        taskindex=10
    else:
        taskindex=11
    return taskindex

cfiles=[]

title=sys.argv[1]
for i in range(2,len(sys.argv)):
    cfiles.append(sys.argv[i])

filelist=os.listdir(os.getcwd())
if 'plots' not in filelist:
    os.mkdir('plots')

y=[]
pltlabels=[]
colors=[(170/255.0,0,45/255.0),(125/255.0,0,125.0/255.0)]

def readData(filename,data):
    lineIndex=0
    data.append([])
    with open(filename) as readData:
        for line in readData:
            lineIndex+=1
            buf=[x if x!='' else '0' for x in line.split('\t')]
            if len(buf)-1>i and lineIndex>5:
                data[len(data)-1].append(float(buf[i]))

for filename in filelist:
    if filename in cfiles:
            print filename
            with open(filename) as readCaption:
                readCaption.readline()
                readCaption.readline()
                parsR=readCaption.readline()
                names=readCaption.readline()
            datanames=names.split('\t')
            pars=parsR.split('\t')
            n=len(datanames)
            L=pars[0]
            np=pars[1]
            parity=pars[2]
            data=[]
            for i in range(0,n-1):
                readData(filename,data)
            y.append(data)
            pltlabels.append('J='+pars[3][0:4]+' g='+pars[4][0:4])
            #pltlabels.append('$\\alpha=$'+parity+'  E='+pars[6])
            #pltlabels.append('N='+pars[1])
            #ptitle='J='+pars[3]+'\t g='+pars[4]

markers=['o','<','>','v','s','^','*']
cols=['b',(200/255.0,0,45/255.0)]

fs=32
ls=28
for j in range(0,n-1):
    if datanames[j][0:4]!='Bulk':
        plt.figure(figsize=(12,10))
        plt.tick_params(labelsize=ls)
        print datanames[j]
        for i in range(0,len(y)):
            x=range(0,len(y[i][j]))
            if tasknum(datanames[j])!=7 and tasknum(datanames[j])!=4:
                cplot,=plt.semilogy(x,y[i][j],markers[i],ms=11.5)
            else:
                cplot,=plt.plot(x,y[i][j],markers[i],ms=11.5)
            cplot.set_label(pltlabels[i])
            #plt.xlabel('Distance $|i-j|$',fontsize=fs)
            #plt.xlabel('Subsystem size $S$',fontsize=fs)
        plt.xlabel('Site $i$',fontsize=fs)
        plt.ylabel(labellist[tasknum(datanames[j])],fontsize=fs)
        plt.xlim(xmin=0,xmax=100)
        plt.ylim(ymin=-.05,ymax=1.05)
        #plt.xticks([10,20,30])
        #ax=plt.gca()
        #ax.set_xticklabels(['10','20','30'])
        plt.legend(loc=2,numpoints=1,fontsize=ls)
        #plt.title(ptitle,fontsize=24)
        tnames=datanames[j].replace('.','_')
        if datanames[j]=='Local density':# or datanames[j]=='Interchain pairwise correlation':
            plt.savefig('thesis_plots/'+title+tnames.replace(' ','_')+'.pdf',bbox_inches='tight')
        plt.close()
            
            
