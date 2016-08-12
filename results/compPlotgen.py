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
datanames=[]
pltlabels=[]
n=[]
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
            datanames.append(names.split('\t'))
            pars=parsR.split('\t')
            n.append(len(datanames[len(datanames)-1]))
            L=pars[0]
            np=pars[1]
            parity=pars[2]
            data=[]
            for i in range(0,n[len(n)-1]-1):
                readData(filename,data)
            y.append(data)
            #pltlabels.append('J='+pars[3][0:4]+' g='+pars[4][0:4])
            #pltlabels.append('$\\alpha=$'+parity+'  E='+pars[6][:7])
            pltlabels.append('E='+pars[6][:7])
            #pltlabels.append('N='+pars[1])
            #ptitle='J='+pars[3]+'\t g='+pars[4]

markers=['o','<','>','v','s','^','*']
cols=['b','g']
markersize=11.5

fs=32
ls=28
plt.figure(figsize=(13,10))
plt.tick_params(labelsize=ls)

targetName='Intrachain correlation'

for i in range(0,len(y)):
    for j in range(0,n[i]-1): 
        if datanames[i][j]==targetName:
            x=range(0,len(y[i][j]))
            if tasknum(datanames[i][j])!=7 and tasknum(datanames[i][j])!=4:
                cplot,=plt.semilogy(x,y[i][j],markers[i],ms=markersize,color=colors[i])
            else:
                cplot,=plt.plot(x,y[i][j],markers[i],ms=markersize,color=colors[i])
            cplot.set_label(pltlabels[i])
#plt.xlabel('Distance $|i-j|$',fontsize=
#plt.xlabel('Subsystem size $S$',fontsize=fs)
plt.xlabel('Site $i$',fontsize=fs)
plt.ylabel(labellist[tasknum(targetName)],fontsize=fs)
#plt.ylabel('Entanglement Entropy',fontsize=fs)
#plt.xlim(xmin=0,xmax=100)
plt.ylim(ymin=-.05,ymax=0.65)
#plt.xticks([10,20,30])
#ax=plt.gca()
#ax.set_xticklabels(['10','20','30'])
plt.legend(loc=9,numpoints=1,fontsize=ls)
#plt.title(ptitle,fontsize=24)
tnames=targetName.replace('.','_')
plt.savefig('thesis_plots/'+title+tnames.replace(' ','_')+'.pdf',bbox_inches='tight')
#plt.show()
plt.close()


fig,ax1=plt.subplots(figsize=(12,10))
ax1.tick_params(labelsize=ls)
ax1.set_xlabel('$|i-j|$',fontsize=fs)
ax2=ax1.twinx()
ax2.tick_params(labelsize=ls)
for i in range(0,len(y)):
    for j in range(0,n[i]-1):
        if datanames[i][j]=='Interchain pairwise correlation' or datanames[i][j]=='Interchain hopping correlation':
            x=range(0,len(y[i][j]))
            if datanames[0][j]=='Interchain pairwise correlation':
                cplot,=ax1.semilogy(x,y[i][j],markers[i],ms=markersize,color=colors[0])
                cplot.set_label(pltlabels[i])
            else:
                ax2.semilogy(x,y[i][j],markers[i],ms=markersize,color=colors[1])

ax1.set_ylabel('$\\left|\\langle a^\dagger_i b^\dagger_i a_j^{} b_j^{} \\rangle \\right|$',fontsize=fs,color=colors[0])
ax2.set_ylabel('$\\left|\\langle a^\dagger_i b^{}_i b_j^{\dagger} a_j \\rangle \\right|$',fontsize=fs,color=colors[1])
#ax1.legend(loc=3,numpoints=1,fontsize=ls)
#ax1.set_xlim(xmax=100)
#ax2.set_xlim(xmax=100)
for yt in ax1.get_yticklabels():
    yt.set_color(colors[0])
for yt in ax2.get_yticklabels():
    yt.set_color(colors[1])
#plt.savefig('thesis_plots/'+title+'pairCorrelations'+'.pdf',bbox_inches='tight')
#plt.show()
plt.close()
            
