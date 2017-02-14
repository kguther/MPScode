#import matplotlib 
#from matplotlib import rc

#rc('text',usetex=True)
#rc('text.latex',preamble='\usepackage{color} \definecolor{jred}{RGB}{170,0,45} \definecolor{gviolet}{RGB}{125,0,125}')


import numpy as np
import matplotlib.pyplot as plt
import os
import sys
labellist=['$\\left|\\langle a^\dagger_i a_j^{} \\rangle \\right|$','$\\left|\\langle a^\dagger_i b^{}_i b_1^{\dagger} a_1 \\rangle \\right|$','$\\left|\\langle n^{a}_i n_0^{a} \\rangle \\right|$','$\\left|\\langle n^{a}_i n_0^{b} \\rangle \\right|$','$\\left|\\langle n^{a}_i \\rangle \\right|$','$\\left|\\langle a^\dagger_i b^\dagger_i a_1^{} b_1^{} \\rangle \\right|$','$\\left|\\langle \\right|\\rangle$','S','$\\left|\\langle n^{a}_i n^{b}_i \\rangle\\right|$','$\\left|\\langle n^{b}_i\\rangle\\right|$','$\\langle a_i^{\dagger} a_{i+1}^{\dagger} a_j a_{j+1} \\rangle$','other']

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
    elif n=="Bulk intrachain superconducting corrleation":
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
#colors=[(174/255.0,28.0/255.0,97.0/255.0),(21.0/255.0,84.0/255.0,193.0/255.0)]
colors=[(52.0/255.0,137.0/255.0,197.0/255.0), (137.0/255.0,199.0/255.0,58.0/255.0)]
fs=32
ls=28

def readData(filename,data):
    lineIndex=0
    data.append([])
    with open(filename) as readData:
        for line in readData:
            lineIndex+=1
            buf=[x if x!='' else '0' for x in line.split('\t')]
            if len(buf)-1>i and lineIndex>5:
                data[len(data)-1].append(float(buf[i]))

def plotArray(y,n):
    markers=['o','s','>','v','s','^','*']
    dfs=12
    plt.figure(figsize=(12,10))
    plt.tick_params(labelsize=ls)
    #targetName=['Intrachain correlation']
    #targetName=['Interchain hopping correlation', 'Interchain pairwise correlation']
    targetName=['Local density']
    k=0
    for i in range(0,len(y)):
        for j in range(0,n[i]-1): 
            if datanames[i][j] in targetName:
                x=range(0,len(y[i][j]))
                if tasknum(datanames[i][j])in [1,5]:
                    cplot,=plt.semilogy(x,y[i][j],markers[k],color=colors[k],ms=dfs)
                else:
                    if tasknum(datanames[i][j]) in [0,6,10]:
                        cplot,=plt.semilogy(x,y[i][j],markers[k],ms=dfs,color=colors[k])
                    else:
                        cplot,=plt.plot(x,y[i][j],markers[k],ms=dfs,color=colors[k])
                #cplot.set_label(labellist[tasknum(datanames[i][j])])
                cplot.set_label(pltlabels[i])
                k+=1
                #plt.xlabel('$|i-j|$',fontsize=fs)
                plt.xlabel('Site $i$',fontsize=fs)
    if(len(targetName)==1):
        plt.ylabel(labellist[tasknum(targetName[0])],fontsize=fs)
    if(len(targetName)!=1 or len(y)!=1):
        plt.legend(loc=3,numpoints=1,fontsize=ls)
    tnames=targetName[0].replace('.','_')
    plt.savefig('../../draft/plots/'+title+tnames.replace(' ','_')+'.pdf',bbox_inches='tight')
    plt.show()
    plt.close()

def plotDeg(y,n):
    marker='o'
    cols=['w',colors[0]]
    marksize=[18,12]
    dfs=11.5
    rim=[3,0]

    plt.figure(figsize=(12,10))
    plt.tick_params(labelsize=ls)

    targetName='Intrachain correlation'

    for i in range(0,len(y)):
        for j in range(0,n[i]-1): 
            if datanames[i][j] in targetName:
                x=range(0,len(y[i][j]))
                cplot,=plt.semilogy(x,y[i][j],marker,color=cols[i],ms=marksize[i],mew=rim[i])
                cplot.set_label(pltlabels[i])
    plt.xlabel('Site $i$',fontsize=fs)
    #plt.xlabel('$|i-j|$', fontsize=fs)
    plt.ylabel(labellist[tasknum(targetName)],fontsize=fs)
    plt.legend(loc=9,numpoints=1,fontsize=ls)
    plt.xlim(xmax=83)
    tnames=targetName.replace('.','_')
    plt.savefig('../../draft/plots/'+title+tnames.replace(' ','_')+'.pdf',bbox_inches='tight')
    plt.show()
    plt.close()

def plotMulti(y,n):
    ca='Interchain pairwise correlation'
    cb='Interchain hopping correlation'
    #ca='Total density correlation'
    #cb='Total magnetization correlation'
    #ca='Local density'
    #cb='Local density B'

    fig,ax1=plt.subplots(figsize=(12,12))
    ax1.tick_params(labelsize=ls)
    if(ca!='Local density'):
        ax1.set_xlabel('$|i-j|$',fontsize=fs)
    else:
        ax1.set_xlabel('$Site i$', fontsize=fs)
    ax2=ax1.twinx()
    ax2.tick_params(labelsize=ls)
    for i in range(0,len(y)):
        for j in range(0,n[i]-1):
            if datanames[i][j]==ca or datanames[i][j]==cb:
                x=range(0,len(y[i][j]))
                if datanames[0][j]==ca:
                    if(ca=='Interchain pairwise correlation'):
                        cplot,=ax1.loglog(x,y[i][j],markers[0],ms=marksize[1],color=colors[0])
                    else:
                        cplot,=ax1.plot(x,y[i][j],markers[0],ms=marksize[1],color=colors[0])
                else:
                    if(ca=='Interchain pairwise correlation'):
                        ax2.loglog(x,y[i][j],markers[1],ms=marksize[1],color=colors[1])
                    else:
                        ax2.plot(x,y[i][j],markers[1],ms=marksize[1],color=colors[1])

    if(ca=='Interchain pairwise correlation'):
        ax1.set_ylabel('$\\left|\\langle a^\dagger_i b^\dagger_i a_j^{} b_j^{} \\rangle \\right|$',fontsize=fs,color=colors[0])
        ax2.set_ylabel('$\\left|\\langle a^\dagger_i b^{}_i b_j^{\dagger} a_j \\rangle \\right|$',fontsize=fs,color=colors[1])
    if(ca=='Total density correlation'):
        ax1.set_ylabel('$\\langle (n^a_i + n^b_i) (n^a_j +n^b_j) \\rangle$',fontsize=fs,color=colors[0])
        ax2.set_ylabel('$\\langle (n^a_i - n^b_i) (n^a_j -n^b_j) \\rangle$',fontsize=fs,color=colors[1])
    if(ca=='Local density'):
        ax1.set_ylabel('$\\langle n^a_i \\rangle$',fontsize=fs,color=colors[0])
        ax2.set_ylabel('$\\langle n^b_i \\rangle$',fontsize=fs,color=colors[1])
    for yt in ax1.get_yticklabels():
        yt.set_color(colors[0])
    for yt in ax2.get_yticklabels():
        yt.set_color(colors[1])
    plt.savefig('../../draft/plots'+title+'multiPlot'+'.pdf',bbox_inches='tight')
    plt.show()
plt.close()

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
            pltlabels.append('$\\alpha=$'+parity+'  E='+pars[6][:9])
            #pltlabels.append('E='+pars[6][:9])
            #pltlabels.append('N='+pars[1])
            #ptitle='J='+pars[3]+'\t g='+pars[4]

#plotArray(y,n)
plotDeg(y,n)
