import numpy as np
import matplotlib.pyplot as plt
import os

filelist=os.listdir(os.getcwd())
if 'plots' not in filelist:
    os.mkdir('plots')

taskname="run_1"

writepd=False

labellist=['$\\left|\\langle a^\dagger_i a_0^{} \\rangle \\right|$','$\\left|\\langle a^\dagger_i b^{}_i a_0^{\dagger} b_0 \\rangle \\right|$','$\\left|\\langle n^{a}_i n_0^{a} \\rangle \\right|$','$\\left|\\langle n^{a}_i n_0^{b} \\rangle \\right|$','$\\left|\\langle n^{a}_i \\rangle \\right|$','$\\left|\\langle a^\dagger_i b^\dagger_i a_0^{} b_0^{} \\rangle \\right|$','$\\left|\\langle \\right|\\rangle$']

def tasknum(n):
    if n=="Intrachain correlation" or n=="Bulk correlation function":
        taskindex=0
    elif n=="Interchain hopping correlation" or n=="Bulk interchain hopping correlation":
        taskindex=1
    elif n=="Intrachain density correlation" or n=="Bulk intrachain density correlation":
        taskindex=2
    elif n=="Interchain density correlation" or n=="Bulk interchain density correlation":
        taskindex=3
    elif n=="Local density" or n=="Local density B":
        taskindex=4
    elif n=="Interchain pairwise correlation" or n=="Bulk interchain pairwise correlation":
        taskindex=5
    elif n=="Second Order":
        taskindex=6
    else:
        taskindex=7
    return taskindex

if writepd:
    with open('phasediagram'+'.txt','w') as pd:
        pd.write('J\tg\tdensity fluctuation\n')

for filename in filelist:
    if filename[0:len(taskname)]==taskname:
        with open(filename) as readCaption:
            readCaption.readline()
            readCaption.readline()
            pars=readCaption.readline()
            names=readCaption.readline()
        datanames=names.split('\t')
        parVals=pars.split('\t')
        data=np.loadtxt(filename,skiprows=5)
        n=len(data[0])
        data=data.transpose()
        L=len(data[0])
        for i in range(0,n):
            x=range(2,L+2)
            bCheck=datanames[i].split(' ')
            if bCheck[0]=='Bulk':
                x=range(0,L/2)
            tasklabel=labellist[tasknum(datanames[i])]
            if (tasknum(datanames[i])==4 and writepd):
                phase=min(data[i])/max(data[i])
                point=parVals[0]+'\t'+parVals[1]+'\t'+str(phase)+'\n'
                with open('phasediagram.txt','a') as pd:
                    pd.write(point)
            plt.figure()
            if(bCheck[0]=='Bulk'):
                plt.loglog(x,abs(data[i][x]),'o')
            else:
                if(tasknum(datanames[i])!=4):
                    plt.semilogy(x,abs(data[i]),'o')
                else:
                    plt.plot(x,data[i],'o')
            plt.xlabel('distance i')
            plt.ylabel(tasklabel)
            tname=datanames[i].replace('.','_')
            plt.title('J='+parVals[0]+' g='+parVals[1]+' E='+parVals[2]+' $\\Delta E$='+parVals[3])
            plt.savefig('plots/'+filename[0:len(filename)-4]+'_'+tname.replace(' ','_')+'.pdf')
            plt.close()
