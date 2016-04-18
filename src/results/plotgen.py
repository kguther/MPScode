import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import scipy.optimize as so
import scipy as sy

def f(x,a,c,d):
    return c+d/(x**a)

filelist=os.listdir(os.getcwd())
if 'plots' not in filelist:
    os.mkdir('plots')

taskname=sys.argv[1]

writeK=False
writepd=False
newpd=False
defaultLegs=['Ground state', '1st excited State']

labellist=['$\\left|\\langle a^\dagger_i a_0^{} \\rangle \\right|$','$\\left|\\langle a^\dagger_i b^{}_i a_0^{\dagger} b_0 \\rangle \\right|$','$\\left|\\langle n^{a}_i n_0^{a} \\rangle \\right|$','$\\left|\\langle n^{a}_i n_0^{b} \\rangle \\right|$','$\\left|\\langle n^{a}_i \\rangle \\right|$','$\\left|\\langle a^\dagger_i b^\dagger_i a_0^{} b_0^{} \\rangle \\right|$','$\\left|\\langle \\right|\\rangle$','S','$\\left|\\langle n^{a}_i n^{b}_i \\rangle\\right|$','$\\left|\\langle n^{b}_i\\rangle\\right|$','$\\langle a_i^{\dagger} a_{i+1}^{\dagger} a_0 a_1 \\rangle$','other']

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

for filename in filelist:
    if filename[0:len(taskname)]==taskname:
        if filename[(len(filename)-6):(len(filename)-4)]!='ES' and filename[(len(filename)-6):(len(filename)-4)]!='_2':
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
            if filename==filelist[0]:
                if writeK:
                    with open('decay_pars_L_'+L+'_N_'+np+'_p_'+parity+'.txt','w') as kp:
                        kp.write('J\tg\t')
                if writepd and newpd:
                    with open('phasediagram_L_'+L+'_N_'+np+'_p_'+parity+'.txt','w') as pd:
                        pd.write('J\tg\tdensity fluctuation\n')
            for i in range(0,n-1):
                data=[]
                dataB=[]
                lineIndex=0
                with open(filename) as readData:
                    for line in readData:
                        lineIndex+=1
                        buf=[x if x!='' else '0' for x in line.split('\t')]
                        if len(buf)-1>i and lineIndex>5:
                            data.append(float(buf[i]))
                filenameB=filename[:(len(filename)-4)]+'_state_2.txt'
                print filenameB
                lineIndex=0
                if filenameB in filelist:
                    with open(filenameB) as readEnergy:
                        readEnergy.readline()
                        readEnergy.readline()
                        enListR=readEnergy.readline()
                    enList=enListR.split('\t')
                    defaultLegs=[pars[6],enList[6]]
                    with open(filenameB) as readData:
                        for line in readData:
                            lineIndex+=1
                            buf=[x if x!='' else '0' for x in line.split('\t')]
                            if len(buf)-1>i and lineIndex>5:
                                dataB.append(float(buf[i]))
                if len(dataB)>0:
                    excitedState=True
                else:
                    excitedState=False
                x=range(0,len(data))
                tasklabel=labellist[tasknum(datanames[i])]
                bCheck=datanames[i].split(' ')
                if bCheck[0]=='Bulk':
                    if writeK:
                        with open('decay_pars_L_'+L+'_N_'+np+'_p_'+parity+'.txt','a') as kp:
                            kp.write('exponent of '+datanames[i]+'\t')
                        if tasknum(datanames[i])==5:
                            xeff=range(int(L)/10,len(data))
                            p0=sy.array([1,1,1])
                            fpars, acc=so.curve_fit(f,xeff,data[int(L)/10:len(data)],p0)
                if (tasknum(datanames[i])==4 and writepd):
                    phase=min(data)/max(data)
                    point=pars[3]+'\t'+pars[4]+'\t'+str(phase)+'\n'
                    with open('phasediagram_L_'+L+'_N_'+np+'_p_'+parity+'.txt','a') as pd:
                        pd.write(point)
                plt.figure()
                if bCheck[0]=='Bulk':
                    if tasknum(datanames[i])==10:
                        plt.plot(x,data,'o')
                        if excitedState:
                            plt.plot(x,dataB,'o')
                            plt.legend(defaultLegs)
                    else:
                        plt.loglog(x,map(abs,data),'o')
                        if excitedState:
                            plt.loglog(x,map(abs,dataB),'o')
                            plt.legend(defaultLegs)
                    if writeK and tasknum(datanames[i])==5:
                        plt.loglog(xeff,f(xeff,fpars[0],fpars[1],fpars[2]),'o')
                else:
                    if (tasknum(datanames[i])!=4 and tasknum(datanames[i])!=7 and tasknum(datanames[i])!=9 and tasknum(datanames[i])!=8):
                        plt.semilogy(x,map(abs,data),'o')
                        if excitedState:
                            plt.semilogy(x,map(abs,dataB),'o')
                            plt.legend(defaultLegs)
                    else:
                        plt.plot(x,data,'o')
                        if excitedState:
                            plt.plot(x,dataB,'o')
                            plt.legend(defaultLegs)
                plt.xlabel('distance i')
                plt.ylabel(tasklabel)
                tname=datanames[i].replace('.','_')
                plt.title('J='+pars[3]+' g='+pars[4]+' W= '+pars[5]+' E='+pars[6]+' $(\\Delta E)^2$='+pars[7].strip())
                plt.savefig('plots/'+filename[0:len(filename)-4]+'_'+tname.replace(' ','_')+'.pdf')
                plt.close()

            
