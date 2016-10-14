import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import scipy.optimize as so
import scipy as sy

def f(x,a,c,d):
    return c+d/(x**a)

def fluctuation(array):
    rho=np.mean(array)
    buf=array
    for i in range(0,len(buf)-1):
        buf[i]=(array[i]-rho)/(array[i+1]-rho)
    return mean(buf)

filelist=os.listdir(os.getcwd())
if 'plots' not in filelist:
    os.mkdir('plots')

taskname=sys.argv[1]

readDensities=True
writeDeg=False
defaultLegs=['Ground state', '1st excited State']

savePlot=False
firstfile=True

labellist=['$\\left|\\langle a^\dagger_i a_0^{} \\rangle \\right|$','$\\left|\\langle a^\dagger_i b^{}_i b_0^{\dagger} a_0 \\rangle \\right|$','$\\left|\\langle n^{a}_i n_0^{a} \\rangle \\right|$','$\\left|\\langle n^{a}_i n_0^{b} \\rangle \\right|$','$\\left|\\langle n^{a}_i \\rangle \\right|$','$\\left|\\langle a^\dagger_i b^\dagger_i a_0^{} b_0^{} \\rangle \\right|$','$\\left|\\langle \\right|\\rangle$','S','$\\left|\\langle n^{a}_i n^{b}_i \\rangle\\right|$','$\\left|\\langle n^{b}_i\\rangle\\right|$','$\\langle a_i^{\dagger} a_{i+1}^{\dagger} a_j a_{j+1} \\rangle$','$\\langle b_i^{\dagger} b_{i+1}^{\dagger} a_j a_j \\rangle$','$\\langle (n^a_i - n^b_i) (n^a_0 -n^b_0) \\rangle$','other']

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
    elif n=="Total magnetization correlation":
        taskindex=12
    elif n=="Bulk interchain superconducting corrleation":
        taskindex=11
    else:
        taskindex=13
    return taskindex

for filename in filelist:
    if filename[0:len(taskname)]==taskname:
        if filename[(len(filename)-6):(len(filename)-4)]!='ES' and filename[(len(filename)-12):(len(filename)-4)]!='_state_2':
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
            if firstfile:
                if writeDeg:
                    with open('SB_degs_L_'+L+'_N_'+np+'.txt','w') as dp:
                        dp.write('L\tE_0\tE_1\tdE_0\tdE_1\n')
            densityOutput=[]
            densityFileName='densities/local_densities_'+filename
            if readDensities:
                with open(densityFileName,'w') as dfn:
                    dfn.write('Local Densites for L='+L+' N='+np+' p='+parity+' J='+pars[3]+' g='+pars[4]+'\n')
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
                
                if readDensities:
                    if datanames[i]=="Local density":
                        densityOutput.append(range(0,len(data)))
                        densityOutput.append(data)
                    if datanames[i]=="Local density B":
                        densityOutput.append(data)
                    if i==n-2:
                        with open(densityFileName,'a') as dfn:
                            for k in range(0,len(densityOutput)):
                                for point in densityOutput[k]:
                                    dfn.write(str(point)+'\t')
                                dfn.write('\n')

                if i==0 and writeDeg:
                    refName=filename.partition('_run_')
                    refNameLast=refName[2].partition('_L_')
                    run=refNameLast[0]
                    degdata=run+'\t'+pars[6]+'\t'+enList[6]+'\t'+pars[7].strip()+'\t'+enList[7].strip()+'\n'
                    with open('SB_degs_L_'+L+'_N_'+np+'.txt','a') as dp:
                        dp.write(degdata)
                x=range(0,len(data))
                xB=range(0,len(dataB))
                tasklabel=labellist[tasknum(datanames[i])]
                bCheck=datanames[i].split(' ')
                plt.figure()
                if bCheck[0]=='Bulk':
                    if tasknum(datanames[i])==10:
                        plt.plot(x,data,'o')
                        if excitedState:
                            plt.plot(xB,dataB,'o')
                            plt.legend(defaultLegs)
                    else:
                        if tasknum(datanames[i]) in [0,1,5]:
                            plt.semilogy(x,map(abs,data),'o')
                        else:
                            plt.loglog(x,map(abs,data),'o')
                            if excitedState:
                                plt.loglog(xB,map(abs,dataB),'o')
                                plt.legend(defaultLegs)
                else:
                    if (tasknum(datanames[i])!=4 and tasknum(datanames[i])!=7 and tasknum(datanames[i])!=9 and tasknum(datanames[i])!=8):
                        plt.semilogy(x,map(abs,data),'o')
                        if excitedState:
                            plt.semilogy(xB,map(abs,dataB),'o')
                            plt.legend(defaultLegs)
                    else:
                        plt.plot(x,data,'o')
                        if excitedState:
                            plt.plot(xB,dataB,'o')
                            plt.legend(defaultLegs)
                plt.xlabel('distance i')
                plt.ylabel(tasklabel)
                tname=datanames[i].replace('.','_')
                plt.title('J='+pars[3]+' g='+pars[4]+' W= '+pars[5]+' E='+pars[6]+' $(\\Delta E)^2$='+pars[7].strip())
                if savePlot:
                    plt.savefig('plots/'+filename[0:len(filename)-4]+'_'+tname.replace(' ','_')+'.pdf')
                if datanames[i]=='Local density':
                    plt.show()
                plt.close()
                firstfile=False

                                
