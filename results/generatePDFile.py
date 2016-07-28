import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import scipy.optimize as so
import scipy as sy

def fluctuation(array):
    rho=np.mean(array)
    buf=np.empty(len(array))
    for i in range(0,len(buf)-1):
        buf[i]=abs((array[i]-rho))
    return np.mean(buf)

def relativeFluctuation(arrayA,arrayB):
    buf=[min([arrayA[i],arrayB[i]]) for i in range(0,len(arrayA))]
    return np.mean(buf)

def geometricMean(array):
    buf=1.0
    for i in range(1,len(array)-1):
        buf*=array[i]
    return buf**(1/float(len(array)))

def readParameters(filename,pars):
    with open(filename) as readCaption:
        readCaption.readline()
        readCaption.readline()
        parsR=readCaption.readline()
        names=readCaption.readline()
    pars.append(parsR.split('\t'))
    return names.split('\t')

def readData(filename,data):
    lineIndex=0
    data.append([])
    with open(filename) as readData:
        for line in readData:
            lineIndex+=1
            buf=[x if x!='' else '0' for x in line.split('\t')]
            if len(buf)-1>i and lineIndex>5:
                data[len(data)-1].append(float(buf[i]))

taskname=sys.argv[1]
pdName='pd_'+taskname.rstrip('_')+'.txt'

with open(pdName,'w') as pd:
    pd.write('rho\talpha\tJ\tg\tgs degeneracy\tdensity fluctuation\tgreens function revival\tgreens function minimum\tsdw-par\tcdw-par\tentanglement entropy parameter\tenergy variance\taccuracy\n')
points=[]

counter=0

filelist=os.listdir(os.getcwd())

for filename in filelist:
    if filename[0:len(taskname)]==taskname:
        if filename[(len(filename)-6):(len(filename)-4)]!='ES' and filename[(len(filename)-11):(len(filename)-4)]!='state_2':
            #print filename
            pars=[]
            cdwParA=[]
            cdwParB=[]
            phase=[]
            revival=[]
            cmin=[]
            datanames=readParameters(filename,pars)
            n=len(datanames)
            if(pars[0][2]=='1'):
                conjugateParity='-1'
            else:
                conjugateParity='1'
            prefix=filename.split('_p_')[0]
            postfix=filename.split('_W_')[1]
            conjugateFilename=prefix+'_p_'+conjugateParity+'_W_'+postfix
            #print conjugateFilename
            if conjugateFilename in filelist:
                readParameters(conjugateFilename,pars)
            for i in range(0,n-1):
                data=[]
                readData(filename,data)
                if conjugateFilename in filelist:
                    readData(conjugateFilename,data)
                #copy the relevant data into containers
                if (datanames[i]=="Local density"):
                    densA=data
                if (datanames[i]=="Local density B"):
                    densB=data
                if (datanames[i]=="Intrachain correlation"):
                    cor=data
                if (datanames[i]=="Entanglement Entropy"):
                    S=data
            if len(pars)==len(densA):
                sdw=[]
                entanglement=[]
                for j in range(0,len(pars)):
                    #get the characteristics for each phase from the containers
                    if len(densA[j])!=0:
                        phase.append(min(densA[j])/max(densA[j]))
                        cdwParA.append(fluctuation(densA[j]))
                    else:
                        phase.append(0.0)
                        cdwParA.append(0.0)
                    if len(densB[j])!=0:
                        cdwParB.append(fluctuation(densB[j]))
                    else:
                        cdwParB.append(0.0)
                    if len(densA[j])!=0 and len(densB[j])!=0:
                        sdw.append(relativeFluctuation(densA[j],densB[j]))
                    else:
                        sdw.append(0.0)
                    if len(cor[j])!=0:
                        if abs(cor[j][0])>1e-12:
                            revival.append(cor[j][len(cor[j])-2]/cor[j][0])
                            if abs(cor[j][len(cor[j])-2])>1e-12:
                                if np.min(cor[j][::2])>np.min(cor[j][1:len(cor[j])-3:2]) or True:
                                    cmin.append(np.min(cor[j][::2])/cor[j][len(cor[j])-2])
                                else:
                                    cmin.append(0.5)
                            else:
                                cmin.append(1.0)
                        else:
                            revival.append(0.0)
                            cmin.append(1.0)
                    else:
                        revival.append(0.0)
                        cmin.append(0.0)
                    entanglement.append(np.mean(S[j]))      
                if conjugateFilename in filelist:
                    degeneracy=abs(float(pars[0][6])-float(pars[1][6]))     
                else:
                    degeneracy=10.0
                                
                cPoint=[]
                if len(pars)>1:
                    #only for the W=0.1 grid since there, the alpha=-1 data has delta E=0 due to a bug
                    if float(pars[0][7])>float(pars[1][7]):
                        j=0
                    else:
                        j=1
                else:
                    j=0
                cdw=(cdwParA[j]+cdwParB[j])/2
                acc=np.sqrt(float(pars[j][7]))/max([abs(float(pars[j][6])),1e-7])
                cPoint.append(str(float(pars[j][1])/float(pars[j][0]))+'\t'+pars[j][2]+'\t'+pars[j][3]+'\t'+pars[j][4]+'\t'+str(degeneracy)+'\t'+str(phase[j])+'\t'+str(revival[j])+'\t'+str(cmin[j])+'\t'+str(sdw[j])+'\t'+str(cdw)+'\t'+str(entanglement[j])+'\t'+pars[j][7]+'\t'+str(acc)+'\n')
                with open(pdName,'a') as pd:
                    for wPoint in cPoint:
                        pd.write(wPoint)

            
            if conjugateFilename in filelist:
                filelist.remove(conjugateFilename)
