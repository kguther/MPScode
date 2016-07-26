import numpy as np
import os
import sys
import matplotlib.pyplot as plt

filelist=os.listdir(os.getcwd())

taskname=sys.argv[1]
J=sys.argv[2]
g=sys.argv[3]

files=[]
dE=[]
variances=[]
stdErr=[]

lenghts=[35,38,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92,95,98,101,104,107]
valid=[]

def getAverage(L):
    deltaE=[]
    energies=[]
    result=[]
    errBuff=[]
    err=[]
    files=[]
    for filename in filelist:
        if filename[0:len(taskname)]==taskname and filename[(len(filename)-6):(len(filename)-4)]!='ES':
            refName=filename.partition('_L_')
            if refName[2][:len(str(L))]==str(L):
                files.append(filename)
                with open(filename) as readEnergy:
                    readEnergy.readline()
                    readEnergy.readline()
                    eList=readEnergy.readline().split('\t')
                    if eList[3]==J and eList[4]==g:
                        valid.append(True)
                    else:
                        valid.append(False)
                    energies.append(float(eList[6]))
                    errBuff.append(float(eList[7]))
    for i in range(0,len(files)):
        refName=files[i]
        partedRefName=refName.partition('_p_')
        for j in range(0,len(files)):
            compName=files[j]
            partedCompName=compName.partition('_p_')
            if partedRefName[0]==partedCompName[0] and partedRefName[2][1:]==partedCompName[2][2:] and valid[i]:
                deltaE.append((energies[i]-energies[j])**2)
                err.append((errBuff[i]+errBuff[j])/2)

    averaged=sum(deltaE)/max(len(deltaE),1)
    errAveragedBuf=[]
    for val in deltaE:
        errAveragedBuf.append((val-averaged)**2)
    errAveraged=np.sqrt(sum(errAveragedBuf))/max(len(errAveragedBuf),1)
    print averaged
    print errAveraged
    result.append(averaged)
    result.append(errAveraged)
    return result

for l in lenghts:
    print l
    buf=getAverage(l)
    variances.append(buf[1])
    dE.append(buf[0])

plt.errorbar(lenghts,dE,yerr=variances,fmt='o')
plt.xlabel('Chain length')
plt.ylabel('Mean quadratic energy splitting')
plt.title('Energy gap, J=g=-2')
plt.xlim([34,109])
#plt.ylim([0.003,0.018])
#plt.savefig('thesis_plots/cdw_energy_scaling.pdf')
plt.show()
