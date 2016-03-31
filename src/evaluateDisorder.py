import numpy as np
import os
import sys

filelist=os.listdir(os.getcwd())

taskname=sys.argv[1]

energies=[]
files=[]
deltaE=[]

for filename in filelist:
    if filename[0:len(taskname)]==taskname and filename[(len(filename)-6):(len(filename)-4)]!='ES':
        files.append(filename)
        with open(filename) as readEnergy:
            readEnergy.readline()
            readEnergy.readline()
            eList=readEnergy.readline()
        energies.append(eList[6])
        files.append(filename)

for i in range(0,len(files)):
    refName=files[i]
    partedRefName=refName.partition('_p_')
    for j in range(0,len(files)):
        compName=files[j]
        partedCompName=compName.partition('_p_')
        if partedRefName[0]==partedCompName[0] and partedRefName[2][1:]==partedCompName[2][1:] and partedRefName[2][0]!=partedCompName[2][0]:
            deltaE.append((energies[i]-energies[j])**2)

print deltaE
averaged=sum(deltaE)/max(len(deltaE),1)
print averaged
