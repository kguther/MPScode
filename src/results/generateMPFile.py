import os
import sys

taskname=sys.argv[1]
target=sys.argv[2]
filelist=os.listdir(os.getcwd())

JValues=[x*1/3.0 for x in range(-12,7)]
gValues=[x*1/3.0 for x in range(-15,12)]
alphas=[-1,1]

req=[[x,y,z] for x in JValues for y in gValues for z in alphas]
av=[]
miss=[]

for filename in filelist:
    if filename[0:len(taskname)]==taskname:
        if filename[(len(filename)-6):(len(filename)-4)]!='ES' and filename[(len(filename)-11):(len(filename)-4)]!='state_2':
            with open(filename) as readPars:
                readPars.readline()
                readPars.readline()
                parsR=readPars.readline()
                namesR=readPars.readline()
            pars=parsR.split('\t')
            if len(namesR)!=0:
                av.append([pars[3],pars[4],pars[2]])
            else:
                miss.append(pars[3]+'\t'+pars[4]+'\t'+pars[2]+'\n')

for x in req:
    missing=1
    for y in av:
        if abs(x[0]-float(y[0]))<1e-2 and abs(x[1]-float(y[1]))<1e-2 and x[2]==int(y[2]):
            missing=0
    if missing==1:
        miss.append(str(x[0])+'\t'+str(x[1])+'\t'+str(x[2])+'\n')

with open(target,'w') as mp:
    for point in miss:
        mp.write(point)
