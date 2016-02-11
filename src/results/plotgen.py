import numpy as np
import matplotlib.pyplot as plt
import os

filelist=os.listdir(os.getcwd())
if 'plots' not in filelist:
    os.mkdir('plots')

taskname="run_1"

writepd=False

labellist=['$\\left|\\langle a^\dagger_i a_0^{} \\rangle \\right|$','$\\left|\\langle a^\dagger_i b^{}_i a_0^{\dagger} b_0 \\rangle \\right|$','$\\left|\\langle n^{a}_i n_0^{a} \\rangle \\right|$','$\\left|\\langle n^{a}_i n_0^{b} \\rangle \\right|$','$\\left|\\langle n^{a}_i \\rangle \\right|$','$\\left|\\langle a^\dagger_i b^\dagger_i a_0^{} b_0^{} \\rangle \\right|$']

def tasknum(n):
  if n=="Intrachain correlation":
    taskindex=0
  elif n=="Interchain hopping correlation":
    taskindex=1
  elif n=="Intrachain density correlation":
    taskindex=2
  elif n=="Interchain density correlation":
    taskindex=3
  elif n=="Local density":
    taskindex=4
  elif n=="Interchain pairwise correlation":
    taskindex=5
  else:
    taskindex=6
  return taskindex

if writepd:
    with open('phasediagram'+'.txt','w') as pd:
        pd.write('J\tg\tdensity fluctuation\n')

for filename in filelist:
    if filename[0:len(taskname)]==taskname:
        with open(filename) as readCaption:
            readCaption.readline()
            pars=readCaption.readline()
            names=readCaption.readline()
        datanames=names.split('\t')
        parVals=pars.split('\t')
        data=np.loadtxt(filename,skiprows=3)
        n=len(data[0])
        data=data.transpose()
        L=len(data[0])
        x=range(2,L+2)
        for i in range(0,n):
            tasklabel=labellist[tasknum(datanames[i])]
            if (tasknum(datanames[i])==4 && writepd):
                phase=min(data[i])/max(data[i])
                point=parVals[0]+'\t'+parVals[1]+'\t'+str(phase)+'\n'
                with open('phasediagram.txt','a') as pd:
                    pd.write(point)
            plt.figure()
            plt.plot(x,data[i],'o')
            plt.xlabel('site i')
            plt.ylabel(tasklabel)
            plt.title('J='+parVals[0]+' g='+parVals[1]+' E='+parVals[2]+' $\\Delta E$='+parVals[3])
            plt.savefig('plots/'+filename[0:len(filename)-4]+'_'+datanames[i].replace(' ','_')+'.pdf')
            plt.close()
