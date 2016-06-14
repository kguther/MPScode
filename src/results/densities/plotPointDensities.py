import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
import scipy as sy
import sys
import os

def symfunc(x,a,b):
    return b*np.cosh(a*(x))

def asymfunc(x,a,b,c):
    return b*np.cos(kf*x+c)+a

def effAbs(x):
    ref=x[mid]
    redX=[]
    for y in x:
        redX.append(abs(ref-y))
    return redX

taskname=sys.argv[1]
L=sys.argv[2]
J=sys.argv[3]
g=sys.argv[4]

def makeName(prefix):
    newName=prefix
    newName+=taskname+'_J_'+J+'_g_'+g
    newName+='.pdf'
    return newName

filelist=os.listdir(os.getcwd())

N=[]
datalegendCollection=[[],[],[],[]]
fulldatCollection=[[],[],[],[]]
for filename in filelist:
    if filename[0:len(taskname)]==taskname:
        with open(filename) as rd:
            parRead=rd.readline()
            getJPre=(parRead.split('J='))[1]
            getJ=(getJPre.split(' g=')[0])
            getg=(getJPre.split(' g=')[1])
            getg=getg.strip()
            getParPre=(parRead.split(' J='))[0]
            getPar=(getParPre.split('p='))[1]
            getNPre=(parRead.split(' p='))[0]
            getN=(getNPre.split('N='))[1]
            getLPre=(parRead.split(' N='))[0]
            getL=(getLPre.split('L='))[1]
        if getJ==J and getg==g and getL==L:
            data=np.loadtxt(filename,skiprows=1)
            if getPar=='1':
                if float(getN)%2==0:
                    k=0
                else:
                    k=1
            else:
                if float(getN)%2==0:
                    k=2
                else:
                    k=3
            fulldatCollection[k].append(data)
            datalegendCollection[k].append("N="+getN+' $\\alpha$='+getPar)
            N.append(float(getN))
            compname=filename

scope=plt.figure
markers=['o','<']
print len(fulldatCollection[k])
subplotcoords=[221,222,223,224]
for j in [0,1,2,3]:
    detail=[]
    plt.subplot(subplotcoords[j])
    i=0
    for dataset in fulldatCollection[j]:
        plt.plot(dataset[0],dataset[1],markers[0],ms=9)
        #plt.plot(dataset[0],dataset[2],markers[1],ms=9)
        detail.append(datalegendCollection[j][i]+" $, \\langle n^a_i \\rangle$")
        #detail.append(datalegends[j][i]+" $, \\langle n^b_i \\rangle$")
        i+=1
        plt.legend(detail,loc=8, fontsize=10)
        plt.xlabel('$i$')
        plt.ylabel('$\\langle n_i \\rangle$')
        #uses the last input name as new filename
plt.savefig('../plots/'+compname[16:(len(compname)-4)]+'_comp_densities.pdf')


sym=plt.figure()
markers=['o','<','>','*','d','s']
for j in range(0,4):
    i=0
    plt.subplot(subplotcoords[j])
    for dataset in fulldatCollection[j]:
        plt.plot(dataset[0],dataset[1]+dataset[2],markers[i],ms=9)
        i+=1
        plt.legend(datalegendCollection[j],loc=8, fontsize=10)
        plt.xlabel('$i$')
        plt.ylabel('$\\langle n_i^a \\rangle + \\langle n_i^b \\rangle $')
        prefix='symmetric_combination'
        symname=makeName(prefix)
plt.savefig(symname)

symfit=plt.figure()
for j in range(0,4):
    plt.subplot(subplotcoords[j])
    for dataset in fulldatCollection[j]:
        p0=sy.array([0.1,0.7])
        fpars,acc=so.curve_fit(symfunc,dataset[0],dataset[1]+dataset[2],p0)
        plt.plot(dataset[0],dataset[1]+dataset[2],'o')
        plt.plot(dataset[0],symfunc(dataset[0],fpars[0],fpars[1]))
    plt.xlabel('$i$')
    plt.ylabel('$\\langle n_i^a \\rangle + \\langle n_i^b \\rangle $')
    datalegendFit=[]
    for entry in datalegendCollection[j]:
        datalegendFit.append(entry)
        datalegendFit.append('a='+str(fpars[1])+' b='+str(fpars[0]))
    plt.legend(datalegendFit,loc=8, fontsize=10)
plt.title('Fit with $a\cosh(b(i-\\frac{L}{2})) + \mathrm{const.}$')
prefix='fit_to_symmetric'
symfitname=makeName(prefix)
#plt.savefig(symfitname)

asym=plt.figure()
for j in range(0,4):
    plt.subplot(subplotcoords[j])
    for dataset in fulldatCollection[0]:
        plt.plot(dataset[0],dataset[1]-dataset[2],'o',ms=9)
    plt.legend(datalegendCollection[j],loc=8, fontsize=10)
    plt.xlabel('$i$')
    plt.ylabel('$\\langle n_i^a \\rangle - \\langle n_i^b \\rangle $')
prefix='asymmetric_combination'
asymname=makeName(prefix)
plt.savefig(asymname)

asymfit=plt.figure()
i=0
datalegendFit=[]
for dataset in fulldatCollection[0]:
    kf=np.pi/float(np.max(dataset[0]))*N[i]
    offset=np.mean(dataset[1]-dataset[2])
    scale=np.max(dataset[1]-dataset[2])-offset
    p0=sy.array([offset,scale,0])
    fpars,acc=so.curve_fit(asymfunc,dataset[0],dataset[1]-dataset[2],p0)
    plt.plot(dataset[0],dataset[1]-dataset[2],'o')
    plt.plot(dataset[0],asymfunc(dataset[0],fpars[0],fpars[1],fpars[2]))
    datalegendFit.append(datalegendCollection[0][i])
    datalegendFit.append('a='+str(fpars[1])+' b='+str(fpars[0]))
    i+=1
plt.xlabel('$i$')
plt.ylabel('$\\langle n_i^a \\rangle - \\langle n_i^b \\rangle $')
plt.title('Fit with $a\cos(k_f x + c) + b$')
#plt.legend(datalegendFit,loc=8, fontsize=10)
prefix='fit_to_asymmetric'
asymfitname=makeName(prefix)
#plt.savefig(asymfitname)
