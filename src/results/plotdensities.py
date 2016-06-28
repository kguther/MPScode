import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
import scipy as sy
import sys

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


def makeName(prefix):
    newName=prefix
    for filename in sys.argv[1:]:
        addname=filename.split('/')
        newName+='_'+addname[1][:len(addname[1])-4]
    newName+='.pdf'
    return newName

N=[]
datalegend=[]
fulldat=[]
for filename in sys.argv[1:]:
        with open(filename) as rd:
            parRead=rd.readline()
            getJPre=(parRead.split('J='))[1]
            getJ=(getJPre.split(' g=')[0])
            getg=(getJPre.split(' g=')[1])
            getParPre=(parRead.split(' J='))[0]
            getPar=(getParPre.split('p='))[1]
            getNPre=(parRead.split(' p='))[0]
            getN=(getNPre.split('N='))[1]
        data=np.loadtxt(filename,skiprows=1)
        fulldat.append(data)
        datalegend.append("J="+getJ+" g="+getg+" $\\alpha$="+getPar)
        #datalegend.append("N="+getN)
        N.append(float(getN))

scope=plt.figure
i=0
markers=['o','<']
detail=[]
for dataset in fulldat:
    plt.plot(dataset[0],dataset[1],markers[0],ms=9)
    plt.plot(dataset[0],dataset[2],markers[1],ms=9)
    detail.append(datalegend[i]+" $, \\langle n^a_i \\rangle$")
    detail.append(datalegend[i]+" $, \\langle n^b_i \\rangle$")
    i+=1
plt.legend(detail,loc=8)
plt.xlabel('$i$')
plt.ylabel('$\\langle n_i \\rangle$')
#uses the last input name as new filename
plt.savefig('plots/'+filename[26:(len(filename)-4)]+'_comp_densities.pdf')

sym=plt.figure()
i=0
markers=['o','<','>','*']
for dataset in fulldat:
    plt.plot(dataset[0],dataset[1]+dataset[2],markers[i],ms=9)
    i+=1
plt.legend(datalegend,loc=9)
plt.xlabel('$i$')
plt.ylabel('$\\langle n_i^a \\rangle + \\langle n_i^b \\rangle $')
prefix='densities/symmetric_combination'
symname=makeName(prefix)
print symname
plt.savefig(symname)

symfit=plt.figure()
for dataset in fulldat:
    p0=sy.array([0.1,0.7])
    fpars,acc=so.curve_fit(symfunc,dataset[0],dataset[1]+dataset[2],p0)
    plt.plot(dataset[0],dataset[1]+dataset[2],'o')
    plt.plot(dataset[0],symfunc(dataset[0],fpars[0],fpars[1]))
plt.xlabel('$i$')
plt.ylabel('$\\langle n_i^a \\rangle + \\langle n_i^b \\rangle $')
datalegendFit=[]
for entry in datalegend:
    datalegendFit.append(entry)
    datalegendFit.append('a='+str(fpars[1])+' b='+str(fpars[0]))
plt.title('Fit with $a\cosh(b(i-\\frac{L}{2})) + \mathrm{const.}$')
plt.legend(datalegendFit,loc=9)
prefix='densities/fit_to_symmetric'
symfitname=makeName(prefix)
#plt.savefig(symfitname)

asym=plt.figure()
for dataset in fulldat:
    plt.plot(dataset[0],dataset[1]-dataset[2],'o',ms=9)
plt.legend(datalegend,loc=9)
prefix='densities/asymmetric_combination'
asymname=makeName(prefix)
plt.xlabel('$i$')
plt.ylabel('$\\langle n_i^a \\rangle - \\langle n_i^b \\rangle $')
plt.savefig(asymname)

asymfit=plt.figure()
i=0
datalegendFit=[]
for dataset in fulldat:
    kf=np.pi/float(np.max(dataset[0]))*N[i]
    offset=np.mean(dataset[1]-dataset[2])
    scale=np.max(dataset[1]-dataset[2])-offset
    p0=sy.array([offset,scale,0])
    fpars,acc=so.curve_fit(asymfunc,dataset[0],dataset[1]-dataset[2],p0)
    plt.plot(dataset[0],dataset[1]-dataset[2],'o')
    plt.plot(dataset[0],asymfunc(dataset[0],fpars[0],fpars[1],fpars[2]))
    datalegendFit.append(datalegend[i])
    datalegendFit.append('a='+str(fpars[1])+' b='+str(fpars[0]))
    i+=1
plt.xlabel('$i$')
plt.ylabel('$\\langle n_i^a \\rangle - \\langle n_i^b \\rangle $')
plt.title('Fit with $a\cos(k_f x + c) + b$')
#plt.legend(datalegendFit,loc=9)
prefix='densities/fit_to_asymmetric'
asymfitname=makeName(prefix)
#plt.savefig(asymfitname)
