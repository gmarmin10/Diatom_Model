'''
Created on Jul 9, 2018

@author: Keisuke
'''

from pylab import *
from Claquin2002Nlimited import *
from Claquin2002Plimited import *
from Claquin2002Llimited import *
from FigSet import *
from matplotlib import *
from matplotlib.pyplot import figure, plot, xlabel, ylabel, xlim, ylim, show, xticks, yticks,legend,stackplot, margins, title
from sf_si import *
import matplotlib.patches as mpat
from numpy import *

#rcParams.update({'text.usetex':True})   #to call real latex
#rcParams.update({'text.latex.preamble':['\\usepackage[greek,english]{babel}']})
#rcParams.update({'font.serif': 'Times New Roman'})
# rc('text', usetex=True)
# rc('font', family='serif')
# rc('font', serif='Times New Roman')
# rcParams.update({'text.latex.preamble':['\\usepackage[greek,english]{babel}']})

rcParams.update({'mathtext.default': 'regular' })

d = claquin2002N()          #import Claquin data
dD = d.D                    #assigning columns to a single array/variable
dNC = d.N/d.C
dSiC = d.Si/d.C
dSiN = d.Si/d.N
dSi = d.Si


D =arange(0.01,1+0.01,0.001)     #growth rate


Qc = 9.408 #(pmol cell-1)
NCmin =0.03 #(molN molC)
Anc = 0.1
NC = NCmin + Anc*D

VsiCell=0.14 #(pmol Si cell-1 d-1)
Vsi = 0.016#(pmol Si mol C-1 s-1)
SiC = Vsi/D
SiN = SiC/NC

#parameters for thalassiosira from passsache
Vmax=VsiCell #(pmol Si cell-1 d-1)
Km=0.023  #(uM)
Vsi=VsiCell/Qc
Si = VsiCell/D 
Si_con = arange(0,1+1e-10,0.001) #(uM) Si concentration

V0=(Vmax*(Si_con))/(Km+Si_con)

 
Xlabel = '$\mu$ (d$^{-1}$)'

###plotting###
#data colors=teal
#model curves=wine-#882255
#allocation: biosynthesis=,C storage=, N storage=, P storage=, photo=,

figure(1)
plot(dD,dNC,'o',color='#882255',markersize=11,zorder=1)
plot(D,NC,color='teal',zorder=-1)
xlabel(Xlabel)
xlim(0,1)
xticks(fontsize=18)
yticks(fontsize=18)
ylim(0,0.125)
ylabel('N:C (mol mol$^{-1}$)')

modelleg=plot(D,NC,'-',color='teal', label="Model",zorder=-1)
dataleg=plot(dD,dNC,'o',color='#882255',markersize=11,label='Claquin et al., 2002')
legend(loc='lower right',fontsize='x-small',frameon=False)

sf('NC-Nlim_K',500)

figure(2)
plot(dD,dSiC,'o',color='#882255',markersize=11,zorder=1)
plot(D,SiC,color='teal',zorder=-1)
xlim(0,1)
ylim(0,0.2)
xlabel(Xlabel)
ylabel('Si:C (mol mol$^{-1}$)')
xticks(fontsize=18)
yticks(fontsize=18)
sf('SiC-Nlim_K',500)

figure(3)
plot(dD,dSiN,'o',color='#882255',markersize=11,zorder=1)
plot(D,SiN,color='teal',zorder=-1)
xlim(0,1)
ylim(0,5)
ylabel('Si:N (mol mol$^{-1}$)')
xlabel(Xlabel)
xticks(fontsize=18)
yticks(fontsize=18)
sf('SiN-Nlim_K',500)

#figure(10)
#plot(dD,d.Si,'o')

# figure (4)
# #Names = ["Oth","Pho","Bio"]
# Colors = ['#9BC2E6','#FFD966',"#A9D08E"]
# stackplot(Dd,Cconst_protein_plot,Cphoto_plot,Cbiosynth_plot,colors=Colors)
# legend(loc=4)
# ylim(0,70)

figure(4)
plot(dD,d.Si,'o', color='#882255',markersize=11,zorder=1)
plot(D,Si,color='teal',zorder=-1)
ylim(0,2)
xlim(0,1)
ylabel('Si ($\mu$mol cell$^{-1}$)')
xlabel(Xlabel)
xticks(fontsize=18)
yticks(fontsize=18)

sf('Si_Nlim_K',500)

##michaelis menten dynamics uptake of Si
figure (5)
plot(Si_con,V0, '-', color='teal')    # , markersize=11
xlim(0,0.5)
xlabel('Si ($\mu$M)')
ylabel('Uptake (pmol Si cell$^{-1}$d$^{-1}$)')
#ylim(0,0.08)
xticks(fontsize=18)
yticks(fontsize=18)
sf('Si_uptake_K',500)


###for nitrogen allocation figures:
Numbertoarray=ones(size(D))
NCessential_plot=NCmin*Numbertoarray
NCvariable=Anc*D

Ntot = NCessential_plot + NCvariable #+ NC_sto

figure(6)
Names = ["Protein min","Protein"]
Colors = ['#33BBEE','#CC6677']
stackplot([],[],[],colors=Colors[::-1],labels=Names[::-1])
stackplot(D,NCessential_plot,NCvariable,colors=Colors)
eleg=mpat.Patch(color='#4B0082', label="Essential",alpha=0.75)
pleg=mpat.Patch(color='#CC6677', label="Synthetic Protein",alpha=0.75)
legend(handles=[eleg,pleg],loc='upper center',bbox_to_anchor=(0.45,-0.35), ncol=2,fontsize='x-small',frameon=False)
xlabel('$\mu$ (d$^{-1}$)')                       #copied from 73
ylabel('(mol N mol C$^{-1}$)')  #copied from 73
title("N Allocation")
xlim(0,1)
#yticks(arange(0,125,step=25))
#ylim(top=100+1e-20)

sf('NC allocation_K',500)

#for C allocation
Cessential = 0.416*Numbertoarray*100
YcnProtein = 4.49
C_protein_min = NCmin*Numbertoarray*YcnProtein*100
C_protein_bio = Anc*D*YcnProtein*100
C_sto=100 - C_protein_min - C_protein_bio - Cessential
C_sto[C_sto<0] = 0
#print (Cpro_plot)


figure(7)
Names = ["Essential","Protein min","Protein bio","Storage"]
Colors = ['#4B0082','#33BBEE','#CC6677','#DDCC77']
stackplot([],[],[],[],[],colors=Colors[::-1],labels=Names[::-1])
stackplot(D,Cessential,C_protein_min,C_protein_bio,C_sto,colors=Colors)
bleg=mpat.Patch(color='#DDCC77',label='Storage',alpha=0.75)
pmin=mpat.Patch(color='#33BBEE',label='Essential Protein',alpha=0.75)
legend(handles=[bleg,pmin],loc='upper center',bbox_to_anchor=(0.30,-0.35), ncol=2,fontsize='x-small',frameon=False)
#legend(loc='upper left')
xlabel('$\mu$ (d$^{-1}$)')                       #copied from 73
ylabel('($\%$ of total C)')  #copied from 73
title('C allocation')
xlim(0,1)
yticks(arange(0,125,step=25))
ylim(0,top=100+1e-20)

sf('C allocation_K',500)

figure(8)
Names = ["Protein min","Protein","Storage"]
Colors = ['#9966FF','#CC6677','#FFC000']
stackplot([],[],[],colors=Colors[::-1],labels=Names[::-1])
stackplot(D,NCessential_plot/Ntot*100,NCvariable/Ntot*100,colors=Colors)
eleg=mpat.Patch(color='#4B0082', label="Essential",alpha=0.75)
pleg=mpat.Patch(color='#CC6677', label="Synthetic Protein",alpha=0.75)
#legend(handles=[eleg,pleg],loc='upper center',bbox_to_anchor=(0.65,-0.35), ncol=2,fontsize='x-small',frameon=False)
xlabel('$\mu$ (d$^{-1}$)')                       #copied from 73
ylabel('N allocation ($\%$)')  #copied from 73
xlim(0,1)



print(nanmean(d.C))
print(d.C)


show()



##### R2 calculations #####
D =arange(0.01,1+0.01,0.2)     #growth rate
dD = dD[2:]                   #assigning columns to a single array/variable
dNC = dNC[2:]
dSiC =dSiC[2:]
dSiN = dSiN[2:]
dSi = dSi[2:]

Qc = 9.408 #(pmol cell-1)
NCmin =0.03 #(molN molC)
Anc = 0.1
NC = NCmin + Anc*D

VsiCell=0.14 #(mol Si cell-1 d-1)
Vsi = 0.016#(mol Si mol C-1 s-1)
SiC = Vsi/D
SiN = SiC/NC

#parameters for thalassiosira from passsache
Vmax=0.073 #(h-1)
Km=1.4  #(uM)
Vsi=VsiCell/Qc
Si = VsiCell/D 

zipobject=zip(dSi,Si)

diff1=[]
for dSi_i,Si_i in zipobject:
    diff=dSi_i-Si_i
    diff1.append(diff)

sq1=[]
for i in diff1:
    sq=i**2
    sq1.append(sq)

sum_sq=sum(sq1)

n=len(dSi)
meanSi=sum(dSi)/n
nmeanSi=n*meanSi

diff2=[]
for i in dSi:
    diff= i-nmeanSi
    diff2.append(diff)
    
sumofdiff=sum(diff2)
Sq=sumofdiff**2
R2=1-(sum_sq/Sq)


print("r-squared_Si_cell=", R2)

zipobject=zip(dNC,NC)

diff1=[]
for dNC_i,NC_i in zipobject:
    diff=dNC_i-NC_i
    diff1.append(diff)

sq1=[]
for i in diff1:
    sq=i**2
    sq1.append(sq)

sum_sq=sum(sq1)

n=len(dNC)
meanNC=sum(dNC)/n
nmeanNC=n*meanNC

diff2=[]
for i in dNC:
    diff= i-nmeanNC
    diff2.append(diff)
    
sumofdiff=sum(diff2)
Sq=sumofdiff**2
R3=1-(sum_sq/Sq)
print("r-squared_N:C=", R3)


zipobject=zip(dSiN,SiN)

diff1=[]
for dSiN_i,SiN_i in zipobject:
    diff=dSiN_i-SiN_i
    diff1.append(diff)

sq1=[]
for i in diff1:
    sq=i**2
    sq1.append(sq)

sum_sq=sum(sq1)

n=len(dSiN)
meanSiN=sum(dSiN)/n
nmeanSiN=n*meanSiN

diff2=[]
for i in dSiN:
    diff= i-nmeanSiN
    diff2.append(diff)
    
sumofdiff=sum(diff2)
Sq=sumofdiff**2
R4=1-(sum_sq/Sq)
print("r-squared_Si:N=", R4)

zipobject=zip(dSiC,SiC)

diff1=[]
for dSiC_i,SiC_i in zipobject:
    diff=dSiC_i-SiC_i
    diff1.append(diff)

sq1=[]
for i in diff1:
    sq=i**2
    sq1.append(sq)

sum_sq=sum(sq1)

n=len(dSiC)
meanSiC=sum(SiC)/n
nmeanSiC=n*meanSiC

diff2=[]
for i in dSiC:
    diff= i-nmeanSiC
    diff2.append(diff)
    
sumofdiff=sum(diff2)
Sq=sumofdiff**2
R5=1-(sum_sq/Sq)
print("r-squared_Si:N=", R5)

