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

rcParams.update({'mathtext.default': 'regular' })
#retrieving Claquin data
d = claquin2002N()          #import Claquin data
dD = d.D                    #assigning columns to a single array/variable
dNC = d.N/d.C
dSiC = d.Si/d.C
dSiN = d.Si/d.N
dSi = d.Si

#Model equations
mu =arange(0.01,1+0.01,0.001)     #growth rate
Qc = 9.408 #(pmol cell-1)
Qp_max = Qc/65    #mol P/cell        Quigg, King,finkel(?) all have ratios below 100        0.00979 from Inomura 2020 #NEED A NEW VALUE FOR THIS!!
Pro_essential =0.03 #(molN molC)
PCmin = 0.00072025     #perry et al (molP/molC)
Pro_synth = 0.1       #mol N mol C-1 d
NC = Pro_essential + Pro_synth*mu
x_RNA = 4.23E-3*4.49    #4.23e-3 (mol P d mol C-1) 4.49 (mol C molN-1)# supplemental material of Inomura 2020
RNA = NC*x_RNA*mu                           
P_sto = (Qp_max/Qc)-PCmin-RNA        
PC = PCmin + RNA +P_sto


VsiCell=0.14 #(pmol Si cell-1 d-1)
Vsi=VsiCell/Qc
Si = VsiCell/mu 
SiC = Vsi/mu
SiN = SiC/NC
SiP = SiC/PC

#figure plotting
Xlabel = '$\mu$ (d$^{-1}$)'

figure(1)
plot(dD,dNC,'o',color='#882255',markersize=11,zorder=1)
plot(mu,NC,color='teal',zorder=-1)
xlabel(Xlabel)
xlim(0,1)
xticks(fontsize=18)
yticks(arange(0,0.180,0.025),fontsize=18)
ylim(0,top=0.175)
ylabel('N:C (mol mol$^{-1}$)')

sf('NC-Nlim_withp',500)

figure(2)
plot(dD,dSiC,'o',color='#882255',markersize=11,zorder=1)
plot(mu,SiC,color='teal',zorder=-1)
xlim(0,1)
ylim(0,0.2)
xlabel(Xlabel)
ylabel('Si:C (mol mol$^{-1}$)')
xticks(fontsize=18)
yticks(fontsize=18)
sf('SiC-Nlim_withp',500)

figure(3)
plot(dD,dSiN,'o',color='#882255',markersize=11,zorder=1)
plot(mu,SiN,color='teal',zorder=-1)
xlim(0,1)
ylim(0,5)
ylabel('Si:N (mol mol$^{-1}$)')
xlabel(Xlabel)
xticks(fontsize=18)
yticks(fontsize=18)
sf('SiN-Nlim__withp',500)

figure(10)
plot(mu,SiP,color='teal',zorder=-1)
xlim(0,1)
ylim(0,50)
ylabel('Si:P (mol mol$^{-1}$)')
xlabel(Xlabel)
xticks(fontsize=18)
yticks(fontsize=18)
sf('SiP-Nlim__withp',500)


figure(4)
plot(dD,d.Si,'o', color='#882255',markersize=11,zorder=1)
plot(mu,Si,color='teal',zorder=-1)
ylim(0,5)
xlim(0,1)
ylabel('Si (mol cell$^{-1}$)')
xlabel(Xlabel)
xticks(fontsize=18)
yticks(fontsize=18)
modelleg=plot(mu,Si,'-',color='teal', label="Model",zorder=-1)
dataleg=plot(dD,d.Si,'o',color='#882255',markersize=11,label='Claquin et al., 2002')
legend(loc='upper right',fontsize='x-small',frameon=False)

sf('Si_Nlim__withp',500)


###for nitrogen allocation figures:
Numbertoarray=ones(size(mu))
NCessential_plot=Pro_essential*Numbertoarray
NCvariable=Pro_synth*mu

Ntot = NCessential_plot + NCvariable #+ NC_sto

figure(6)
Names = ["Protein min","Protein"]
Colors = ['#33BBEE','#CC6677']
stackplot([],[],[],colors=Colors[::-1],labels=Names[::-1])
stackplot(mu,NCessential_plot,NCvariable,colors=Colors)
eleg=mpat.Patch(color='#4B0082', label="Essential",alpha=0.75)
pleg=mpat.Patch(color='#CC6677', label="Synthesis",alpha=0.75)
xlabel('$\mu$ (d$^{-1}$)')                      
ylabel('(mol N mol C$^{-1}$)')  
yticks(fontsize=18)
xticks(arange(0,1.25,step=0.25),fontsize=18)
xlim(0,1)
ylim(0,0.20)

sf('NC allocation__withp',500)

#for C allocation
Cessential = 0.416*Numbertoarray*100
YcnProtein = 4.49
C_protein_min = Pro_essential*Numbertoarray*YcnProtein*100
C_protein_bio = Pro_synth*mu*YcnProtein*100
C_sto=100 - C_protein_min - C_protein_bio - Cessential
C_sto[C_sto<0] = 0



figure(7)
Names = ["Essential","Protein min","Protein bio","Storage"]
Colors = ['#4B0082','#33BBEE','#CC6677','#DDCC77']
stackplot([],[],[],[],[],colors=Colors[::-1],labels=Names[::-1])
stackplot(mu,Cessential,C_protein_min,C_protein_bio,C_sto,colors=Colors)
bleg=mpat.Patch(color='#DDCC77',label='Storage',alpha=0.75)
pmin=mpat.Patch(color='#33BBEE',label='Essential Protein',alpha=0.75)
xlabel('$\mu$ (d$^{-1}$)')                       
ylabel('($\%$ of total C)')  
xlim(0,1)
xticks(fontsize=18)
yticks(arange(0,125,step=25),fontsize=18)
ylim(0,top=100+1e-20)

sf('C allocation__withp',500)

# setup for p allocation figure
Numbertoarray=ones(size(mu))

Pessential = PCmin*Numbertoarray
PCvariable = RNA
P_sto = P_sto*Numbertoarray
P_sto[P_sto<0] = 0 

figure(9)
Names = ["Essential","Storage","RNA"]
Colors = ['#4B0082','#DDCC77','#CC6677']
stackplot([],[],[],[],colors=Colors[::-1],labels=Names[::-1])
stackplot(mu,Pessential,P_sto,PCvariable,colors=Colors)
bleg=mpat.Patch(color='#DDCC77',label='Storage',alpha=0.75)
pmin=mpat.Patch(color='#33BBEE',label='Essential Protein',alpha=0.75)
xlabel('$\mu$ (d$^{-1}$)')                       
ylabel('(mol P mol C$^{-1}$)')  
xlim(0,1)
ylim(0,top=0.020)
yticks(fontsize=18)
xticks(arange(0,1.25,step=0.25),fontsize=18)
sf('P allocation__withp',500)

NP=NC/PC

figure(11)
plot(mu,NP,color='teal',zorder=-1)
ylabel('N:P (mol mol$^{-1}$)')
xlabel(Xlabel)
xlim(0,1)
ylim(0,top=10)
xticks(fontsize=18)
yticks(arange(0,12,step=2),fontsize=18)

sf('NP_withp', 500)


print(nanmean(d.C))
print(d.C)


show()

