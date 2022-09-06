'''

@author: garmin, Keisuke
'''

from pylab import *
from Claquin2002Nlimited import *
from Claquin2002Plimited import *
from Claquin2002Llimited import *
from Claquin2002PlimitedV import *
from FigSet import *
from matplotlib import *
from matplotlib.pyplot import figure, plot, xlabel, ylabel, xlim, ylim, show, yticks, xticks, legend,stackplot, title
import matplotlib.patches as mpat

rcParams.update({'mathtext.default': 'regular' })

#retrieving data points from Claquin et al 2002
d = claquin2002P()
dD = d.D
dNC = d.N/d.C
dSiC = d.Si/d.C
dSiN = d.Si/d.N

# Model equations
mu =arange(0.01,1+0.01,0.01)

V = -113.18*mu + 578.86 #(m3 cell-1)        #y-int from table 1 in Claquin volume         linear interpolation
Qc = -16.552*mu + 25.409 #(pmol cell-1)     #y-int from table 1 in Claquin carbon per cell    linear interpolation
VsiCell = 0.14   #(mol Si cell-1 d-1)
Pro_essential =0.03 #(molN molC)
Pro_synth = 0.1
NCsto = 0.035
NC = Pro_essential + Pro_synth*mu + NCsto
PCmin = 0.00072025
x_RNA = 4.23E-3*4.49                    #mol P *d* mol N-1: Values from inomura 2020 supplemental
RNA = NC*x_RNA*mu     
PC = PCmin + RNA 

Si_max=0.95
Vsi = VsiCell/Qc
Si = VsiCell/mu  
Si[Si>Si_max]=Si_max
SiC = Si/Qc
SiN = SiC/NC
SiP = SiC/PC
#Figure making
Xlabel = '$\mu$ (d$^{-1}$)'

figure(1)
plot(dD,dNC,'o', color='#882255',markersize=11,zorder=1)
plot(mu,NC,color='#999933',zorder=-1)
xlabel(Xlabel)
ylabel('N:C (mol mol$^{-1}$)')
xlim(0,1)
xticks(fontsize=18)
yticks(arange(0,0.180,0.025),fontsize=18)
ylim(0,top=0.175)


figure(2)
plot(dD,dSiC,'o', color='#882255',markersize=11,zorder=1)
plot(mu,SiC,color='#999933',zorder=-1)
xlim(0,1)
ylim(0,0.2)
xlabel(Xlabel)
ylabel('Si:C (mol mol$^{-1}$)')
xticks(fontsize=18)
yticks(fontsize=18)


figure(3)
plot(dD,dSiN,'o', color='#882255',markersize=11,zorder=1)
plot(mu,SiN,color='#999933',zorder=-1)
ylim(0,5)
ylabel('Si:N (mol mol$^{-1}$)')
xlabel(Xlabel)
xlim(0,1)
ylim(0,4)
xticks(fontsize=18)
yticks(fontsize=18)


figure(10)
plot(mu,SiP,color='#999933',zorder=-1)
xlim(0,1)
ylim(0,top=60)
ylabel('Si:P (mol mol$^{-1}$)')
xlabel(Xlabel)
xticks(fontsize=18)
yticks(fontsize=18)



figure(4)
plot(dD,d.Si,'o', color='#882255',markersize=11,zorder=1)
plot(mu,Si,color='#999933',zorder=-1)
ylim(0,1.5)
xlim(0,1)
ylabel('Si (pmol cell$^{-1}$)')
xlabel(Xlabel)
xticks(fontsize=18)
yticks(fontsize=18)
modelleg=plot(mu,Si,'-',color='#999933', label="Model",zorder=-1)
dataleg=plot(dD,d.Si,'o',color='#882255',markersize=11,label='Claquin et al., 2002')
legend(loc='upper right',fontsize='x-small',frameon=False)


Numbertoarray=ones(size(mu))
NCessential_plot=Pro_essential*Numbertoarray
NCvariable=Pro_synth*mu
NC_sto=NCsto*Numbertoarray

Cessential = 0.416*Numbertoarray*100
YcnProtein = 4.49
C_protein_min = Pro_essential*Numbertoarray*YcnProtein*100
C_protein_bio = Pro_synth*mu*YcnProtein*100

C_sto = 100 - C_protein_min - C_protein_bio - Cessential
C_sto[C_sto<0] = 0

Ntot = NCessential_plot + NCvariable + NC_sto

figure(6)
Names = ["Protein min","Storage","Protein"]
Colors = ['#33BBEE','#DDCC77','#CC6677']
stackplot([],[],[],colors=Colors[::-1],labels=Names[::-1])
stackplot(mu,NCessential_plot,NC_sto,NCvariable,colors=Colors)
eleg=mpat.Patch(color='#4B0082', label="Essential",alpha=0.75)
pleg=mpat.Patch(color='#CC6677', label="Synthetic Protein",alpha=0.75)
xlabel('$\mu$ (d$^{-1}$)')                       
ylabel('(mol N mol C$^{-1}$)')  
xlim(0,1)
ylim(top=0.20)
xticks(fontsize=18)
yticks(fontsize=18)

figure(7)
Names = ["Essential","Protein min","Protein","C Storage"]
Colors = ['#4B0082','#33BBEE','#CC6677','#DDCC77']
stackplot([],[],[],[],[],colors=Colors[::-1],labels=Names[::-1])
stackplot(mu,Cessential,C_protein_min,C_protein_bio,C_sto,colors=Colors)
bleg=mpat.Patch(color='#DDCC77',label='Storage',alpha=0.75)
pmin=mpat.Patch(color='#33BBEE',label='Essential Protein',alpha=0.75)
xlabel('$\mu$ (d$^{-1}$)')                       
ylabel('($\%$ of total C)')  
xlim(0,1)
yticks(arange(0,125,step=25),fontsize=18)
xticks(fontsize=18)
ylim(0,top=100+1e-20)


# setup for p allocation figure
Numbertoarray=ones(size(mu))

Pessential = PCmin*Numbertoarray
PCvariable = RNA 

figure(8)
Names = ["Essential","RNA"]
Colors = ['#4B0082','#CC6677']
stackplot([],[],[],colors=Colors[::-1],labels=Names[::-1])
stackplot(mu,Pessential,PCvariable,colors=Colors)
bleg=mpat.Patch(color='#DDCC77',label='Storage',alpha=0.75)
pmin=mpat.Patch(color='#CC6677',label='Essential Protein',alpha=0.75)
xlabel('$\mu$ (d$^{-1}$)')                       
ylabel('(mol P mol C$^{-1}$)')  
xlim(0,1)
ylim(top=0.020)

yticks(fontsize=18)
xticks(arange(0,1.25,step=0.25),fontsize=18)


NP=NC/PC

figure(9)
plot(mu,NP,color='#999933',zorder=-1)
ylabel('N:P (mol mol$^{-1}$)')
xlabel(Xlabel)
xlim(0,1)
ylim(0,top=100)
xticks(fontsize=18)
yticks(arange(0,120,step=20),fontsize=18)


print(nanmean(d.C))
print(d.C)

show()
