'''
Created on Jul 9, 2018

@author: Keisuke
'''

from pylab import *
from Claquin2002Nlimited import *
from Claquin2002Plimited import *
from Claquin2002Llimited import *
from Claquin2002PlimitedV import *
from FigSet import *
from matplotlib import *
from matplotlib.pyplot import figure, plot, xlabel, ylabel, xlim, ylim, show, yticks, xticks, legend,stackplot, title
from sf_si2 import *
import matplotlib.patches as mpat

# rc('text', usetex=True)
# rc('font', family='serif')
# rc('font', serif='Times New Roman')
# rcParams.update({'text.latex.preamble':['\\usepackage[greek,english]{babel}']})

rcParams.update({'mathtext.default': 'regular' })

d = claquin2002P()
dD = d.D
dNC = d.N/d.C
dSiC = d.Si/d.C
dSiN = d.Si/d.N

D =arange(0.01,1+0.01,0.01)

V = -113.18*D + 578.86 #(m3 cell-1)        #y-int from table 1 in Claquin volume         linear interpolation
Qc = -16.552*D + 25.409 #(pmol cell-1)     #y-int from table 1 in Claquin carbon per cell    linear interpolation
VsiCell = 0.15   #(mol Si cell-1 d-1)
NCmin =0.03 #(molN molC)
Anc = 0.1
NCstoMax = 0.035
NC = NCmin + Anc*D + NCstoMax

Vsi = 0.016   #(mol Si mol C-1 d-1)
Vsi = 0.008
Vsi = VsiCell/Qc
SiC = Vsi/D
SiN = SiC/NC

Si = VsiCell/D  
Si_con = arange(0,1+1e-10,0.001) #(uM) Si concentration
Vmax=VsiCell #(mol Si cell-1 d-1)
Km=0.023  #(uM)
V0=(Vmax*(Si_con))/(Km+Si_con)

Xlabel = '$\mu$ (d$^{-1}$)'

###plotting###
#data colors=teal
#model curves=wine-#882255
#allocation: biosynthesis=,C storage=, N storage=, P storage=, photo=,

figure(1)
plot(dD,dNC,'o', color='#882255',markersize=11,zorder=1)
plot(D,NC,color='#999933',zorder=-1)
xlabel(Xlabel)
ylabel('N:C (mol mol$^{-1}$)')
xlim(0,1)
ylim(top=0.175)
#ylim(0,0.175)
xticks(fontsize=18)
yticks(fontsize=18)

modelleg=plot(D,NC,'-',color='#999933', label="Model",zorder=-1)
dataleg=plot(dD,dNC,'o',color='#882255',markersize=11,label='Claquin et al., 2002')
legend(loc='upper left',fontsize='x-small',frameon=False)

sf('NC-Plim_K',500)

figure(2)
plot(dD,dSiC,'o', color='#882255',markersize=11,zorder=1)
plot(D,SiC,color='#999933',zorder=-1)
xlim(0,1)
ylim(0,0.1)
xlabel(Xlabel)
ylabel('Si:C (mol mol$^{-1}$)')
xticks(fontsize=18)
yticks(fontsize=18)

sf('SiC-Plim_K',500)

figure(3)
plot(dD,dSiN,'o', color='#882255',markersize=11,zorder=1)
plot(D,SiN,color='#999933',zorder=-1)
ylim(0,1)
ylabel('Si:N (mol mol$^{-1}$)')
xlabel(Xlabel)
xlim(0,1)
xticks(fontsize=18)
yticks(fontsize=18)

sf('SiN-Plim_K',500)

figure(4)
plot(dD,d.Si,'o', color='#882255',markersize=11,zorder=1)
plot(D,Si,color='#999933',zorder=-1)
ylim(0,2)
xlim(0,1)
ylabel('Si (mol cell$^{-1}$)')
xlabel(Xlabel)
xticks(fontsize=18)
yticks(fontsize=18)

sf('Si-Plim_K',500)

figure (5)
plot(Si_con,V0, '-', color='#999933')    # , markersize=11
xlim(0,0.5)
xlabel('Si ($\mu$M)')
ylabel('Uptake (pmol Si cell$^{-1}$d$^{-1}$)')
#ylim(0,0.08)
xticks(fontsize=18)
yticks(fontsize=18)
sf('Si_uptake_K',500)

Numbertoarray=ones(size(D))
NCessential_plot=NCmin*Numbertoarray
NCvariable=Anc*D
NC_sto=NCstoMax*Numbertoarray


# yticks(arange(0,125,step=25))
# ylim(top=100+1e-20)



#for C allocation
# Cessential_plot=NCmin*Numbertoarray*Qc*100
# Cpro_plot=100*Anc*D*Qc*((100-Cessential_plot)/100)
# C_sto=100-Cessential_plot-Cpro_plot
# C_sto[C_sto<0] = 0
# print (Cpro_plot)

Cessential = 0.416*Numbertoarray*100
YcnProtein = 4.49
C_protein_min = NCmin*Numbertoarray*YcnProtein*100
C_protein_bio = Anc*D*YcnProtein*100
#C_N_sto = NC_sto*2*100
#print(C_N_sto)
C_sto = 100 - C_protein_min - C_protein_bio - Cessential
C_sto[C_sto<0] = 0
#C_N_sto = 100 - C_protein_min - C_protein_bio - Cessential - C_sto
#NC_sto = C_N_sto/2/100

Ntot = NCessential_plot + NCvariable + NC_sto

figure(6)
Names = ["Protein min","Storage","Protein"]
Colors = ['#33BBEE','#DDCC77','#CC6677']
stackplot([],[],[],colors=Colors[::-1],labels=Names[::-1])
stackplot(D,NCessential_plot,NC_sto,NCvariable,colors=Colors)
eleg=mpat.Patch(color='#4B0082', label="Essential",alpha=0.75)
pleg=mpat.Patch(color='#CC6677', label="Synthetic Protein",alpha=0.75)
#legend(handles=[eleg,pleg],loc='upper center',bbox_to_anchor=(0.65,-0.35), ncol=2,fontsize='x-small',frameon=False)
xlabel('$\mu$ (d$^{-1}$)')                       #copied from 73
ylabel('(mol N mol C$^{-1}$)')  #copied from 73
#title("N Allocation")
xlim(0,1)

sf('NC allocation_no_legend_K',500)

figure(7)
Names = ["Essential","Protein min","Protein","C Storage"]
Colors = ['#4B0082','#33BBEE','#CC6677','#DDCC77']
stackplot([],[],[],[],[],colors=Colors[::-1],labels=Names[::-1])
stackplot(D,Cessential,C_protein_min,C_protein_bio,C_sto,colors=Colors)
bleg=mpat.Patch(color='#DDCC77',label='Storage',alpha=0.75)
pmin=mpat.Patch(color='#33BBEE',label='Essential Protein',alpha=0.75)
#legend(handles=[bleg,pmin],loc='upper center',bbox_to_anchor=(0.30,-0.35), ncol=2,fontsize='x-small',frameon=False)
#legend(loc='upper left')
xlabel('$\mu$ (d$^{-1}$)')                       #copied from 73
ylabel('($\%$ of total C)')  #copied from 73
#title('C allocation')
xlim(0,1)
yticks(arange(0,125,step=25))
ylim(0,top=100+1e-20)

sf('C allocation_no-legend_K',500)

figure(8)
Names = ["Protein min","Storage","Protein"]
Colors = ['#9966FF','#FFC000','#CC6677']
stackplot([],[],[],colors=Colors[::-1],labels=Names[::-1])
stackplot(D,NCessential_plot/Ntot*100,NC_sto/Ntot*100,NCvariable/Ntot*100,colors=Colors)
eleg=mpat.Patch(color='#4B0082', label="Essential",alpha=0.75)
pleg=mpat.Patch(color='#CC6677', label="Protein",alpha=0.75)
#legend(handles=[eleg,pleg],loc='upper center',bbox_to_anchor=(0.65,-0.35), ncol=2,fontsize='x-small',frameon=False)
xlabel('$\mu$ (d$^{-1}$)')                       #copied from 73
ylabel('N allocation ($\%$)')  #copied from 73
xlim(0,1)


print(nanmean(d.C))
print(d.C)

show()


##### R2 calculations #####
D =arange(0.01,1+0.01,0.2)     #growth rate
dD = dD[2:8]                   #assigning columns to a single array/variable
dNC = dNC[2:8]
dSiC =dSiC[2:8]
dSiN = dSiN[2:8]
dSi = d.Si[2:8]

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
meanNC=sum(NC)/n
nmeanNC=n*meanNC

diff2=[]
for i in dSi:
    diff= i-nmeanNC
    diff2.append(diff)
    
sumofdiff=sum(diff2)
Sq=sumofdiff**2
R2=1-(sum_sq/Sq)


print("r-squared N:C=", R2)