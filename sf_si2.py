'''
Created on Nov 1, 2021
Specifically for the Diatom model N limitation figures

'''
from pylab import * 
from matplotlib.pyplot import savefig

def sf(figName,Dpi):
    First_part="C:/Users/19046/Documents/General_Research/Diatom_Model/figures/Si/Plim/"
    #Second_part=savefolder+"\\"
    Figure_name=str(figName)
    Last_part=".tif"
    savefig(First_part+Figure_name+Last_part,dpi=Dpi)
