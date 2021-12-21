'''
Created on Jul 9, 2018

@author: Keisuke
'''
from numpy import *
from matplotlib import *
from pylab import *

class Claquin2002Pv:
    def __init__(self,D,C,N,Si,V):
        self.D = D
        self.C = C
        self.N = N
        self.Si = Si
        self.V = V

def claquin2002Pv():
    d = genfromtxt('Claquin2002Plimited.csv',delimiter=',')
    
    Data = Claquin2002Pv(d[:,0],d[:,8],d[:,9],d[:,10],d[:,7])
    
    return Data