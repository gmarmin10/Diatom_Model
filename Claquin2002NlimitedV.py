'''
Created on Jul 9, 2018

@author: Keisuke
'''
from matplotlib import *
from numpy import *
from pylab import *

class Claquin2002Nv:
    def __init__(self,D,C,N,Si,V):
        self.D = D
        self.C = C
        self.N = N
        self.Si = Si
        self.V = V

def claquin2002Nv():
    d = genfromtxt('Claquin2002Nlimited.csv',delimiter=',')
    
    Data = Claquin2002Nv(d[:,0],d[:,8],d[:,9],d[:,10],d[:,7])
    
    return Data