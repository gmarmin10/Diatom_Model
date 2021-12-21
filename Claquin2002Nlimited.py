'''
Created on Jul 9, 2018

@author: Keisuke
'''
from numpy import *
from matplotlib import *
from pylab import *

class Claquin2002N:
    def __init__(self,D,C,N,Si,V):
        self.D = D
        self.C = C
        self.N = N
        self.Si = Si
        self.V = V

def claquin2002N():
    d = genfromtxt('Claquin2002Nlimited.csv',delimiter=',')
    
    Data = Claquin2002N(d[:,0],d[:,1],d[:,2],d[:,3],d[:,7])
    
    return Data