'''
Created on Jul 9, 2018

@author: Keisuke
'''

from pylab import *

class Claquin2002L:
    def __init__(self,D,C,N,Si,V):
        self.D = D
        self.C = C
        self.N = N
        self.Si = Si
        self.V = V

def claquin2002L():
    d = genfromtxt('Claquin2002LightLimited.csv',delimiter=',')
    
    Data = Claquin2002L(d[:,0],d[:,1],d[:,2],d[:,3],d[:,7])
    
    return Data