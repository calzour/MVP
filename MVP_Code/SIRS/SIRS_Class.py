import copy
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from scipy import ndimage
import warnings
import random

class SIRS():
    def __init__(self,p1,p2,p3,lattice_size=50):
        '''
        p1 is probability of S=>I,
        p2 is probability of I=>R
        p3 is probability of R=>S
        '''
        self.lattice_size=lattice_size
        self.p1=p1
        self.p2=p2
        self.p3=p3
        self.N=int(lattice_size**2.)
        #Make the lattice with random numbers of S I R cells
        self.lattice=np.random.choice([-1,0,1],size=(lattice_size,lattice_size))

    def RandPoint(self):
        '''
        Generates a random coordinate
        within the lattice
        '''
        x=np.random.randint(0,self.lattice_size)
        y=np.random.randint(0,self.lattice_size)
        return (x,y)

    def FindNN(self,x,y):
        '''
        Finds the nearest neighbours of
        randomly pick spin at coordinate x,y
        '''
        nn_xl=(x-1)%self.lattice_size
        nn_xr=(x+1)%self.lattice_size
        nn_yu=(y+1)%self.lattice_size
        nn_yd=(y-1)%self.lattice_size
        nn = [self.lattice[nn_xl,y],self.lattice[nn_xr,y],
        self.lattice[x,nn_yu],self.lattice[x,nn_yd]]
        return(nn)

    @staticmethod
    def Acceptor(p):
        '''
        Static method to allow
        for probability calculations
        '''
        a=np.random.uniform(0,1)
        if p>a:return 1
        else: return 0

    def Rules(self,x,y,nnList):
        '''
        Susceptable: -1
        Infected: 0
        Recovered: 1
        '''
        #Counter takes in the list of nearest neighbours so that number of alive neighbours can be counted
        cnt=Counter(nnList)
        if self.lattice[x,y]== -1 and cnt[0]>=1 and SIRS.Acceptor(self.p1)==True:
            self.lattice[x,y]=0
        elif self.lattice[x,y]==0 and SIRS.Acceptor(self.p2)==True:
            self.lattice[x,y]=1
        elif self.lattice[x,y]==1 and SIRS.Acceptor(self.p3)==True:
            self.lattice[x,y]=-1

    def Update(self):
        x,y = self.RandPoint()
        nn=self.FindNN(x,y)
        self.Rules(x,y,nn)
        return self.lattice

    def Calc_Infected(self):
        '''
        self.array==0 sets elements to true for zero vales
        and np.count_nonzero will count all the true values
        which gives the number of cells which are infected
        even though infected is set to be 0
        '''
        infected=np.count_nonzero(self.lattice==0)
        #print('Number of infected cells: %.i'%(infected))
        return infected

    def Get_Variance(self,infected,infected_sq):
        av_infd=np.average(infected)
        av_infd_sq=np.average(infected_sq)
        variance=(av_infd_sq - av_infd**2. )/self.N
        return variance

    def FracImmune(self,fraction):
        '''
        Immune: 2
        randomly selects points in the
        lattice to become immune until
        a passed in fraction has been
        converted to have immunity
        '''
        indices=[]
        L=list(range(self.lattice_size))
        while len(indices)<int(self.N*fraction):
            x=random.choice(L)
            y=random.choice(L)
            item=[x,y]
            if item not in indices:
                indices.append(item)
        for i in indices:
            self.lattice[i[0],i[1]]=2
