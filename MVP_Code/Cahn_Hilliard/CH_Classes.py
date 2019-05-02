import copy
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from scipy import ndimage
import warnings
import random


class CH():
    def __init__(self,dt,dx,lattice_size,phi):
        self.dt=dt
        self.dx=dx
        self.lattice_size=lattice_size
        self.N=lattice_size*lattice_size
        self.phi=phi
        self.order_lattice=np.full((self.lattice_size,self.lattice_size),phi)
        self.chem_lattice=np.zeros((self.lattice_size,self.lattice_size))
        self.a =float(0.1)
        self.k = float(0.1)
        self.M = float(0.1)

    def PBC(self,index):
        '''
        Ensure periodic boundary conditions
        '''
        return index%(self.lattice_size)

    def Noise(self):
        '''
        Apply noise to the order lattice
        '''
        for i in range(self.lattice_size):
            for j in range(self.lattice_size):
                self.order_lattice[i,j]=self.order_lattice[i,j]+ random.uniform(-0.1,0.1)

    def UpdateChemSite(self,i,j):
        '''
        Updates the Chemical potential site
        '''
        self.chem_lattice[i,j] = -self.a*(self.order_lattice[i,j]) + self.a*(self.order_lattice[i,j])**3. - (self.k/self.dx**2.)*(self.order_lattice[self.PBC(i+1),j]+self.order_lattice[self.PBC(i-1),j]+self.order_lattice[i,self.PBC(j+1)] + self.order_lattice[i,self.PBC(j-1)]-4*self.order_lattice[i,j])

    def CalcChemPot(self):
        '''
        Updates the chemical potential site
        '''
        for i in range(self.lattice_size):
            for j in range(self.lattice_size):
                self.UpdateChemSite(i,j)

    def UpdateOrderLattice(self):
        '''
        Updates the order lattice sites
        '''
        CONSTANT=(self.dt * self.M)/(self.dx**2.0)
        for i in range(self.lattice_size):
            for j in range(self.lattice_size):
                self.order_lattice[i,j] = self.order_lattice[i,j] + CONSTANT*(self.chem_lattice[self.PBC(i+1),j]+self.chem_lattice[self.PBC(i-1),j]+self.chem_lattice[i,self.PBC(j+1)]+self.chem_lattice[i,self.PBC(j-1)] - 4*self.chem_lattice[i,j])

    def Update(self):
        '''
        Run the updating of the lattices
        '''
        self.CalcChemPot()
        self.UpdateOrderLattice()

    def DotProduct(self,list1,list2):
        '''
        Takes the dot product
        '''
        dotlist=[]
        if len(list1) == len(list2):
            for i in range(len(list1)):
                dotlist.append(list1[i]*list2[i])
                return sum(dotlist)

    def DelSq(self,array,i,j):
        '''
        Performs Del Squared operation on two elements in an array
        '''
        axisx=(array[self.PBC(i+1),j]-array[self.PBC(i-1),j])/(2.0*self.dx)
        axisy=(array[i,self.PBC(j-1)]-array[i,self.PBC(j-1)])/(2.0*self.dx)
        dotproduct=self.DotProduct([axisx,axisy],[axisx,axisy])
        return dotproduct

    def FreeEnergyDensity(self,i,j):
        '''
        Calculates the free energy density
        '''
        fed=(-self.a/2.0) * (self.order_lattice[i,j])**2.0 + (self.a/4.0)*(self.order_lattice[i,j])**4.0 + (self.k/2.0)*(self.DelSq(self.order_lattice,i,j))
        return fed

    def LatticeFreeEnergy(self):
        '''
        Calculates the free energy for the entire lattice
        '''
        fe_list=[]
        for i in range(self.lattice_size):
            for j in range(self.lattice_size):
                fe_list.append(self.FreeEnergyDensity(i,j))
        return sum(fe_list)
