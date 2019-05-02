import numpy as np
import sys
import time
import matplotlib.animation as animation
import matplotlib.pyplot as plt


class Lattice():
    #Constructor
    def __init__(self,lattice_size=50,J=1.0,kb=1.0):
        self.J=J
        self.kb=kb
        self.lattice_size=lattice_size

    def MakeLattice(self):
        '''
        Creates a random initital lattice
        where values are either -1 or 1
        #line of code bellow seems to slow down the animation for some reason
        lattice=np.random.choice((-1,1),(self.lattice_size,self.lattice_size))
        '''
        lattice= np.random.uniform(-1,1,size=(self.lattice_size,self.lattice_size))
        for i,row in enumerate(lattice):
            for j,number in enumerate(row):
                if number<0:
                    lattice[i,j]=np.sign(-1)
                elif number>0:
                    lattice[i,j]=np.sign(+1)
        return lattice

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
        nn = [[nn_xl,y],[nn_xr,y],[x,nn_yu],[x,nn_yd]]
        return(nn)

    def LocalEnergy(self,x,y,nn,lattice):
        '''
        Calculates the energy of the nearest
        neighbours for a single point in the
        lattice and return the total energy
        of the configuration.
        '''
        point = lattice[x,y]
        E=0
        for i,coord in enumerate(nn):
            nx=coord[0]
            ny=coord[1]
            E+=point*lattice[nx,ny]
        E=-self.J*E
        return E

    def Acceptor(self,lattice,energy_i,energy_f,T,x,y,x2=None,y2=None,type='G'):
        '''
        Metroplois algorithm
        '''
        del_E = energy_f - energy_i
        if del_E <= 0:                              #Accepts switch
            pass
        elif del_E >0:                              #Rejects or accepts based on probability
            dart = np.random.uniform(0,1)
            if dart <= np.exp(-(del_E/(self.kb*T))):
                pass
            elif type=='G':
                lattice[x,y]=lattice[x,y]*(-1.0)
            elif type=='K':
                lattice[x,y]=(lattice[x,y])*(-1.0)
                lattice[x2,y2]=(lattice[x2,y2])*(-1.0)
        return lattice


    def Dynamics(self,lattice,T,type='G'):
        if type=='G':
            '''
            Glauber Dynamics
            '''
            x,y = self.RandPoint()
            nn=self.FindNN(x,y)
            energy_i =self.LocalEnergy(x,y,nn,lattice)
            lattice[x,y]=(lattice[x,y])*(-1.0)          #Flips the spin to test
            energy_f = self.LocalEnergy(x,y,nn,lattice) #re-calculates the energy
            lattice= self.Acceptor(lattice,energy_i,energy_f,T,x,y)
            return lattice
        if type=='K':
            '''
            Kawasaki Dynamics-
            Chose the two consecutive spin flips as
            I was running into a lot of problems
            with the attempt to swap the spins and
            this method ended up being easier to implement
            which is why I chose this method.

            From feedback, I don't need to flip the spin, test then flip back if
            not needed, I could just do calcuation of new energy without having
            to actually flip the spin 
            '''
            x,y=self.RandPoint()
            x2,y2 = self.RandPoint()
            #If both same points, chooses another point
            while lattice[x, y] == lattice[x2, y2]:
                x2,y2 = self.RandPoint()
            nn,nn2=self.FindNN(x,y),self.FindNN(x2,y2)
            if [x2,y2] in nn:
                energy_i=self.LocalEnergy(x,y,nn,lattice)+self.LocalEnergy(x2,y2,nn2,lattice) + 4.0*self.J*-1.0
                lattice[x,y],lattice[x2,y2]=(lattice[x,y])*(-1.0),(lattice[x2,y2])*(-1.0)
                energy_f =self.LocalEnergy(x,y,nn,lattice) + self.LocalEnergy(x2,y2,nn2,lattice) +4.0*self.J*-1.0
                lattice=self.Acceptor(lattice,energy_i,energy_f,T,x,y,x2,y2,type='K')
            else:
                energy_i=self.LocalEnergy(x,y,nn,lattice)+self.LocalEnergy(x2,y2,nn2,lattice)
                lattice[x,y],lattice[x2,y2]=(lattice[x,y])*(-1.0),(lattice[x2,y2])*(-1.0)
                energy_f=self.LocalEnergy(x,y,nn,lattice) + self.LocalEnergy(x2,y2,nn2,lattice)
                lattice=self.Acceptor(lattice,energy_i,energy_f,T,x,y,x2,y2,type='K')
            return lattice
