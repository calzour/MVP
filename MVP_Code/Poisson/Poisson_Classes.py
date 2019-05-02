import copy
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from scipy import ndimage
import warnings
import random
import pickle
from datetime import date

class PC():
    def __init__(self,dt,dx,lattice_size,psi):
        self.dt=dt
        self.dx=dx
        self.lattice_size=lattice_size
        self.N=lattice_size*lattice_size*lattice_size
        self.e=float(1.0)
        self.psi_lattice=np.full((self.lattice_size,self.lattice_size,self.lattice_size),psi)
        #self.previous_psi_lattice=np.full((self.lattice_size,self.lattice_size,self.lattice_size),psi)
        self.rho_lattice=np.zeros((self.lattice_size,self.lattice_size,self.lattice_size))
        self.halfway=int(self.lattice_size/2.0)


    def PointCharge(self):
        for i,value in enumerate(self.rho_lattice):
            for j in value:
                for k in j:
                    self.rho_lattice[self.halfway][self.halfway][self.halfway]=1.0
        return self.rho_lattice


    def UpdatePsi1(self):
        '''
        Changes element one at a time
        Guass-Siedal
        '''
        for i in range(1,self.lattice_size-1):
            for j in range(1,self.lattice_size-1):
                for k in range(1,self.lattice_size-1):
                    '''
                    self.psi_lattice[i][j][k]=(1.0/6.0)*(self.psi_lattice[i-1][j][k] + self.psi_lattice[i+1][j][k] + self.psi_lattice[i][j-1][k] + self.psi_lattice[i][j+1][k] + self.psi_lattice[i][j][k-1] + self.psi_lattice[i][j][k+1] + self.rho_lattice[i][j][k])
                    '''
                    a=self.psi_lattice[i-1][j][k] +self.psi_lattice[i+1][j][k]
                    b=self.psi_lattice[i][j-1][k] +self.psi_lattice[i][j+1][k]
                    c=self.psi_lattice[i][j][k-1] +self.psi_lattice[i][j][k+1] +self.rho_lattice[i][j][k]
                    self.psi_lattice[i][j][k]=(1.0/6.0)*(a+b+c)
        return self.psi_lattice


    def UpdatePsi2(self):
        '''
        Update the psi lattice all at once
        Jacobi
        '''
        temp_lattice=np.zeros((self.lattice_size,self.lattice_size,self.lattice_size))
        for i in range(1,self.lattice_size-1):
            for j in range(1,self.lattice_size-1):
                for k in range(1,self.lattice_size-1):
                    temp_lattice[i][j][k]=(1.0/6.0)*(self.psi_lattice[i-1][j][k] + self.psi_lattice[i+1][j][k] + self.psi_lattice[i][j-1][k] + self.psi_lattice[i][j+1][k] + self.psi_lattice[i][j][k-1] + self.psi_lattice[i][j][k+1] + self.rho_lattice[i][j][k])
        self.psi_lattice=temp_lattice
        return self.psi_lattice

    def Update1(self,threshold):
        '''
        Guass-Siedal
        '''
        self.previous_psi_lattice = np.copy(self.psi_lattice)
        self.psi_lattice = np.copy(self.UpdatePsi1())
        diff_array=np.absolute(self.previous_psi_lattice - self.psi_lattice)
        diff=np.sum(diff_array)
        while diff > threshold:
            self.previous_psi_lattice = np.copy(self.psi_lattice)
            self.psi_lattice = np.copy(self.UpdatePsi1())
            diff_array=np.absolute(self.previous_psi_lattice - self.psi_lattice)
            diff=np.sum(diff_array)
        print(diff)
        return self.psi_lattice


    def Update2(self,threshold):
        '''
        Jacobi
        '''
        self.previous_psi_lattice = np.copy(self.psi_lattice)
        self.psi_lattice = np.copy(self.UpdatePsi2())
        diff_array=np.absolute(self.previous_psi_lattice - self.psi_lattice)
        diff=np.sum(diff_array)
        while diff > threshold:
            self.previous_psi_lattice = np.copy(self.psi_lattice)
            self.psi_lattice = np.copy(self.UpdatePsi1())
            diff_array=np.absolute(self.previous_psi_lattice - self.psi_lattice)
            diff=np.sum(diff_array)
        print(diff)
        return self.psi_lattice

    def PlotHeatMap(self):
        rho_lattice=np.copy(self.rho_lattice)
        value_rho=rho_lattice[:][:][self.halfway]
        x = y = np.linspace(0,self.lattice_size-1,self.lattice_size)
        X,Y=np.meshgrid(x,y)
        plt.contourf(X,Y,value_rho,30)
        plt.xlabel('x-axis')
        plt.ylabel('y-axis')
        plt.title('Contour plot of Rho Lattice')
        plt.colorbar();
        plt.savefig('Graphs/RhoContourPlot',fmt='png')
        plt.clf()

    def GetEField(self,psi_lattice,tag):
        '''
        Only call this method once the
        Poisson Equation has been solved
        '''
        E_field=np.gradient(psi_lattice)
        with open('E_Field_'+tag+str(date.today())+'_'+str(self.lattice_size),'wb') as f:
            pickle.dump(E_field,f)
        norm=1.0*np.sum(E_field[:])
        u=E_field[0]/norm
        v=E_field[1]/norm
        q=E_field[2]/norm
        return (u,v,q)

    def GetBField(self,psi_lattice):
        bi, bj = np.gradient(psi_lattice[:][:][1])
        B_field=np.array([bi,-bj,np.zeros((self.lattice_size,self.lattice_size))])
        normB=-1*np.sqrt((B_field[0])**2.0 + (B_field[1])**2.0)
        B_fieldnorm=B_field/normB
        return B_field, normB


    def LogLogPlot(self,psi_lattice):
        r_list=[]
        pot_list=[]
        for i in range(self.lattice_size):
            for j in range(self.lattice_size):
                for k in range(self.lattice_size):
                    r=np.sqrt(((i-self.halfway)**2.0)+((j-self.halfway)**2.0)+((k-self.halfway)**2.0))
                    r_list.append(np.log(r))
                    pot=psi_lattice[i][j][k]
                    pot_list.append(np.log(pot))
        return (r_list,pot_list)

    def OverrelaxUpdate(self,w):
        '''
        Changes element one at a time
        Guass-Siedal
        '''
        for i in range(1,self.lattice_size-1):
            for j in range(1,self.lattice_size-1):
                for k in range(1,self.lattice_size-1):
                    '''
                    self.psi_lattice[i][j][k]=(1.0/6.0)*(self.psi_lattice[i-1][j][k] + self.psi_lattice[i+1][j][k] + self.psi_lattice[i][j-1][k] + self.psi_lattice[i][j+1][k] + self.psi_lattice[i][j][k-1] + self.psi_lattice[i][j][k+1] + self.rho_lattice[i][j][k])
                    '''
                    a=self.psi_lattice[i-1][j][k] +self.psi_lattice[i+1][j][k]
                    b=self.psi_lattice[i][j-1][k] +self.psi_lattice[i][j+1][k]
                    c=self.psi_lattice[i][j][k-1] +self.psi_lattice[i][j][k+1] +self.rho_lattice[i][j][k]
                    self.psi_lattice[i][j][k]= (1-w)*self.psi_lattice[i][j][k] +(1.0/6.0)*(a+b+c)*w
        return self.psi_lattice

    def OverRelax(self,n,threshold,psi):
        '''
        Guass-Siedal
        '''
        self.psi_lattice=np.full((self.lattice_size,self.lattice_size,self.lattice_size),psi)
        w_list=[]
        iterations_list=[]
        self.rho_lattice=self.PointCharge()
        for w in np.linspace(1.,1.98,n):
            iterations=0
            print(w)
            self.psi_lattice=np.full((self.lattice_size,self.lattice_size,self.lattice_size),psi)
            w_list.append(w)
            self.previous_psi_lattice = np.copy(self.psi_lattice)
            self.psi_lattice = np.copy(self.OverrelaxUpdate(w))
            diff_array=np.absolute(self.previous_psi_lattice - self.psi_lattice)
            diff=np.sum(diff_array)
            while diff > threshold:
                iterations+=1
                self.previous_psi_lattice = np.copy(self.psi_lattice)
                self.psi_lattice = np.copy(self.OverrelaxUpdate(w))
                diff_array=np.absolute(self.previous_psi_lattice - self.psi_lattice)
                diff=np.sum(diff_array)
            iterations_list.append(iterations)
            print(diff)
        return (self.psi_lattice, w_list, iterations_list)
