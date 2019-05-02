import numpy as np
import sys
import time
import copy
from multiprocessing import Queue, Process
import Ising_Model_Classes as imc
from Ising_Model_Classes import Lattice
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from ClassAnimator import Animator


class IsingDataClass():
    def __init__(self,lattice_size,kb=1.0,J=1.0):
        '''
        Initialises the class with the lattice size and constants
        '''
        self.lattice_size=lattice_size
        self.N=lattice_size*lattice_size
        self.kb=kb
        self.J=J

    def TotalMag(self,lattice):
        '''
        Sums every spin in the lattice
        '''
        M=0                      #initialise M for magnitisation
        for i in lattice:
            for j in i:
                M+=j
        return M

    def AvM(self,M_list):
        '''
        Calculates the average value from a list
        '''
        avM=np.average(M_list)
        return avM

    def Susceptability(self,M_list):
        '''
        Calculates the susceptability
        '''
        chi= (1.0/(self.N * self.kb))*np.var(M_list)
        return chi

    def TotalE(self,lattice):
        '''
        Iterates through the entire list and
        calculates the energy of the lattice
        '''
        E=0
        for i,row in enumerate(lattice):
            for j,value in enumerate(row):
                nn=Lattice.FindNN(self,i,j)
                E+=Lattice.LocalEnergy(self,i,j,nn,lattice)
        return(E/2.0)

    def HeatCapacity(self,Total_E_list,T):
        '''
        Calculates the heat capacity from a list
        '''
        C= (1.0/(self.N * self.kb *T**2.0))*(np.var(Total_E_list))
        return C

    def ErrorSTD(self,list):
        '''
        Calculates the standard error on the mean
        '''
        error=np.std(list)/np.sqrt(len(list))
        return error

    def varianceOfList(self,nplist,T):
        '''
        Calculates the variance of a list of numbers
        '''
        listaverage = np.mean(nplist)               #calcs average of a list
        averageXAllSquared = (listaverage **2)      #squares the average value of the list
        tempList = []                               #empty list
        for i in range(len(nplist)):                #loops through the range of the list
            tempList.append( (nplist[i]) **2)       #creates a list where each value is the origonal value squared
        xSquaredAveraged = np.mean(tempList)        ##calcs the average of the squared values
        variance = (1.0/((self.N)**2.0)*self.kb * T)*(xSquaredAveraged - averageXAllSquared )       #subtract the difference
        return(variance, listaverage)               #returns the variance and the average of the origonal list

    def bootstrapErrors(self,List,T):
        '''
        Calculates the bootstrap errors by sampling randomly from a list
        '''
        #calculate parameters
        noSamples = 100
        listSize = len(List)
        sampleList = [] #empty list for sampling
        obsList = [] #list for variances of sample lists
        variance = 0.0 #the error that the function returns
        dummy = 0.0     #dummy variable for variance function
        for p in range(noSamples): #create 100 random samples
            sampleList = np.random.choice( List, listSize )         #randomly selects a listSize number of elements from the List
            tempVariance, dummy = self.varianceOfList(sampleList,T) #calcs variance of random sample and saves to a list
            obsList.append(tempVariance)
        #after 100 times, find the variance of this new list
        varianceSq, dummy = self.varianceOfList(obsList,T)          #calcs the variance of the randomly sample variances
        variance = np.sqrt(varianceSq)                              #returns the actualy variance (takes sqrt) (maybe not needed?)
        return(variance)
