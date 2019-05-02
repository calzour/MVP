import sys
import numpy as np
import Ising_data as id
import matplotlib.pyplot as plt
import Ising_Model_Classes as imc
from datetime import date
import os

def main():
    '''
    Run code in following format:
    python3 < Ising_Data_Main.py><lattice_size><G or K><tag>
    where G=Glauber Dynamics and K=Kawasaki
    Note that T has been removed for this file compared to the Main_Ising file,
    this is because I want to hard-code the temperature range for data taking
    <tag> is a string which will be added to the graph savefig names to distinguish them if wanted
    '''
    lattice_size=int(sys.argv[1])
    type=str(sys.argv[2])
    tag=str(sys.argv[3])

    #Number of sweeps
    sweeps=10100
    M_list=[]                   #list for containing total magnetisation of a microstate
    M_av_list=[]
    T_list=[]                   #List containing the temperatures when magnetisation was taken
    chi_list=[]
    Energy_list=[]
    Energy_av_list=[]
    heatcap_list=[]
    M_errors=[]
    E_errors=[]
    chi_errors=[]
    c_errors=[]

    #instance the class
    ising=imc.Lattice(lattice_size)             #not used
    idata=id.IsingDataClass(lattice_size)

    #Create the lattice
    if type=='G':
        lattice=np.ones((lattice_size,lattice_size))
        print(lattice)
    elif type=='K':
        '''
        Makes a lattice with half ones and half -ones
        '''
        lattice1=np.ones((lattice_size,lattice_size))
        #print(lattice1)
        lattice2=np.full(((lattice_size,lattice_size)),-1.0)
        #print(lattice2)
        lattice=np.append(lattice1, lattice2,axis=0)
        lattice=lattice[1::2]
        print(lattice)

    for T in np.linspace(1.0,3.0,10):       #start, end, number of steps
        M_list=[]
        Energy_list=[]
        T_list.append(T)
        print(T)
        for i in range(sweeps):
            if i%500 == 0:
                print(i)
            for j in range(lattice_size*lattice_size):
                ising.Dynamics(lattice,T,type)                 #updates the lattice/applies dynamics
            if i%10 ==0 and i>100:
                M_list.append(idata.TotalMag(lattice))
                Energy_list.append(idata.TotalE(lattice))

        #Make measurements for <M> and the suceptability
        M_av_list.append(abs(idata.AvM(M_list)))                    #each append adds my <M> value
        M_errors.append(idata.ErrorSTD(M_list))                     #Appends the standard deviation of the magnetisation
        print('Average M is: '+str(idata.AvM(M_list)))

        #Energy measurements
        Energy_av_list.append(idata.AvM(Energy_list))               #Calcs and adds average energy to a list <E>
        E_errors.append(idata.ErrorSTD(Energy_list))

        #Calculates the susceptability and heat capacity plus their errors
        chi_list.append(idata.Susceptability(M_list))           #Want to also calculate the suceptability here too
        heatcap_list.append(idata.HeatCapacity(Energy_list,T))
        chi_errors.append(idata.bootstrapErrors(M_list,1.0))
        c_errors.append(idata.bootstrapErrors(Energy_list,T**2.0))


    #makes a new folder for putting data into if one doesnt exist already
    try:
        os.makedirs(tag)
    except OSError:
        if not os.path.isdir(tag):
            raise+type
    
    #M vs. T
    plt.plot(T_list,np.abs(M_av_list))
    plt.title('Magnetisation vs. Temperature')
    plt.ylabel('Average Magnetisation')
    plt.xlabel('Temperature')
    plt.errorbar(T_list, np.abs(M_av_list),yerr=M_errors,fmt='k-',capsize=2)
    plt.savefig(tag+'/Mag_v_T_'+str(lattice_size)+'_'+str(date.today())+'_'+tag+'_'+type ,format='pdf')
    #plt.show()
    plt.clf()

    #Chi vs. T
    plt.plot(T_list,chi_list)
    plt.title('Susceptability vs. Temperature')
    plt.ylabel('Susceptability')
    plt.xlabel('Temperature')
    plt.errorbar(T_list,chi_list,yerr=chi_errors,fmt='k-',capsize=2)
    plt.savefig(tag+'/Sus_v_T_'+str(lattice_size)+'_'+str(date.today())+'_'+tag+'_'+type,format='pdf')
    #plt.show()
    plt.clf()

    #Energy vs T
    plt.plot(T_list,Energy_av_list)
    plt.title('Temperature vs. Energy')
    plt.ylabel('Average Energy')
    plt.xlabel('Temperature')
    plt.errorbar(T_list,Energy_av_list,yerr=E_errors,fmt='k-',capsize=2)
    plt.savefig(tag+'/E_v_T_'+str(lattice_size)+'_'+str(date.today())+'_'+tag+'_'+type,format='pdf')
    #plt.show()
    plt.clf()


    #Heat capacity vs T
    plt.plot(T_list,heatcap_list)
    plt.title('Temperature vs. Heat Capacity')
    plt.ylabel('Heat Capacity')
    plt.xlabel('Temperature')
    plt.errorbar(T_list,heatcap_list,yerr=c_errors,fmt='k-',capsize=2)
    plt.savefig(tag+'/C_v_T_'+str(lattice_size)+'_'+str(date.today())+'_'+tag+'_'+type,format='pdf')
    #plt.show()
    plt.clf()


if __name__ == "__main__":
    main()
