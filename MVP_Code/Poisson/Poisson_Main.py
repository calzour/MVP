import os
import sys
import time
import copy
import pickle
import numpy as np
from datetime import date
import matplotlib.pyplot as plt
from collections import Counter
from ClassAnimator import Animator
from Poisson_Classes import PC as PCc
import matplotlib.animation as animation
from multiprocessing import Queue, Process
from matplotlib.animation import FuncAnimation

def main():

    lattice_size=int(sys.argv[1])
    threshold=float(sys.argv[2])
    #mode=str(sys.argv[3])

    #threshold=0.00001
    #lattice_size=15
    psi=0.0
    dt,dx=1.0,1.0

    halfway=int(lattice_size/2.0)

    pc=PCc(dt,dx,lattice_size,psi)
    rho_lattice=pc.PointCharge()


    pc.PointCharge()
    pc.PlotHeatMap()


    psi_lattice=pc.Update1(threshold)
    with open('psi_lattice1_'+str(date.today())+'_'+str(lattice_size),'wb') as f:
        pickle.dump(psi_lattice,f)


    u,v,q= pc.GetEField(psi_lattice,'Gauss')

    #Save data to output file
    fileu1=open('UE_field_Gauss'+str(date.today())+'_'+str(lattice_size),'w')
    filev1=open('VE_field_Gauss'+str(date.today())+'_'+str(lattice_size),'w')
    fileq1=open('QE_field_Gauss'+str(date.today())+'_'+str(lattice_size),'w')
    for i in range(lattice_size):
        for j in range(lattice_size):
            for k in range(lattice_size):
                fileu1.write('{} {} {} {}\n'.format(i,j,k,(u[i,j,k])))
                filev1.write('{} {} {} {}\n'.format(i,j,k,(v[i,j,k])))
                fileq1.write('{} {} {} {}\n'.format(i,j,k,(q[i,j,k])))
    fileu1.close()
    filev1.close()
    fileq1.close()

    #plot the psi lattice after updating
    plt.imshow(psi_lattice[halfway][:][:])
    plt.quiver(q[halfway][:][:],v[halfway][:][:],angles='xy')
    plt.xlabel('$E_{x}$ Field')
    plt.ylabel('$E_{y}$ Field')
    plt.title('Gauss Electric Field Plot')
    plt.savefig('Graphs/E_Field_Gauss'+str(date.today())+'_'+str(lattice_size),fmt='png')
    plt.clf()

    r_list,pot_list=pc.LogLogPlot(psi_lattice)

    with open('r_list1'+str(date.today())+'_'+str(lattice_size),'wb') as f:
        pickle.dump(r_list,f)
    with open('pot_list1'+str(date.today())+'_'+str(lattice_size),'wb') as f:
        pickle.dump(pot_list,f)

    plt.scatter(r_list,pot_list,s=2,c='k')
    plt.xlabel('log(distance)')
    plt.ylabel('log(psi)')
    plt.title('Gauss-Seidal Method')
    plt.savefig('Graphs/loglogplot_Gauss'+str(date.today())+'_'+str(lattice_size),fmt='png')
    plt.clf()






    ############################################################################
    psi_lattice= pc.Update2(threshold)
    with open('psi_lattice2'+str(date.today())+'_'+str(lattice_size),'wb') as f:
        pickle.dump(psi_lattice,f)

    #Save data to output file
    fileu2=open('UE_field_Jacobi'+str(date.today())+'_'+str(lattice_size),'w')
    filev2=open('VE_field_Jacobi'+str(date.today())+'_'+str(lattice_size),'w')
    fileq2=open('QE_field_Jacobi'+str(date.today())+'_'+str(lattice_size),'w')
    for i in range(lattice_size):
        for j in range(lattice_size):
            for k in range(lattice_size):
                fileu2.write('{} {} {} {}\n'.format(i,j,k,(u[i,j,k])))
                filev2.write('{} {} {} {}\n'.format(i,j,k,(v[i,j,k])))
                fileq2.write('{} {} {} {}\n'.format(i,j,k,(q[i,j,k])))
    fileu2.close()
    filev2.close()
    fileq2.close()


    u,v,q= pc.GetEField(psi_lattice,'Jacobi')


    #plot the psi lattice after updating
    plt.imshow(psi_lattice[halfway][:][:])
    plt.quiver(q[halfway][:][:],v[halfway][:][:],angles='xy')
    plt.xlabel('$E_{x}$ Field')
    plt.ylabel('$E_{y}$ Field')
    plt.title('Jacobi Electric Field Plot')
    plt.savefig('Graphs/E_Field_Jacobi'+str(date.today())+'_'+str(lattice_size),fmt='png')
    plt.clf()

    r_list,pot_list=pc.LogLogPlot(psi_lattice)

    with open('r_list2'+str(date.today())+'_'+str(lattice_size),'wb') as f:
        pickle.dump(r_list,f)
    with open('pot_list2'+str(date.today())+'_'+str(lattice_size),'wb') as f:
        pickle.dump(pot_list,f)


    plt.scatter(r_list,pot_list,s=2,c='k')
    plt.xlabel('log(distance)')
    plt.ylabel('log(psi)')
    plt.title('Jacobi Method')
    plt.savefig('Graphs/loglogplot_Jacobi'+str(date.today())+'_'+str(lattice_size),fmt='png')
    plt.clf()


    ############################################################################

    #reinstance class
    pc=PCc(dt,dx,lattice_size,psi)
    rho_lattice=pc.rho_lattice
    rho_lattice[:][:][1]=1.0  #make a wire
    #Rerun the updateing of the psi lattice
    psi_lattice=pc.Update1(threshold)
    #print(rho_lattice)
    plt.imshow(psi_lattice[:][:][1])
    plt.clf()

    B_field,normB=pc.GetBField(psi_lattice)
    plt.imshow(psi_lattice[:][:][1])
    plt.quiver(B_field[0]/normB,B_field[1]/normB,angles='xy',scale=None)
    plt.xlabel('$B_{x}$ Field')
    plt.ylabel('$B_{y}$ Field')
    plt.title('Gauss Magnetic Field Plot')
    plt.savefig('Graphs/B_Field_Gauss'+str(date.today())+'_'+str(lattice_size),fmt='png')
    plt.show()


    ############################################################################
    pc=PCc(dt,dx,lattice_size,psi)
    rho_lattice=pc.PointCharge()
    n=20
    psi_lattice,w_list,iterations_list = pc.OverRelax(n,threshold,psi)
    plt.plot(w_list,iterations_list)
    plt.xlabel('w')
    plt.ylabel('Iterations before Convergence')
    plt.savefig('Graphs/RelaxationThing',fmt='png')
    plt.show()











if __name__=='__main__':
    main()
