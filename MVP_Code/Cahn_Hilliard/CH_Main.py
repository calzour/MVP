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
from CH_Classes import CH as CHc
import matplotlib.animation as animation
from multiprocessing import Queue, Process
from matplotlib.animation import FuncAnimation

def main():
    lattice_size=int(sys.argv[1])
    phi=float(sys.argv[2])    #set to zero
    sweeps=int(sys.argv[3])
    mode=str(sys.argv[4])
    dt,dx=1,1


    ch=CHc(dt,dx,lattice_size,phi)
    ch.Noise()
    lattice=ch.order_lattice


    if mode=='viz1':

        #Initiates a Queue to hold different lattices to be animated
        lattice_queue = Queue()
        lattice_queue.put( (copy.deepcopy(lattice)) )       #Copies the first random lattice to the queue
        animator = Animator(lattice_queue,sweeps)                  #Instances and initiates the Animator class
        animator_proc = Process(target=animator.animate)    #Sets up a process with animate as the target
        animator_proc.start()
        #Sweeps through the update rules for the lattice
        for i in range(sweeps):
            for j in range(ch.lattice_size):
                ch.Update()
            lattice_queue.put( copy.deepcopy(ch.order_lattice) )         #copies new lattice to Queue to be animated

    if mode=='viz2':
        #Update plot function which sweeps the array and simulates the SIRS model
        def UpdatePlot(*args):
            image.set_array(ch.order_lattice)
            ch.Update()
            return image,

        #Animates the simulation
        chImage = plt.figure()
        image = plt.imshow(ch.order_lattice,animated=True)
        model = FuncAnimation(chImage,UpdatePlot,interval=10,blit=True)
        plt.show()

    if mode =='data':

        fe_data=[]
        time_list=[]
        phi=0
        sweeps=60000
        #sweeps only needs to be about 50000
        for i in range(sweeps):
            ch.Update()
            if i>=500 and i %500 == 0:
                fe_data.append(ch.LatticeFreeEnergy())
                time_list.append(i)
                print(i)
        plt.plot(time_list,fe_data)
        plt.xlabel('Time')
        plt.ylabel('Free Energy')
        plt.title('phi = '+str(phi))
        plt.savefig('Graphs/Phi_0_Cahn',fmt='png')
        plt.show()

        fe_data=[]
        time_list=[]
        phi=0.5
        sweeps=60000

        ch=CHc(dt,dx,lattice_size,phi)
        ch.Noise()
        lattice=ch.order_lattice

        #sweeps only needs to be about 50000
        for i in range(sweeps):
            ch.Update()
            if i>=500 and i %500 == 0:
                fe_data.append(ch.LatticeFreeEnergy())
                time_list.append(i)
                print(i)
        plt.plot(time_list,fe_data)
        plt.xlabel('Time')
        plt.ylabel('Free Energy')
        plt.title('phi = '+str(phi))
        plt.savefig('Graphs/Phi_0_5_Cahn',fmt='png')
        plt.show()

        fe_data=[]
        time_list=[]
        phi= -0.5
        sweeps=60000

        ch=CHc(dt,dx,lattice_size,phi)
        ch.Noise()
        lattice=ch.order_lattice

        #sweeps only needs to be about 50000
        for i in range(sweeps):
            ch.Update()
            if i>=500 and i %500 == 0:
                fe_data.append(ch.LatticeFreeEnergy())
                time_list.append(i)
                print(i)
        plt.plot(time_list,fe_data)
        plt.xlabel('Time')
        plt.ylabel('Free Energy')
        plt.title('phi = '+str(phi))
        plt.savefig('Graphs/Phi_minus0_5_Cahn',fmt='png')
        plt.show()






if __name__=='__main__':
    main()
