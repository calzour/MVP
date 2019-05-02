import numpy as np
import sys
import time
import copy
from multiprocessing import Queue, Process
import Ising_Model_Classes as imc
import Ising_data as id
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from ClassAnimator import Animator

def main():
    '''
    Run code in following format:
    python3 < Main_Ising.py><lattice_size><T><G or K>
    where G=Glauber Dynamics and K=Kawasaki
    '''
    lattice_size=int(sys.argv[1])
    T=float(sys.argv[2])
    type=str(sys.argv[3])
    #Number of sweeps
    sweeps=10100
    #instance the class
    ising=imc.Lattice(lattice_size)
    #Create the lattice
    lattice=ising.MakeLattice()
    #Initiates a Queue to hold different lattices to be animated
    lattice_queue = Queue()
    lattice_queue.put( (copy.deepcopy(lattice)) )       #Copies the first random lattice to the queue
    animator = Animator(lattice_queue)                  #Instances and initiates the Animator class
    animator_proc = Process(target=animator.animate)    #Sets up a process with animate as the target
    animator_proc.start()                               #Starts the animation
    #Perform dynamics
    for i in range(sweeps):
        for j in range(lattice_size*lattice_size):
            ising.Dynamics(lattice,T,type)
        lattice_queue.put( copy.deepcopy(lattice) )         #copies new lattice to Queue to be animated
if __name__ == "__main__":
    main()
