import sys
import time
import copy
import numpy as np
from GoL_Class import GoL as GoL
import matplotlib.pyplot as plt
from collections import Counter
from ClassAnimator import Animator
import matplotlib.animation as animation
from multiprocessing import Queue, Process
from matplotlib.animation import FuncAnimation

def main():
    '''
    Run code in following format:
    python3 <GoL_Main.py><lattice_size><R or A or B or G><number of sweeps><viz or data>
    where R=Random, A=Absorbing, B=Blinker, G=Glider
    '''

    #Set up the sys arguments
    lattice_size=int(sys.argv[1])
    state=str(sys.argv[2])
    sweeps=int(sys.argv[3])
    type=str(sys.argv[4])



    #Instance the class
    GoLClass=GoL(lattice_size,state)
    lattice=GoLClass.lattice

    if type=='viz':
        #Update plot function which sweeps the array and simulates the GOL
        def UpdatePlot(*args):
            image.set_array(GoLClass.lattice)
            GoLClass.Update_Lattice()
            return image,

        #Animating the simulation. plt.imshow creates much faster moving animations so it is used here.
        GLImage = plt.figure()
        image = plt.imshow(GoLClass.lattice,animated=True)
        model = FuncAnimation(GLImage,UpdatePlot,interval=10,blit=True)
        plt.show()

    if type=='data':
        #Get data for CoM and velocity
        if state=='G':
            time=0
            #required to make it work
            sweeps=150
            for i in range(sweeps):
                GoLClass.Update_Lattice()
                time_list,CoM_list,time=GoLClass.Get_CoM(time)

            print('Velocity of the Glider is: %.3f'%(GoLClass.Get_Velocity(time_list,CoM_list)))
            plt.plot(time_list,CoM_list)
            plt.ylabel('CoM position')
            plt.xlabel('Time')
            plt.show()





if __name__=='__main__':
    main()
