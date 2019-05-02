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
from SIRS_Class import SIRS as sirsc
import matplotlib.animation as animation
from multiprocessing import Queue, Process
from matplotlib.animation import FuncAnimation

def main():
    '''
    Run code in following format:
    python3 <SIRS_Main.py><lattice_size><p1><p2><p3><number of sweeps><mode><state>

    state can be: user, waves, absord or balanced
    mode can be: viz1, viz2, phase, slice or immune
    '''
    #Take in sys arguments
    lattice_size=int(sys.argv[1])
    p1=float(sys.argv[2])
    p2=float(sys.argv[3])
    p3=float(sys.argv[4])
    sweeps=int(sys.argv[5])
    mode=str(sys.argv[6])
    if str(sys.argv[7])=='user':
        pass
    elif str(sys.argv[7])=='waves':
        print('Setting probabilities to allow for waves.')
        p1=0.8
        p2=0.095
        p3=0.01
    elif str(sys.argv[7])=='absorb':
        #should result in everything being infected
        print('Setting probabilities to allow for all cells infected.')
        p1=0.
        p2=0.1
        p3=1.
    elif str(sys.argv[7])=='balanced':
        print('Setting probabilities to allow for dynamic equilibrium.')
        p1=0.5
        p2=0.5
        p3=0.5
    else:
        pass

    #Instance the class
    sirs=sirsc(p1,p2,p3,lattice_size)
    lattice=sirs.lattice

################################################################################

    #Worse animation
    if mode=='viz1':
        #Update plot function which sweeps the array and simulates the SIRS model
        def UpdatePlot(*args):
            image.set_array(sirs.lattice)
            sirs.Update()
            return image,

        #Animates the simulation
        sirsImage = plt.figure()
        image = plt.imshow(sirs.lattice,animated=True)
        model = FuncAnimation(sirsImage,UpdatePlot,interval=10,blit=True)
        plt.show()

################################################################################

    #Contour method of visulisation (worse method, runs in parallel though)
    if mode=='viz2':

        #Initiates a Queue to hold different lattices to be animated
        lattice_queue = Queue()
        lattice_queue.put( (copy.deepcopy(lattice)) )       #Copies the first random lattice to the queue
        animator = Animator(lattice_queue,sweeps)                  #Instances and initiates the Animator class
        animator_proc = Process(target=animator.animate)    #Sets up a process with animate as the target
        animator_proc.start()
        #Sweeps through the update rules for the lattice
        for i in range(sweeps):
            for j in range(sirs.N):
                sirs.Update()
            lattice_queue.put( copy.deepcopy(sirs.lattice) )         #copies new lattice to Queue to be animated

################################################################################

    if mode=='phase':
        #set p2 to a constant of 0.5
        p2=0.5
        p1=0.0
        p3=0.0
        #100 sweeps to equilibriate, 1000 for taking measurements
        sweeps=1100
        av_infdlist=[]
        var_infd=[]
        steps=0.05
        p1list=np.arange(p1,1.,steps)
        p3list=np.arange(p3,1.,steps)

        tag=mode+str(sweeps)
        try:
            os.makedirs(tag)
        except OSError:
            if not os.path.isdir(tag):
                raise+type

        #The main loop for data collection
        for a in np.arange(p1, 1.,steps):
            #print('p1: %.3f'%(a))
            for s in np.arange(p3,1.,steps):
                #print('p3: %.3f \n'%(s))
                infected=[]
                infected_sq=[]
                #Reinstance class with updated probabilities
                sirs=sirsc(a,p2,s,lattice_size)
                print('p1: %.3f'%(sirs.p1))
                print('p3: %.3f\n'%(sirs.p3))
                for i in range(sweeps):
                    for j in range(sirs.N):
                        sirs.Update()
                    if i%10==0 and i>100:
                        inf=sirs.Calc_Infected()
                        infected.append(inf)
                        infected_sq.append(inf**2.0)
                av_infd=(np.average(infected))/sirs.N
                av_infdlist.append(av_infd)
                var_infd.append(sirs.Get_Variance(infected,infected_sq))

        #Pickling
        with open(tag+'/p1list2_'+str(sweeps)+'_'+str(steps)+'_'+str(date.today()),'wb') as f:
            pickle.dump(p1list,f)
        with open(tag+'/p3list2_'+str(sweeps)+'_'+str(steps)+'_'+str(date.today()),'wb') as f:
            pickle.dump(p3list,f)

        #reshape list to the correct shape of array for contour plotting
        av_infdarray=np.array(av_infdlist).reshape(len(p1list),len(p1list))
        var_infdarray=np.array(var_infd).reshape(len(p1list),len(p1list))
        with open(tag+'/av_infdarray2_'+str(sweeps)+'_'+str(steps)+'_'+str(date.today()),'wb') as f:
            pickle.dump(av_infdarray,f)
        with open(tag+'/var_infdarray2_'+str(sweeps)+'_'+str(steps)+'_'+str(date.today()),'wb') as f:
            pickle.dump(var_infdarray,f)

        #Plots
        X,Y = np.meshgrid(p1list,p3list)
        plt.contourf(X,Y,av_infdarray,30)
        plt.xlabel("P1 (S->I)")
        plt.ylabel("P3 (R->S)")
        plt.title("Phase diagram of " +'Average Infected Fraction')
        plt.colorbar();
        plt.savefig(tag+'/PhaseDiagram'+str(steps)+'_steps'+str(date.today()),format='pdf')
        plt.clf()

        X,Y = np.meshgrid(p1list,p3list)
        plt.contourf(X,Y,var_infdarray,30)
        plt.xlabel("P1 (S->I)")
        plt.ylabel("P3 (R->S)")
        plt.title("Phase diagram of " +'Variance')
        plt.colorbar();
        plt.savefig(tag+'/VarianceDiagram'+str(steps)+'_steps'+str(date.today()),format='pdf')
        plt.show()

################################################################################

    if mode=='slice':
        #constants
        p2=0.5
        p3=0.5
        p1=0.2
        steps=0.005
        sweeps=10100
        #lists
        p1list=[]
        av_infdlist=[]
        var_infd=[]

        #make a folder if one doesn't exist to keep outputs tidy
        tag=mode+str(steps)
        try:
            os.makedirs(mode+str(steps))
        except OSError:
            if not os.path.isdir(tag):
                raise+type

        #main loop for data analysis
        for a in np.arange(p1,0.5,steps):
            print('Changed p1 to: %.3f'%(a))
            infected=[]
            infected_sq=[]
            #reinstance class
            sirs=sirsc(a,p2,p3,lattice_size)
            p1list.append(float(sirs.p1))
            #begin sweeping
            for i in range(sweeps):
                for j in range(sirs.N):
                    sirs.Update()
                if i%10==0 and i>100:
                    inf=sirs.Calc_Infected()
                    infected.append(inf)
                    infected_sq.append(inf**2.0)
            av_infd=(np.average(infected))/sirs.N
            av_infdlist.append(av_infd)
            var_infd.append(sirs.Get_Variance(infected,infected_sq))


        #Pickling
        with open(tag+'/p1list2_'+str(sweeps)+'_'+str(steps)+'_'+str(date.today()),'wb') as f:
            pickle.dump(p1list,f)
        with open(tag+'/av_infdlist_'+str(sweeps)+'_'+str(steps)+'_'+str(date.today()),'wb') as f:
            pickle.dump(av_infdlist,f)
        with open(tag+'/var_infd_'+str(sweeps)+'_'+str(steps)+'_'+str(date.today()),'wb') as f:
            pickle.dump(var_infd,f)

        #plotting the variance against probability
        plt.plot(p1list, var_infd)
        plt.scatter(p1list, var_infd,c='k',marker='X')
        plt.xlabel('P1 (S=>I)')
        plt.ylabel('Variance (p3=0.5)')
        plt.savefig(tag+'/slice_in_'+str(steps)+'_steps'+str(date.today()),format='pdf')
        plt.show()

################################################################################
    if mode=='immune':
        '''
        Will eventually want a plot of x-axis = fraction of lattice immune
                               against y-axis = average infected.
        '''
        p1=0.5
        p2=0.5
        p3=0.5
        immunefrac=[]
        av_infdlist=[]
        yerrors=[]
        steps=0.005
        sweeps=10100
        for f in np.arange(0.0,1.0,steps):
            sirs=sirsc(p1,p2,p3,lattice_size)
            sirs.FracImmune(f)
            infected=[]
            print('Fraction of immune cells is:%.3f'%(f))
            immunefrac.append(f)
            #instance class
            for i in range(sweeps):
                for j in range(sirs.N):
                    sirs.Update()
                if i%10==0 and i>100:
                    inf=sirs.Calc_Infected()
                    infected.append(inf)
            av_infd=(np.average(infected))/sirs.N
            error=np.std(infected)/sirs.N
            stdm=error/np.sqrt(len(infected))
            yerrors.append(stdm)
            av_infdlist.append(av_infd)

        #make a folder if one doesn't exist to keep outputs tidy
        tag=mode+str(steps)
        try:
            os.makedirs(mode+str(steps))
        except OSError:
            if not os.path.isdir(tag):
                raise+type


        plt.plot(immunefrac,av_infdlist,'k')
        plt.errorbar(immunefrac,av_infdlist,yerr=yerrors,fmt='none',ecolor='r',capsize=4,capthick=1)
        plt.xlabel('$F_{im}$')
        plt.ylabel('Average Infected Fraction')
        plt.savefig(tag+'/AvInfectvsImunefrac'+str(lattice_size)+str(steps)+str(date.today()),format='pdf')
        plt.show()
















if __name__=='__main__':
    main()
