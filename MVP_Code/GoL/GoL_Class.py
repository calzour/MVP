import copy
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from scipy import ndimage
import warnings

warnings.simplefilter('ignore', np.RankWarning)

'''
R=Random, A=Absorbing, B=Blinker, G=Glider
'''


class GoL():
    def __init__(self, lattice_size=50,state='R'):
        self.lattice_size=lattice_size
        self.N=lattice_size**2.0
        self.state=state
        self.CoM_list=[]
        self.time_list=[]

        if state=='R':
            self.lattice=np.random.choice([0,1],size=(lattice_size,lattice_size))

        elif state=='A':
            '''
            Simple absoring state- might improve with beehive if time permits
            '''
            point1=int((np.random.uniform()*lattice_size))
            point2=int((np.random.uniform()*lattice_size))
            self.lattice=np.full((lattice_size,lattice_size),0)    #creates lattice full of -1 values
            for i in range(10):
                self.lattice[(point1+i*3)%self.lattice_size,(point2+i*3)%self.lattice_size]=1       #make a random site the opposite value to be absorbed


        elif state=='B':
            self.lattice=np.full((lattice_size,lattice_size),0)
            bpoint1=int((np.random.uniform()*lattice_size))
            bpoint2=int((np.random.uniform()*lattice_size))
            #Creates 3 points in a row
            for i in range(3):
                bpoint1+=1
                #condtion to catch boundary effects
                if bpoint1>= lattice_size:
                    bpoint1=0
                self.lattice[bpoint1,bpoint2]=1

        elif state=='G':
            self.lattice=np.full((lattice_size,lattice_size),0)
            #start glider at top left
            gpoint1=1
            gpoint2=1
            #Construct the glider
            self.lattice[gpoint1-1,gpoint2]=1
            self.lattice[gpoint1+1,gpoint2]=1
            self.lattice[gpoint1,gpoint2+1]=1
            self.lattice[gpoint1+1,gpoint2+1]=1
            self.lattice[gpoint1+1,gpoint2-1]=1



    def NN(self, i,j):
        '''
        Calculates the 8 nearest neighbours for a point (i, j)
        '''
        iright=(i+1)%self.lattice_size
        ileft=(i-1)%self.lattice_size
        jup=(j+1)%self.lattice_size
        jdown=(j-1)%self.lattice_size

        nn=[self.lattice[i,jup],self.lattice[i,jdown],self.lattice[iright,j],
        self.lattice[ileft,j],self.lattice[iright,jup],self.lattice[iright,jdown],
        self.lattice[ileft,jup],self.lattice[ileft,jdown]]
        return(nn)

    def Rules(self,i,j):
        nnList=self.NN(i,j)
        #Counter takes in the list of nearest neighbours so that number of alive neighbours can be counted
        cnt=Counter(nnList)
        if self.lattice[i,j]==1:
            #Rule for having two or three nearest neighbours meaning cell (i,j) stays alive
            if cnt[1]==2 or cnt[1]==3: return 1
            else: return 0
        #Rule to make dead cells alive
        elif self.lattice[i,j]== 0 and cnt[1]==3: return 1
        else: return 0

    def Update_Lattice(self):
        '''
        Creates a deep copy of the lattice where every element is 1
        Sweeps through the whole lattice once applying rules to each cell
        '''
        lattice_temp=copy.deepcopy(self.lattice)
        for i in range(self.lattice_size):
            for j in range(self.lattice_size):
                #Apply rules to point(i,j)
                cell= self.Rules(i,j)
                #update temporary lattice
                lattice_temp[i,j]=cell
        self.lattice=lattice_temp

    def Get_CoM(self,time):
        x,y=ndimage.measurements.center_of_mass(self.lattice)
        #Ignore boundary conditions
        if 3< x < (self.lattice_size-5) and 3< y < (self.lattice_size-5):
            time+=1
            if time%200==0:
                time=0
            self.time_list.append(time)
            CoM=np.sqrt(x**2 + y**2)
            self.CoM_list.append(CoM)
        return (self.time_list,self.CoM_list,time)

    def Get_Velocity(self,time_list,CoM_list):
        #fit a power of 1 polynomial (i.e. a straight line) to get gradient
        coeffs=np.polyfit(time_list,CoM_list,1)
        return coeffs[0]
