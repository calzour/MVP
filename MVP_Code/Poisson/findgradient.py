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
from scipy.optimize import curve_fit


r_list = pickle.load(open('r_list22019-04-05_20','rb'))
pot_list=pickle.load(open('pot_list22019-04-05_20','rb'))

'''
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
'''

with open('pot_data.dat','w') as writer:
    for r,pot in zip(r_list, pot_list):
        writer.write('{:.5f} {:.5f}\n'.format(r,pot))


r_list=np.array(r_list)
r_list = r_list[np.logical_not(np.isnan(r_list))]
r_list = r_list[np.logical_not(np.isinf(r_list))]

pot_list=np.array(pot_list)
pot_list = pot_list[np.logical_not(np.isnan(pot_list))]
pot_list = pot_list[np.logical_not(np.isinf(pot_list))]
for i in r_list:
    if i == 'NaN':
        i=0.0
    elif i =='inf':
        i=0.0

for i in pot_list:
    if i == 'NaN':
        i=0.0
    elif i =='inf':
        i=0.0

print(len(pot_list))
print(len(r_list))
pot_list=pot_list[:10]
r_list=r_list[:10]


def f(x, A, B): # this is your 'straight line' y=f(x)
    return A*x + B


m,c = curve_fit(f,r_list,pot_list)[0]
print(m)

plt.scatter(r_list,pot_list)
plt.show()
