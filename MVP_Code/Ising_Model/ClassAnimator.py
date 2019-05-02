import time
import copy
from multiprocessing import Queue, Process
import matplotlib.animation as animation
import matplotlib.pyplot as plt


"""
Animator
"""

class Animator():
    def __init__(self, lattice_queue):
        self.queue  = lattice_queue
        data = lattice_queue.get()
        self.fig, self.ax_array = plt.subplots()
        cont = self.ax_array.contourf(data, vmin=-1, vmax=1)

    def update(self, i):
        if (self.queue.empty()):
            time.sleep(0.1)
        else:
            self.ax_array.cla()
            data = self.queue.get()
            cont = self.ax_array.contourf(data, vmin=-1, vmax=1)
            return cont

    def animate(self):
        anim = animation.FuncAnimation(self.fig, self.update,interval=1)
        plt.show()
