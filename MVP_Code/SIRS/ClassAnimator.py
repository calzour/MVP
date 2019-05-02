import time
import copy
from multiprocessing import Queue, Process
import matplotlib.animation as animation
import matplotlib.pyplot as plt


"""
Animator
"""

class Animator():
    def __init__(self, lattice_queue,frames):
        self.queue  = lattice_queue
        self.frames=frames
        data = lattice_queue.get()
        self.fig, self.ax_array = plt.subplots()
        self.mat = self.ax_array.matshow(data, vmin=-1, vmax=1)

    def update(self, i):
        if (self.queue.empty()):
            time.sleep(0.005)
        else:
            data = self.queue.get()
            self.mat.set_data(data)

    def animate(self):
        anim = animation.FuncAnimation(self.fig, self.update,interval=1,frames=self.frames)
        plt.show()
        return
