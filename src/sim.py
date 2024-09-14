import numpy as np
import matplotlib.pyplot as plt
import scipy

from cgen import *

_TYPES = Literal['mks', 'cgs', 'au']

class Simulation:
    def __init__(self, c_system: System, system: _TYPES = 'mks'):
        self.c_system = c_system
        self.system = system

        self.fig, self.ax = plt.subplots(1, 1, subplot_kw={'projection': 'polar'})
        self.fig.tight_layout()
        self.ax.axis('off')

    def draw(self):
        for _ in self.c_system:
            pos = (_.data.position.r_comp, _.data.position.theta_comp)
            self.ax.plot(pos[1], pos[0], marker='o',
                         markersize=_.r/Constants.solar_radius(self.system), color=f'#{_.color}'
            )


if __name__ == '__main__':
    sim = Simulation(c_system=System(), system='au')
    sim.draw()
    plt.show()

