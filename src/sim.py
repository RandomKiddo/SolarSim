"""
This file is licensed by the GNU GPLv3 License
Copyright © 2024 RandomKiddo
"""

import matplotlib.pyplot as plt
import warnings
import argparse
import os

from cgen import *
from collections import defaultdict

_TYPES = Literal['mks', 'cgs', 'au']

SCALE_FACTORS = {'au': 100_000}


class Simulation:
    def __init__(self, c_system: System, system: _TYPES = 'mks'):
        """
        Initializes a new simulation.
        :param c_system: The celestial system to use.
        :param system: The system of units to use. Defaults to 'mks'.
        """
        self.c_system = c_system
        self.system = system

        self.fig, self.ax = plt.subplots(1, 1, subplot_kw={'projection': 'polar'})
        self.fig.tight_layout()
        self.ax.axis('off')

        self.stars = []
        for _ in self.c_system:
            if _.star:
                self.stars.append(_)

        self.initialize_conditions()

        self.paths = defaultdict(lambda: [])
        self.ensure_naming_conventions()

    def ensure_naming_conventions(self) -> None:
        """
        Ensures naming conventions for storing paths of the celestial bodies.
        """
        count = 0
        for _ in self.c_system:
            if not _.star and _.label == '':
                _.label = f'CelestialBody_{count}'
                count += 1

    def draw(self) -> None:
        """
        Draws the current position of all the bodies of the system.
        """
        for _ in self.c_system:
            pos = (SCALE_FACTORS[self.system] * _.data.position.r_comp, _.data.position.theta_comp)
            self.paths[_.label].append(pos)
            z = list(zip(*self.paths[_.label]))
            self.ax.plot(z[1], z[0], marker='o',
                         markersize=_.r / Constants.solar_radius(self.system), color=f'#{_.color}'
                         )

    def eff_star(self) -> CelestialBody:
        """
        Calculates the effective star for a binary (or greater) star system. This uses the "loose addition" of
        celestial bodies, defined in the CelestialBody class.
        :return: The effective star as a CelestialBody instance.
        """
        body = self.stars[0]
        for _ in self.stars[1:]:
            body += _
        return body

    def initialize_conditions(self) -> None:
        """
        Initializes the start conditions of the simulation, including initial position and velocity polar vectors.
        """
        for _ in self.c_system:
            if not _.star:
                _.data.position = PolarVector(r_comp=_.a, theta_comp=0.0)
                star = self.eff_star()
                _.data.velocity = PolarVector(r_comp=0.0,
                                              theta_comp=math.sqrt(
                                                  (1 + _.e) * (_.alpha(focus=star, system=self.system)) / _.a))

    def step(self, dt: float = 0.01) -> None:
        """
        Calculates the next step of the simulation, calculating new positions and velocities of the bodies of the system.
        :param dt: The time step to use. Defaults to 0.01.
        """
        for _ in self.c_system:
            if not math.isclose(_.a, 0.0):
                dtheta = _.data.velocity.theta_comp * dt
                _.data.position.theta_comp += dtheta
                _.data.position.r_comp = _.pos()
                _.data.velocity.r_comp = _.v_r(focus=self.eff_star(), system=self.system)
                _.data.velocity.theta_comp = _.v_theta(focus=self.eff_star(), system=self.system)

    def sim(self, dt: float = 0.01) -> None:
        """
        Runs the simulation.
        """
        if dt >= 0.01:
            warnings.warn(f'Differential time step dt relatively large at {dt}. Using large dt could yield lossy simulation results')
        while True:
            self.step(dt=dt)
            self.draw()
            plt.pause(0.0001)
            self.ax.clear()
            self.ax.axis('off')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='SolarSim', description='A python-based celestial body system simulator')
    parser.add_argument('--dt', type=float, action='store',
                        default=0.01, help='differential time step to use in the simulation')
    parser.add_argument('-u', '--units', type=str, action='store', default='mks', help='system of units to use for the simulation')
    parser.add_argument('-s', '--solar', action='store_true', default=False, help='use the solar system for the simulation')
    parser.add_argument('-p', '--pause', type=float, action='store', default=0.0001,
                        help='amount of time to pause the simulation before updating')
    parser.add_argument('--fp', type=str, action='store', default=None, help='the filepath to a json file representing the system')

    args = parser.parse_args()

    if args.solar:
        sim = Simulation(c_system=System(), system=args.units)
    if args.fp is not None and os.path.exists(args.fp):
        sim = Simulation(c_system=System.read_from_json(args.fp), system=args.units)
    elif args.fp is None:
        warnings.warn('Solar system argument not specified, yet no filepath was provided. Defaulting to solar system.')
        sim = Simulation(c_system=System(), system=args.units)
    else:  # file doesn't exist
        warnings.warn(f'Provided filepath {args.fp} could not be located. Defaulting to solar system.')
        sim = Simulation(c_system=System(), system=args.units)

    sim.sim(dt=args.dt)
