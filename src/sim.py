"""
This file is licensed by the GNU GPLv3 License
Copyright © 2024 RandomKiddo
"""

import matplotlib.pyplot as plt
import warnings
import argparse
import os
import numpy as np
import csv
import ast

from cgen import *
from collections import defaultdict

_TYPES = Literal['mks', 'cgs', 'au']

SCALE_FACTORS = {'au': 100_000}


class Simulation:
    def __init__(self, c_system: System, system: _TYPES = 'mks', cutoff: float = 0.0):
        """
        Initializes a new simulation.
        :param c_system: The celestial system to use.
        :param system: The system of units to use. Defaults to 'mks'.
        :param cutoff: The cutoff distance to not include in simulation.
        """
        self.c_system = c_system
        self.system = system

        self.fig, self.ax = plt.subplots(1, 1, subplot_kw={'projection': 'polar'})
        self.fig.tight_layout()
        self.ax.axis('off')

        if self.c_system is not None:
            self.stars = []
            for _ in self.c_system:
                if _.star:
                    self.stars.append(_)

            self.initialize_conditions()

            self.paths = defaultdict(lambda: [])
            self.ensure_naming_conventions()

            self.dthetas = defaultdict(lambda: [])

            if cutoff > 0.0:
                copy = []
                for _ in self.c_system:
                    if _.a <= cutoff:
                        copy.append(_)
                self.c_system = System(copy)

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
                         markersize=_.r / Constants.solar_radius(self.system), color=f'#{_.color}'  # todo change marker size
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
                v_perihelion = math.sqrt(Constants.G(self.system) * self.eff_star().m * (1 / _.a) * ((1 + _.e) / (1 - _.e)))
                if len(self.dthetas[_.label]) != 0:
                    avg = sum(self.dthetas[_.label]) / len(self.dthetas[_.label])
                else:
                    avg = 100  # arbitrary large value
                if _.data.velocity.theta_comp > v_perihelion:
                    _.data.velocity.theta_comp = v_perihelion
                dtheta = _.data.velocity.theta_comp * dt
                if dtheta > avg*2:  # todo fix, broken for Mercury and alike
                    dtheta = avg
                self.dthetas[_.label].append(dtheta)
                _.data.position.theta_comp += dtheta
                _.data.position.r_comp = _.pos()
                _.data.velocity.r_comp = _.v_r(focus=self.eff_star(), system=self.system)
                _.data.velocity.theta_comp = _.v_theta(focus=self.eff_star(), system=self.system)

    def sim(self, dt: float = 0.01, sp: str = None, pt: float = 0.0001) -> None:
        """
        Runs the simulation. Use a keyboard stroke in the terminal to interrupt.
        :param dt: The differential time step to use. Defaults to 0.01.
        :param sp: The string save path to save the .npz file. Defaults to None.
        :param pt: The plot pause time to use, should be proportional to dt.
        """
        if dt >= 0.01:
            warnings.warn(f'Differential time step dt relatively large at {dt}. Using large dt could yield lossy simulation results')
        try:
            while True:
                self.step(dt=dt)
                self.draw()
                plt.pause(pt)
                self.ax.clear()
                self.ax.axis('off')
        except KeyboardInterrupt:
            if sp is not None:
                with open(sp, 'w', newline='') as f:
                    w = csv.writer(f)
                    w.writerow(self.paths.keys())
                    for i in range(len(self.paths[list(self.paths.keys())[0]])):
                        row = []
                        for ii in list(self.paths.keys()):
                            row.append(self.paths[ii][i])  # todo save metadata
                        w.writerow(row)
                f.close()

    @staticmethod
    def sim_from_csv(fp: str, delimiter: str = ',', system: _TYPES = 'mks', metadata: Any = None) -> None:
        """
        Simulates from a CSV file.
        :param fp: The filepath to the CSV.
        :param delimiter: The CSV delimiter. Defaults to ','.
        :param system: The system of units to use. Defaults to 'mks'.
        :param metadata: The metadata information of the system.
        """
        if not os.path.exists(fp):
            raise FileNotFoundError(f'Could not find .csv file {fp}')
        with open(fp, newline='') as f:
            r = csv.reader(f, delimiter=delimiter)
            faulty_sim = Simulation(c_system=None, system=system)
            paths = defaultdict(lambda: [])
            labels_in_order = []
            for row in r:
                for _ in row:
                    paths[_] = []
                    labels_in_order.append(_)
                break
            try:
                while True:
                    row = next(r)
                    for _ in range(len(row)):
                        lit = ast.literal_eval(row[_])
                        paths[labels_in_order[_]].append(lit)
                        z = list(zip(*paths[labels_in_order[_]]))
                        faulty_sim.ax.plot(z[1], z[0], marker='o')
                    plt.pause(0.0001)
                    faulty_sim.ax.axis('off')
            except KeyboardInterrupt:
                f.close()
            except StopIteration:
                f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='SolarSim', description='A python-based celestial body system simulator')
    parser.add_argument('--dt', type=float, action='store',
                        default=0.01, help='differential time step to use in the simulation')
    parser.add_argument('-u', '--units', type=str, action='store', default='mks', help='system of units to use for the simulation')
    parser.add_argument('-s', '--solar', action='store_true', default=False, help='use the solar system for the simulation')
    parser.add_argument('--fp', type=str, action='store', default=None, help='the filepath to a json file representing the system')
    parser.add_argument('--sp', type=str, action='store', default=None, help='the path to save the simulation paths as a csv file')
    parser.add_argument('-cut', '--cutoff', type=float, action='store', default=0.0, help='the zoom factor to zoom-in by')  # todo
    parser.add_argument('-g', '--gif', type=str, action='store', default=None, help='store the simulation as a gif with the provided filepath')
    parser.add_argument('--pt', type=float, action='store', default=0.0001, help='plot pause timing to use in the simulation')
    parser.add_argument('-v', '--verbose', action='store_true', default=False, help='log more information to the console')

    args = parser.parse_args()

    if args.solar:
        sim = Simulation(c_system=System(), system=args.units, cutoff=args.cutoff)
    elif args.fp is not None and os.path.exists(args.fp):
        sim = Simulation(c_system=System.read_from_json(args.fp), system=args.units, cutoff=args.cutoff)
    elif args.fp is None:
        warnings.warn('Solar system argument not specified, yet no filepath was provided. Defaulting to solar system.')
        sim = Simulation(c_system=System(), system=args.units, cutoff=args.cutoff)
    else:  # file doesn't exist
        warnings.warn(f'Provided filepath {args.fp} could not be located. Defaulting to solar system.')
        sim = Simulation(c_system=System(), system=args.units, cutoff=args.cutoff)

    if args.gif is not None:
        warnings.warn('Gif functionality not currently integrated-- will be added in next update.')
    if args.sp is not None:
        warnings.warn('Save-path functionality working, but reading from csv does not save all metadata.')

    # todo verbose and restart line

    sim.sim(dt=args.dt, sp=args.sp, pt=args.pt)
