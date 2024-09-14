import warnings
import math
import json

from typing import *
from typing_extensions import *
from collections import Iterable

_TYPES = Literal['mks', 'cgs', 'au']


class Constants:
    @staticmethod
    def G(system: _TYPES = 'mks') -> float:
        """
        Returns the value of Newton's gravitational constant, G.
        :param system: The system of units to return the value of (either 'mks', 'cgs', or 'au'). Defaults to 'mks'.
        :return: The float value corresponding to the given system of units.
        """
        options = get_args(_TYPES)
        if system not in options:
            warnings.warn(f'{system} not a recognized for of unit system, defaulting to MKS for G')
        if system == 'mks' or system not in options:
            return 6.6743015e-11
        elif system == 'cgs':
            return 6.6743015e-8
        else:
            return 4*(math.pi**2)

    @staticmethod
    def solar_mass(system: _TYPES = 'mks') -> float:
        """
        Returns the value of solar mass
        :param system: The system of units to return the value of (either 'mks', 'cgs', or 'au'). Defaults to 'mks'.
        :return: The float value corresponding to the given system of units.
        """
        options = get_args(_TYPES)
        if system not in options:
            warnings.warn(f'{system} not a recognized for of unit system, defaulting to MKS for solar_mass')
        if system == 'mks' or system not in options:
            return 1.9885e30
        elif system == 'cgs':
            return 1.9885e33
        else:
            return 1.0

    @staticmethod
    def solar_radius(system: _TYPES = 'mks') -> float:
        """
        Returns the value of solar radius
        :param system: The system of units to return the value of (either 'mks', 'cgs', or 'au'). Defaults to 'mks'.
        :return: The float value corresponding to the given system of units.
        """
        options = get_args(_TYPES)
        if system not in options:
            warnings.warn(f'{system} not a recognized for of unit system, defaulting to MKS for solar_radius')
        if system == 'mks' or system not in options:
            return 6.957e8
        elif system == 'cgs':
            return 6.957e10
        else:
            return 0.0046524726

    @staticmethod
    def earth_mass(system: _TYPES = 'mks') -> float:
        """
        Returns the value of earth mass
        :param system: The system of units to return the value of (either 'mks', 'cgs', or 'au'). Defaults to 'mks'.
        :return: The float value corresponding to the given system of units.
        """
        options = get_args(_TYPES)
        if system not in options:
            warnings.warn(f'{system} not a recognized for of unit system, defaulting to MKS for earth__mass')
        if system == 'mks' or system not in options:
            return 5.972168e24
        elif system == 'cgs':
            return 5.972168e27
        else:
            return 3.0034146856629e-06

    @staticmethod
    def earth_radius(system: _TYPES = 'mks') -> float:
        """
        Returns the value of earth radius
        :param system: The system of units to return the value of (either 'mks', 'cgs', or 'au'). Defaults to 'mks'.
        :return: The float value corresponding to the given system of units.
        """
        options = get_args(_TYPES)
        if system not in options:
            warnings.warn(f'{system} not a recognized for of unit system, defaulting to MKS for earth_radius')
        if system == 'mks' or system not in options:
            return 6378.137
        elif system == 'cgs':
            return 6378137
        else:
            return 0.00004263521

    @staticmethod
    def au(system: _TYPES = 'mks') -> float:
        """
        Returns the value of astronomical unit
        :param system: The system of units to return the value of (either 'mks', 'cgs', or 'au'). Defaults to 'mks'.
        :return: The float value corresponding to the given system of units.
        """
        options = get_args(_TYPES)
        if system not in options:
            warnings.warn(f'{system} not a recognized for of unit system, defaulting to MKS for au')
        if system == 'mks' or system not in options:
            return 1.495979e11
        elif system == 'cgs':
            return 1.495979e13
        else:
            return 1

class PolarVector:
    def __init__(self, r_comp: float = 0.0, theta_comp: float = 0.0) -> None:
        """
        Initializes a new polar vector.
        :param r_comp: The r-component, defaults to 0.0.
        :param theta_comp: The theta-component, defaults to 0.0.
        """
        self.r_comp = r_comp
        self.theta_comp = theta_comp

    def dist(self, other: Self) -> float:
        """
        Returns the distance between two polar vectors.
        :param other: The other polar vector.
        :return: The distance as a float.
        """
        return math.sqrt(self.r_comp**2 + other.r_comp**2 - 2*self.r_comp*other.r_comp*math.cos(other.theta_comp-self.theta_comp))

    def mag_cross(self, u: Self) -> float:
        return (self.r_comp*u.theta_comp)-(self.theta_comp*u.r_comp)

class DataPair:
    def __init__(self, position: PolarVector, velocity: PolarVector):
        """
        Initializes a new data pair of position and velocity of an object in polar coordinates.
        :param position: The position polar vector.
        :param velocity: The velocity polar vector.
        """
        self.position = position
        self.velocity = velocity


class CelestialBody:
    def __init__(self, m: float, r: float, a: float = 0.0, color: str = '000000',
                 data: DataPair = None, e: float = 0.0, star: bool = False, label: str = '') -> None:
        """
        Initializes a new celestial body.
        :param m: The mass of the body.
        :param r: The radius of the body.
        :param a: The distance of the body from its system center. Defaults to 0.0.
        :param color: The color of the body as a hex code (without #). Defaults to black (#000000).
        :param data: The data pairing of the body. Defaults to None.
        :param e: The eccentricity value e of the orbit. Defaults to 0.0
        :param star: Boolean value representing if this body is a star. Defaults to False.
        :param label: The name or label of the celestial body. Defaults to an empty str.
        """
        self.m = m
        self.r = r
        self.a = a
        self.color = color
        self.e = e
        self.star = star
        self.label = label

        if data is not None:
            self.data = data
        else:
            self.data = DataPair(PolarVector(a), PolarVector())

    def force(self, cb: Self, a: float, system: _TYPES = 'mks') -> float:
        """
        Calculates the Newtonian force between this celestial body and another.
        :param cb: The second celestial body.
        :param a: The distance between the bodies. The radius of both bodies should be already taken into account.
        :param system: The system of units to use for G (either 'mks', 'cgs', or 'au'). Defaults to 'mks'.
        :return: The float value of the force between the bodies.
        """
        return Constants.G(system) * self.m * cb.m / (a**2)

    def __eq__(self, other: Self) -> bool:
        """
        Checks if two celestial bodies are equal by using `math.isclose`. Does not check for equality of color, only mass, radius, and distance.
        :param other: The other CelestialBody to compare to.
        :return: True if equal, False otherwise.
        """
        return math.isclose(self.m, other.m) and math.isclose(self.r, other.r) and math.isclose(self.a, other.a)

    def __neq__(self, other: Self) -> bool:
        """
        Checks if two celestial bodies are not equal by using `math.isclose`. Does not check for equality of color, only mass, radius, and distance.
        :param other: The other CelestialBody to compare to.
        :return: True if not equal, False otherwise.
        :param other:
        :return:
        """
        return not __eq__(self, other)

    def pos(self) -> float:
        """
        Returns the current position of the celestial body based on Kepler's elliptical orbit formulae.
        :return: The position float value (radially).
        """
        return (self.a*(1-self.e**2))/(1+self.e*math.cos(self.data.position.theta_comp))

    def alpha(self, focus: Self, system: _TYPES = 'mks') -> float:
        """
        Returns the alpha value of the body with respect to a focus, G(m+M).
        :param focus: The CelestialBody of focus.
        :param system: The system of units to use. Defaults to 'mks'.
        :return: The float alpha value.
        """
        return Constants.G(system) * (self.m + focus.m)

    def mag_H(self) -> float:
        """
        Returns the magnitude of the H vector, defined as r x r dot (r cross r dot).
        :return: The magnitude as a float.
        """
        return self.data.position.mag_cross(self.data.velocity)

    def v_r(self, focus: Self, system: _TYPES = 'mks') -> float:
        """
        Returns the radial velocity at the current position.
        :param focus: The CelestialBody focus of the system.
        :param system: The system of units to use. Defaults to 'mks'
        :return: The radial velocity as a float.
        """
        return (self.e*math.sin(self.data.position.theta_comp))*math.sqrt((self.alpha(focus, system)**2)/(self.mag_H()**2))

    def v_theta(self, focus: Self, system: _TYPES = 'mks') -> float:
        """
        Returns the tangential (theta) velocity at the current position.
        :param focus: The CelestialBody focus of the system.
        :param system: The system of units to use. Defaults to 'mks'/
        :return: The theta velocity as a float.
        """
        return (1 + self.e * math.cos(self.data.position.theta_comp)) * math.sqrt(
            (self.alpha(focus, system) ** 2) / (self.mag_H() ** 2))

    def __add__(self, other: Self) -> Self:
        """
        Defines how to (loosely) add two celestial bodies. This method should only be used for effective bodies.
        This defines a new celestial body of combined mass, no radius, and averaged semi-major axis.
        :param other: The other celestial body to add.
        :return: The new "added" celestial body.
        """
        a = (self.a + other.a) / 2
        return CelestialBody(m=self.m+other.m, r=0.0, a=a, star=self.star)

    def __iadd__(self, other: Self) -> None:
        """
        Defines how to (loosely) add two celestial bodies. This method should only be used for effective bodies.
        This defines a new celestial body of combined mass, no radius, and averaged semi-major axis. Sets the current
        celestial body to the "added" one.
        :param other: The other celestial body to add.
        """
        self.a = (self.a + other.a) / 2
        self.m += other.m
        self.r = 0.0

class System(Iterable):
    def __init__(self, bodies: List[CelestialBody] = None) -> None:
        """
        Initializes a new system of stars and planets.
        :param bodies: A list of celestial bodies of the system.
        """
        if bodies is None:
            bodies = System.solar_system()

        self.bodies = bodies

        self.max_distance = max([bodies[_].a for _ in range(len(bodies))])
        self.max_mass = max([bodies[_].m for _ in range(len(bodies))])

    def __getitem__(self, idx: int) -> CelestialBody:
        """
        Gets the celestial body instance at the given index.
        :param idx: The index.
        :return: The CelestialBody instance.
        """
        return self.bodies[idx]

    def __iter__(self) -> Iterator:
        """
        Returns the iterator for the system.
        :return: The iterator instance.
        """
        self.counter = 0
        return self

    def __next__(self) -> CelestialBody:
        """
        Gets the next item of the iterator.
        :return: The next CelestialBody item.
        """
        self.counter += 1
        if self.counter <= len(self.bodies):
            return self.bodies[self.counter-1]
        else:
            self.counter = 0
            raise StopIteration

    @staticmethod
    def solar_system() -> List[CelestialBody]:
        """
        Creates a list of celestial bodies representing our solar system.
        :return: The list of bodies.
        """
        sys = []
        with open('solar_system.json', 'r') as f:
            data = json.load(f)
            for _ in data:
                body = data[_]
                sys.append(CelestialBody(m=body['m'], r=body['r'], a=body['a'], color=body['color'],
                                         e=body['e'], star=body['star'], label=str(_)))
        f.close()
        return sys

