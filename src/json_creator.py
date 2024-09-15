"""
This file is licensed by the GNU GPLv3 License
Copyright © 2024 RandomKiddo
"""

import json

from cgen import *


def create(fp: str, c_system: System = None, list_system: List[CelestialBody] = None) -> None:
    if c_system is None and list_system is None:
        raise ValueError('One of c_system or list_system required. Both cannot be None.')

    system = []
    if c_system is not None:
        system = c_system.bodies
    else:
        system = list_system

    data = {}
    count = 0
    for _ in system:
        label = _.label if _.label != '' else f'Celestial Body {count}'
        info = {
            'm': _.m, 'r': _.r, 'a': _.a, 'color': _.color, 'e': _.e, 'star': _.star
        }
        data[label] = info

    json_object = json.dumps(data, indent=4)

    with open(fp, 'w') as f:
        f.write(json_object)
    f.close()
