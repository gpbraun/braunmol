"""
core.py

Gabriel Braun, 2026
"""

__all__ = [
    "COLORS",
    "ELEMENTS",
]


import importlib
import json

import matplotlib as mpl
from pymol import cmd

PYMOL_RESOURCES_PATH = importlib.resources.files("braunmol._data")

TWCOLORS_PATH = PYMOL_RESOURCES_PATH.joinpath("colors.json")
ELEMENTS_PATH = PYMOL_RESOURCES_PATH.joinpath("elements.json")

DEFAULT_CONFIG = {
    "fog": 0,
    "orthoscopic": 0,
    "transparency": 0.5,
    "antialias": 3,
    "ambient": 0.5,
    "spec_count": 5,
    "shininess": 50,
    "specular": 1,
    "reflect": 0.1,
    "ray_trace_mode": 1,
    "ray_trace_gain": 4,
    "ray_opaque_background": 0,
}

for command, value in DEFAULT_CONFIG.items():
    cmd.set(command, value)


COLORS = {}

for color_name_base, color_value in json.loads(TWCOLORS_PATH.read_text()).items():
    for shade, hex_code in color_value.items():
        pymol_name = f"tw_{color_name_base}_{shade}"
        COLORS[pymol_name] = mpl.colors.to_rgb(hex_code)
        cmd.set_color(pymol_name, COLORS[pymol_name])

cmd.bg_color("white")


ELEMENTS = json.loads(ELEMENTS_PATH.read_text())

for elem_symbol, elem_props in ELEMENTS.items():
    cmd.alter(f"elem {elem_symbol}", f"vdw={elem_props['vdw']:.2f}")
    cmd.set_color(elem_props["name"], COLORS[elem_props["color"]])
