"""
surfaces.py

Gabriel Braun, 2026
"""

import matplotlib as mpl
import numpy as np
from pymol import cmd

CMAP_COLOR_NUM = 32
"""
Número de samples do colormap.
"""


@cmd.extend
def mep_surface(
    den_cube,
    mep_cube,
    obj_name=None,
    den_level=0.02,
    mep_range=[-0.10, 0.40],
    cmap="Spectral",
    _self=cmd,
):
    """
    Superfície: Potencial Eletrostático Molecular (MEP).

    >>> PyMOL> mep_surface obj_name, den_cube, mep_cube
    """
    obj_name = obj_name or _self.get_object_list()[0]

    _self.isosurface(f"{obj_name}_mep", den_cube, den_level)

    colors = mpl.colormaps[cmap](np.linspace(0, 1, CMAP_COLOR_NUM))
    for i, color in enumerate(colors):
        cmd.set_color(f"__{cmap}_{i}", mpl.colors.to_rgb(color))

    _self.ramp_new(
        f"__{obj_name}_ramp",
        mep_cube,
        color=[f"__{cmap}_{i}" for i in range(CMAP_COLOR_NUM)],
        range=mep_range,
    )
    _self.disable(f"__{obj_name}_ramp")

    _self.color(f"__{obj_name}_ramp", f"{obj_name}_mep")
    _self.set("transparency", 0.5, f"{obj_name}_mep")
