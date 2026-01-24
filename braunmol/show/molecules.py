"""
molecules.py

Gabriel Braun, 2026
"""

from pymol import cmd


@cmd.extend
def show_vdw(
    mol_obj=None,
    spheres="large",
    _self=cmd,
):
    """
    Representação: Esferas de Van der Waals.

    >>> PyMOL> show_vdw [ mol_obj ]
    """
    mol_obj = mol_obj or _self.get_object_list()[0]

    cmd.show("spheres", mol_obj)

    cmd.set("sphere_scale", 0.8, mol_obj)
    cmd.set("sphere_scale", 0.6, f"{mol_obj} and elem H")


@cmd.extend
def show_bas(
    mol_obj=None,
    spheres="large",
    _self=cmd,
):
    """
    Representação: Ball and Stick.

    >>> PyMOL> show_bas [ mol_obj [, shperes ] ]
    """
    mol_obj = mol_obj or _self.get_object_list()[0]

    cmd.show("sticks", mol_obj)
    cmd.show("spheres", mol_obj)

    # configurações: sticks
    cmd.set("stick_radius", 0.10, mol_obj)
    cmd.set("stick_h_scale", 0.80, mol_obj)
    cmd.set("stick_color", "white", mol_obj)

    # configurações: balls
    if spheres == "large":
        cmd.set("sphere_scale", 0.20, mol_obj)
        cmd.set("sphere_scale", 0.17, f"{mol_obj} and elem H")
    if spheres == "small":
        cmd.set("sphere_scale", 0.17, mol_obj)
        cmd.set("sphere_scale", 0.14, f"{mol_obj} and elem H")


@cmd.extend
def add_vdw(
    mol_obj=None,
    vdw_obj=None,
    _self=cmd,
):
    """
    Adição: Esferas de Van der Waals

    >>> PyMOL> ball_and_stick <mol_obj>
    """
    mol_obj = mol_obj or _self.get_object_list()[0]
    vdw_obj = vdw_obj or f"{mol_obj}_vdw"

    _self.create(vdw_obj, mol_obj)

    _self.set("sphere_scale", 1.0, f"{vdw_obj} and elem H")
    _self.set("sphere_scale", 0.8, vdw_obj)

    _self.hide("nonbonded", vdw_obj)
    _self.hide("lines", vdw_obj)
    _self.hide("sticks", vdw_obj)
    _self.set("sphere_transparency", 0.5, vdw_obj)
