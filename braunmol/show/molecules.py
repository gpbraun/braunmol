"""
molecules.py

Gabriel Braun, 2026
"""

from pymol import cmd


@cmd.extend
def show_vdw(
    mol_obj=None,
    _self=cmd,
):
    """
    Representação: Esferas de Van der Waals.

    >>> PyMOL> show_vdw [ mol_obj ]
    """
    mol_obj = mol_obj or _self.get_object_list()[0]

    _self.show("spheres", mol_obj)

    _self.set("sphere_scale", 0.8, mol_obj)
    _self.set("sphere_scale", 0.7, f"{mol_obj} and elem H")


@cmd.extend
def show_bas(
    mol_obj=None,
    _self=cmd,
):
    """
    Representação: Ball and Stick.

    >>> PyMOL> show_bas [ mol_obj ]
    """
    mol_obj = mol_obj or _self.get_object_list()[0]

    _self.show("sticks", mol_obj)
    _self.show("spheres", mol_obj)

    # configurações: sticks
    _self.set("stick_radius", 0.10, mol_obj)
    _self.set("stick_h_scale", 0.80, mol_obj)
    _self.set("stick_color", "white", mol_obj)

    # configurações: balls
    _self.set("sphere_scale", 0.20, mol_obj)
    _self.set("sphere_scale", 0.17, f"{mol_obj} and elem H")


@cmd.extend
def show_ion(
    mol_obj=None,
    spheres=True,
    sticks=False,
    _self=cmd,
):
    """
    Representa: rede de íons.

    >>> PyMOL> show_ionic [ mol_obj [, spheres [, sticks ]]]
    """
    mol_obj = mol_obj or _self.get_object_list()[0]

    _self.set("sphere_scale", 0.20, mol_obj)

    _self.set("stick_radius", 0.075, mol_obj)
    _self.set("stick_h_scale", 1.00, mol_obj)
    _self.set("stick_color", -1, mol_obj)

    if int(spheres):
        _self.show("spheres", mol_obj)
    if int(sticks):
        _self.show("sticks", mol_obj)


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

    _self.set("sphere_scale", 0.8, vdw_obj)
    _self.set("sphere_scale", 0.7, f"{vdw_obj} and elem H")

    _self.hide("nonbonded", vdw_obj)
    _self.hide("lines", vdw_obj)
    _self.hide("sticks", vdw_obj)
    _self.set("sphere_transparency", 0.5, vdw_obj)
