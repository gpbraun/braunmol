"""
cell.py

Gabriel Braun, 2026
"""

import tempfile
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.local_env import (
    CrystalNN,
    MinimumVIRENN,
    ValenceIonicRadiusEvaluator,
)
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import Molecule, Structure
from pymol import cgo, cmd


@dataclass(slots=True)
class CellInfo:
    lattice: object
    repeats: np.ndarray


class CellRegistry:
    __slots__ = ("_data",)

    def __init__(self):
        self._data: dict[str, CellInfo] = {}

    def new(self, mol_obj: str, lattice, repeats: np.ndarray) -> None:
        self._data[mol_obj] = CellInfo(
            lattice=lattice.copy(), repeats=np.asarray(repeats, dtype=int)
        )

    def get(self, mol_obj: str) -> CellInfo:
        return self._data[mol_obj]


def _cell_registry(_self) -> CellRegistry:
    """
    Singleton associado ao objeto cmd do PyMOL.
    """
    reg = getattr(_self, "_cell_registry", None)
    if reg is None:
        reg = CellRegistry()
        setattr(_self, "_cell_registry", reg)
    return reg


def get_pb_frac_coords(frac_coords):
    """
    Retorna: coords fracionais no range [0,1] a partir de coords no range [0,1).
    """
    yield frac_coords

    zero_indices = np.where(np.isclose(frac_coords, 0))[0]

    for r in range(1, len(zero_indices) + 1):
        for indices_to_flip in combinations(zero_indices, r):
            pb_frac_coords = frac_coords.copy()
            pb_frac_coords[list(indices_to_flip)] = 1
            yield pb_frac_coords


def structure_to_molecule(structure, repeats):
    """
    Retorna: "molécula" contendo os átomos de uma supercélula.
    """
    molecule = Molecule([], [])

    for site in structure.make_supercell(repeats).sites:
        for frac_coords in get_pb_frac_coords(site.frac_coords):
            coords = structure.lattice.get_cartesian_coords(frac_coords)
            if hasattr(molecule, "specie"):
                molecule.append(site.specie, coords)
            elif hasattr(molecule, "species"):
                molecule.append(site.species, coords)

    return molecule


@cmd.extend
def load_cell(
    file_path,
    mol_obj=None,
    reps="1-1-1",
    _self=cmd,
    **kwargs,
):
    """
    Carrega: célula unitária (na forma de molécula) a partir de um arquivo '.cif'.

    >>> PyMOL> load_cell file_path [, mol_obj [, reps ]]
    """
    mol_obj = mol_obj or Path(file_path).stem

    repeats = np.fromstring(reps, sep="-", count=3, dtype=int)
    structure = Structure.from_file(file_path, check_cif=False)

    # salva parâmetros relevantes para uso posterior
    reg = _cell_registry(_self)
    reg.new(mol_obj, structure.lattice, repeats)

    # calcula as cargas, cria a supercélula e converte em molécula
    structure = BVAnalyzer().get_oxi_state_decorated_structure(structure)
    molecule = structure_to_molecule(structure, repeats)

    # cria um arquivo temporário e carrega no PyMOL
    tmp = tempfile.NamedTemporaryFile(suffix=".mol")
    molecule.to(tmp.name)
    _self.load(tmp.name, mol_obj, **kwargs)
    tmp.close()

    # ajusta os raios iônicos e cargas formais
    for i, specie in enumerate(molecule.species, 1):
        try:
            ionic_radius = max(float(specie.ionic_radius), 1.0)
            _self.alter(
                f"{mol_obj} and index {i}",
                f"vdw={ionic_radius}; formal_charge={int(specie.oxi_state)}",
            )
        except Exception:
            continue

    # ajusta as ligações
    _self.unbond(mol_obj, mol_obj)
    for (i, si), (j, sj) in combinations(enumerate(molecule.sites, 1), 2):
        if int(si.specie.oxi_state) * int(sj.specie.oxi_state) > 0:
            continue

        try:
            cutoff = si.specie.ionic_radius + sj.specie.ionic_radius
            if si.distance(sj) < cutoff:
                _self.bond(f"{mol_obj} and index {i}", f"{mol_obj} and index {j}")
        except:
            continue

    _self.rebuild(mol_obj)


@cmd.extend
def add_cell(
    mol_obj=None,
    cll_obj=None,
    color="black",
    linewidth=0.02,
    _self=cmd,
):
    """
    Adiciona: estrutura da célula.

    >>> PyMOL> add_cell [ mol_obj [, cll_obj [, color [, linewidth ]]]]
    """
    mol_obj = mol_obj or _self.get_object_list()[0]
    cll_obj = cll_obj or f"{mol_obj}_cell"

    # recupera os parâmetros da rede
    reg = _cell_registry(_self)
    cell_info = reg.get(mol_obj)
    lattice = cell_info.lattice
    repeats = cell_info.repeats

    fracs = np.indices((2, 2, 2)).reshape(3, -1).T
    corners = lattice.get_cartesian_coords(fracs)
    edges = [
        (i, j)
        for i, j in combinations(range(8), 2)
        if np.sum(np.abs(fracs[i] - fracs[j])) == 1
    ]

    idxs = np.indices(repeats).reshape(3, -1).T
    offsets = idxs.dot(lattice.matrix)

    rgb = _self.get_color_tuple(color)
    cell_cgo = []
    for off in offsets:
        for i, j in edges:
            p1 = corners[i] + off
            p2 = corners[j] + off
            cell_cgo.extend([cgo.CYLINDER, *p1, *p2, linewidth, *rgb, *rgb])

    _self.load_cgo(cell_cgo, cll_obj)
    _self.show("cgo", cll_obj)
