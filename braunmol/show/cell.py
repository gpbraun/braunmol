"""
cell.py

Gabriel Braun, 2026
"""

import tempfile
from dataclasses import dataclass
from itertools import combinations, product
from pathlib import Path

import numpy as np
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.graphs import StructureGraph
from pymatgen.analysis.local_env import CrystalNN, MinimumVIRENN
from pymatgen.core.structure import Molecule, Structure
from pymol import cgo, cmd

# =============================================================================
# Registry (lattice + repeats)
# =============================================================================


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
    reg = getattr(_self, "_cell_registry", None)
    if reg is None:
        reg = CellRegistry()
        setattr(_self, "_cell_registry", reg)
    return reg


# =============================================================================
# PyMOL: load_cell
# =============================================================================


@cmd.extend
def load_cell(
    file_path,
    mol_obj=None,
    reps="1-1-1",
    _self=cmd,
    **kwargs,
):
    """
    >>> PyMOL> load_cell file_path [, mol_obj [, reps ]]
    """
    mol_obj = mol_obj or Path(file_path).stem
    repeats = np.fromstring(reps, sep="-", count=3, dtype=int)

    structure = Structure.from_file(file_path, check_cif=False)

    reg = _cell_registry(_self)
    reg.new(mol_obj, structure.lattice, repeats)

    # oxi-states (ajuda CrystalNN em iônico/misto)
    try:
        structure = BVAnalyzer().get_oxi_state_decorated_structure(structure)
    except Exception:
        pass

    # supercélula periódica correta
    sc = structure.copy()
    sc.make_supercell(repeats)

    # estratégia NN
    nn = CrystalNN()
    try:
        nn.get_nn_info(sc, 0)
    except Exception:
        nn = MinimumVIRENN()

    # grafo de ligações (pymatgen pronto)
    sg = StructureGraph.with_local_env_strategy(sc, nn)
    n = len(sc)

    # ------------------------------------------------------------
    # PB images: gera tuplas (0/1) para coordenadas ~0
    # ------------------------------------------------------------
    def pb_images(frac_coords: np.ndarray):
        zeros = np.flatnonzero(np.isclose(frac_coords, 0.0))
        if len(zeros) == 0:
            yield (0, 0, 0)
            return
        for mask in product((0, 1), repeat=len(zeros)):
            img = np.zeros(3, dtype=int)
            img[zeros] = mask
            yield (int(img[0]), int(img[1]), int(img[2]))

    # ------------------------------------------------------------
    # Normaliza uma ligação periódica (i -> j com jimage possivelmente negativa)
    # para uma representação equivalente só com imagens não-negativas (0/1)
    #
    # shift = max(-jimage, 0)
    # i@(shift) -- j@(jimage + shift)
    # ------------------------------------------------------------
    def normalize_image(jimage):
        jx, jy, jz = (int(x) for x in jimage)
        sx = -jx if jx < 0 else 0
        sy = -jy if jy < 0 else 0
        sz = -jz if jz < 0 else 0
        return (sx, sy, sz), (jx + sx, jy + sy, jz + sz)

    # ------------------------------------------------------------
    # 1) arestas base (únicas / não-direcionadas) usando o grafo
    #    key = (i,img_i0,j,img_j0) já normalizado p/ imagens >=0
    # ------------------------------------------------------------
    base_edges = set()
    for i in range(n):
        for cs in sg.get_connected_sites(i):
            j = int(cs.index)
            img_i0, img_j0 = normalize_image(cs.jimage)
            a = (i, img_i0)
            b = (j, img_j0)
            if a <= b:
                base_edges.add((i, img_i0, j, img_j0))
            else:
                base_edges.add((j, img_j0, i, img_i0))

    # ------------------------------------------------------------
    # 2) imagens existentes por site:
    #    - sempre (0,0,0)
    #    - PB images
    #    - imagens mínimas exigidas pelas ligações base
    # ------------------------------------------------------------
    images = [set([(0, 0, 0)]) for _ in range(n)]

    for i, site in enumerate(sc.sites):
        images[i].update(pb_images(site.frac_coords))

    for i, img_i0, j, img_j0 in base_edges:
        images[i].add(img_i0)
        images[j].add(img_j0)

    # ------------------------------------------------------------
    # 3) Molecule + idx_map (site, image) -> índice PyMOL (1-based)
    # ------------------------------------------------------------
    lattice = sc.lattice
    species_list = []
    coords_list = []
    idx_map = {}

    for i, site in enumerate(sc.sites):
        specie = site.specie if hasattr(site, "specie") else site.species
        f0 = site.frac_coords
        for img in sorted(images[i]):
            coords = lattice.get_cartesian_coords(f0 + np.array(img, dtype=float))
            species_list.append(specie)
            coords_list.append(coords)
            idx_map[(i, img)] = len(species_list)

    molecule = Molecule(species_list, coords_list)

    # ------------------------------------------------------------
    # 4) carrega no PyMOL
    # ------------------------------------------------------------
    try:
        _self.delete(mol_obj)
    except Exception:
        pass

    tmp = tempfile.NamedTemporaryFile(suffix=".mol")
    molecule.to(tmp.name)
    _self.load(tmp.name, mol_obj, **kwargs)
    tmp.close()

    # vdw + formal_charge (quando disponíveis)
    for k, specie in enumerate(molecule.species, 1):
        try:
            ionic_radius = max(float(specie.ionic_radius), 1.0)
            _self.alter(
                f"{mol_obj} and index {k}",
                f"vdw={ionic_radius}; formal_charge={int(specie.oxi_state)}",
            )
        except Exception:
            continue

    # ------------------------------------------------------------
    # 5) bonds: replica as ligações base para todas as PB images
    # ------------------------------------------------------------
    _self.unbond(mol_obj, mol_obj)

    bond = _self.bond
    idx_get = idx_map.get
    seen = set()

    for i, img_i0, j, img_j0 in base_edges:
        # replica a ligação para todas as imagens do átomo i
        for img_i in images[i]:
            tx = img_i[0] - img_i0[0]
            ty = img_i[1] - img_i0[1]
            tz = img_i[2] - img_i0[2]
            img_j = (img_j0[0] + tx, img_j0[1] + ty, img_j0[2] + tz)

            if img_j not in images[j]:
                continue

            a = idx_get((i, img_i))
            b = idx_get((j, img_j))
            if a is None or b is None or a == b:
                continue

            if a > b:
                a, b = b, a
            if (a, b) in seen:
                continue
            seen.add((a, b))

            bond(f"{mol_obj} and index {a}", f"{mol_obj} and index {b}")

    _self.rebuild(mol_obj)


# =============================================================================
# PyMOL: add_cell
# =============================================================================


@cmd.extend
def add_cell(
    mol_obj=None,
    cll_obj=None,
    color="black",
    linewidth=0.02,
    _self=cmd,
):
    """
    >>> PyMOL> add_cell [ mol_obj [, cll_obj [, color [, linewidth ]]]]
    """
    mol_obj = mol_obj or _self.get_object_list()[0]
    cll_obj = cll_obj or f"{mol_obj}_cell"

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
