# save_as_tikz.py
# Export a PyMOL selection to TikZ using \molecule, \atom, \bond (2D only).
# Units: Angstroms (Å), with perspective scaling taken from the current PyMOL view.
#
# Optional: export the unit cell (from braunmol.show.cell._cell_registry) as simple paths:
#   \draw[cell] (x1,y1) -- (x2,y2);
#
# Output format:
#   \begin{tikzpicture}
#   \begin{molecule}[angstrom=...cm]
#     \coordinate (C1) at (...,...);
#     ...
#     \draw[cell] (...,...) -- (...,...);   % optional
#     \atom[r=...,<color>]{C1};
#     \bond[r=...,p=...,q=...]{C1}{H1};
#     ...
#   \end{molecule}
#   \end{tikzpicture}

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from pymol import cmd, stored

from braunmol.core import ELEMENTS
from braunmol.show.cell import _cell_registry

DEFAULT_COLOR = "gray"
DEPTH_TOL = 0.10  # camera-space z binning (Å)
CLAMP_FRAC = 0.95  # keep p+q <= frac * projected bond length
STATE = 1  # state to read settings from
CELL_STYLE = "cell"  # TikZ style name used in \draw[<style>] for the unit cell


# ───────────────────────── view ──────────────────────────


@dataclass(frozen=True, slots=True)
class View:
    R: np.ndarray
    o_model: np.ndarray
    o_cam: np.ndarray
    z_near: float
    perspective: bool

    @staticmethod
    def from_pymol() -> "View":
        v = np.asarray(cmd.get_view(), dtype=float)
        return View(
            R=v[:9].reshape(3, 3),
            o_cam=v[9:12],
            o_model=v[12:15],
            z_near=float(v[15]),
            perspective=bool(v[17] < 0),
        )

    def cam_xyz(self, xyz: np.ndarray) -> np.ndarray:
        return (xyz - self.o_model) @ self.R + self.o_cam

    def scale_at(self, cam: np.ndarray) -> float:
        if not self.perspective:
            return 1.0
        z = float(cam[2])
        if abs(z) < 1e-9:
            return 1.0
        return float(self.z_near / (-z))

    def project_xy(self, cam: np.ndarray) -> np.ndarray:
        s = self.scale_at(cam)
        return cam[:2] * s


# ───────────────────── settings collection ─────────────────────


@dataclass(frozen=True, slots=True)
class AtomSettings:
    sphere_scale: float
    stick_radius: float
    stick_h_scale: float


def _collect_atom_settings(selection: str) -> Dict[int, AtomSettings]:
    rows: List[Tuple[int, float, float, float]] = []
    cmd.iterate_state(
        STATE,
        selection,
        "rows.append((index, s.sphere_scale, s.stick_radius, s.stick_h_scale))",
        space={"rows": rows},
        quiet=1,
    )

    out: Dict[int, AtomSettings] = {}
    for idx, sphere_scale, stick_radius, stick_h_scale in rows:
        out[int(idx)] = AtomSettings(
            sphere_scale=float(sphere_scale),
            stick_radius=float(stick_radius),
            stick_h_scale=float(stick_h_scale),
        )
    return out


def _try_get_bond_setting(
    name: str, selection: str, a_idx: int, b_idx: int
) -> Optional[float]:
    get_bond = getattr(cmd, "get_bond", None)
    if get_bond is None:
        return None

    sel1 = f"({selection}) and index {a_idx}"
    sel2 = f"({selection}) and index {b_idx}"

    try:
        val = get_bond(name, sel1, sel2, quiet=1)
    except Exception:
        return None

    if isinstance(val, (list, tuple)):
        if not val:
            return None
        try:
            return float(val[0])
        except Exception:
            return None

    try:
        return float(val)
    except Exception:
        return None


def _tikz_color(sym: str) -> str:
    return str(ELEMENTS.get(sym, {}).get("color", DEFAULT_COLOR)).replace("_", "-")


def _layer_key(z: float) -> int:
    return int(round(float(z) / DEPTH_TOL))


def _begin_molecule_line() -> str:
    ang = getattr(stored, "tikz_angstrom_cm", None)
    try:
        ang_f = float(ang)
    except Exception:
        ang_f = float("nan")

    if np.isfinite(ang_f) and ang_f > 0:
        return f"\\begin{{molecule}}[angstrom={ang_f:.3f}cm]"
    return "\\begin{molecule}"


# ───────────────────── primitives ─────────────────────


@dataclass(slots=True)
class Atom2D:
    idx: int
    sym: str
    xyz: np.ndarray
    vdw: float
    sphere_scale: float
    stick_radius: float
    stick_h_scale: float

    name: str = ""
    color: str = ""
    cam: np.ndarray = field(default_factory=lambda: np.zeros(3))
    proj: np.ndarray = field(default_factory=lambda: np.zeros(2))
    cam_scale: float = 1.0

    def prepare(self, view: View) -> None:
        self.color = _tikz_color(self.sym)
        self.cam = view.cam_xyz(self.xyz)
        self.cam_scale = view.scale_at(self.cam)
        self.proj = self.cam[:2] * self.cam_scale

    @property
    def sphere_r(self) -> float:
        return float(self.vdw * self.sphere_scale * self.cam_scale)

    def coord_line(self) -> str:
        x, y = self.proj
        return f"  \\coordinate ({self.name}) at ({x:.3f},{y:.3f});"

    def atom_line(self) -> str:
        return f"  \\atom[r={self.sphere_r:.3f},{self.color}]{{{self.name}}};"


@dataclass(frozen=True, slots=True)
class Bond2D:
    a: Atom2D
    b: Atom2D
    stick_radius_override: Optional[float] = None

    def _sin_theta(self) -> float:
        d = self.b.cam - self.a.cam
        denom = float(np.linalg.norm(d)) or 1.0
        return float(np.hypot(d[0], d[1]) / denom)

    def pq_offsets(self) -> Tuple[float, float]:
        d2 = self.b.proj - self.a.proj
        L = float(np.linalg.norm(d2))
        if L <= 1e-12:
            return 0.0, 0.0

        sin_t = self._sin_theta()
        p = float(self.a.sphere_r * sin_t)
        q = float(self.b.sphere_r * sin_t)

        max_sum = CLAMP_FRAC * L
        s = p + q
        if s > max_sum and s > 1e-12:
            k = max_sum / s
            p *= k
            q *= k

        return p, q

    def thickness_r(self) -> float:
        base = self.stick_radius_override
        if base is None:
            base = 0.5 * (self.a.stick_radius + self.b.stick_radius)

        if self.a.sym == "H" and self.b.sym != "H":
            base *= self.a.stick_h_scale
        elif self.b.sym == "H" and self.a.sym != "H":
            base *= self.b.stick_h_scale
        elif self.a.sym == "H" and self.b.sym == "H":
            base *= 0.5 * (self.a.stick_h_scale + self.b.stick_h_scale)

        cam_scale = 0.5 * (self.a.cam_scale + self.b.cam_scale)
        return float(base * cam_scale)

    def tikz_line(self) -> str:
        p, q = self.pq_offsets()
        r = self.thickness_r()
        return f"  \\bond[r={r:.3f},p={p:.3f},q={q:.3f}]{{{self.a.name}}}{{{self.b.name}}};"


@dataclass(frozen=True, slots=True)
class CellEdge2D:
    p1: np.ndarray
    p2: np.ndarray
    zmid: float

    def tikz_line(self) -> str:
        x1, y1 = self.p1
        x2, y2 = self.p2
        return f"  \\draw[{CELL_STYLE}] ({x1:.3f},{y1:.3f}) -- ({x2:.3f},{y2:.3f});"


# ───────────────────── unit cell extraction ─────────────────────


def _cell_edges_for_object(obj_name: str, view: View) -> List[CellEdge2D]:
    """
    Builds the same edges as add_cell(), but emits projected 2D segments.

    Requires:
      reg = _cell_registry(cmd)
      cell_info = reg.get(obj_name)
      cell_info.lattice, cell_info.repeats
    """
    reg = _cell_registry(cmd)
    cell_info = reg.get(obj_name)
    if cell_info is None:
        return []

    lattice = getattr(cell_info, "lattice", None)
    repeats = getattr(cell_info, "repeats", None)
    if lattice is None or repeats is None:
        return []

    repeats = tuple(int(x) for x in repeats)

    fracs = np.indices((2, 2, 2)).reshape(3, -1).T
    corners = lattice.get_cartesian_coords(fracs)

    edges = [
        (i, j)
        for i, j in combinations(range(8), 2)
        if int(np.sum(np.abs(fracs[i] - fracs[j]))) == 1
    ]

    idxs = np.indices(repeats).reshape(3, -1).T
    offsets = idxs.dot(lattice.matrix)

    out: List[CellEdge2D] = []
    for off in offsets:
        for i, j in edges:
            P1 = corners[i] + off
            P2 = corners[j] + off

            cam1 = view.cam_xyz(P1)
            cam2 = view.cam_xyz(P2)

            p1 = cam1[:2] * view.scale_at(cam1)
            p2 = cam2[:2] * view.scale_at(cam2)

            zmid = 0.5 * (float(cam1[2]) + float(cam2[2]))
            out.append(CellEdge2D(p1=p1, p2=p2, zmid=zmid))

    return out


# ───────────────────── exporter ─────────────────────


class TikzMolecule2D:
    def __init__(self, selection: str = "all", *, cell: bool = False):
        self.selection = selection
        self.cell = bool(cell)

        self.atoms: List[Atom2D] = []
        self.bonds: List[Bond2D] = []
        self._adj: Dict[int, List[Bond2D]] = {}

        self._cell_obj: Optional[str] = None

    def load(self) -> None:
        model = cmd.get_model(self.selection)
        aset = _collect_atom_settings(self.selection)

        atoms: List[Atom2D] = []
        for a in model.atom:
            idx = int(a.index)
            s = aset.get(idx, AtomSettings(1.0, 0.25, 1.0))
            atoms.append(
                Atom2D(
                    idx=idx,
                    sym=str(a.symbol),
                    xyz=np.asarray(a.coord, dtype=float),
                    vdw=float(getattr(a, "vdw", 1.5)),
                    sphere_scale=s.sphere_scale,
                    stick_radius=s.stick_radius,
                    stick_h_scale=s.stick_h_scale,
                )
            )

        counts: Dict[str, int] = {}
        for at in sorted(atoms, key=lambda x: x.idx):
            counts[at.sym] = counts.get(at.sym, 0) + 1
            at.name = f"{at.sym}{counts[at.sym]}"

        self.atoms = atoms

        idx_to_atom = {a.idx: a for a in self.atoms}
        self.bonds = []
        self._adj = {a.idx: [] for a in self.atoms}

        for b in getattr(model, "bond", []):
            i0, i1 = int(b.index[0]), int(b.index[1])
            a_idx = int(model.atom[i0].index)
            b_idx = int(model.atom[i1].index)

            a_atom = idx_to_atom.get(a_idx)
            b_atom = idx_to_atom.get(b_idx)
            if a_atom is None or b_atom is None:
                continue

            stick_override = _try_get_bond_setting(
                "stick_radius", self.selection, a_idx, b_idx
            )
            bond = Bond2D(a=a_atom, b=b_atom, stick_radius_override=stick_override)
            self.bonds.append(bond)
            self._adj[a_idx].append(bond)
            self._adj[b_idx].append(bond)

        if self.cell:
            objs = cmd.get_object_list(self.selection) or cmd.get_object_list() or []
            self._cell_obj = objs[0] if objs else None

    def export_molecule_env(self) -> str:
        view = View.from_pymol()
        for a in self.atoms:
            a.prepare(view)

        atoms_sorted = sorted(self.atoms, key=lambda a: float(a.cam[2]))

        cell_edges: List[CellEdge2D] = []
        if self.cell and self._cell_obj:
            cell_edges = _cell_edges_for_object(self._cell_obj, view)

        atom_layers: Dict[int, List[Atom2D]] = {}
        for a in atoms_sorted:
            atom_layers.setdefault(_layer_key(a.cam[2]), []).append(a)

        cell_layers: Dict[int, List[CellEdge2D]] = {}
        for e in cell_edges:
            cell_layers.setdefault(_layer_key(e.zmid), []).append(e)

        keys = sorted(set(atom_layers.keys()) | set(cell_layers.keys()))

        lines: List[str] = []
        lines.append(_begin_molecule_line())

        for a in atoms_sorted:
            lines.append(a.coord_line())

        seen_bonds: set[Tuple[str, str]] = set()

        for k in keys:
            # cell behind atoms/bonds in same layer
            for e in cell_layers.get(k, []):
                lines.append(e.tikz_line())

            layer_atoms = atom_layers.get(k, [])
            for a in layer_atoms:
                lines.append(a.atom_line())

            for a in layer_atoms:
                for bond in self._adj.get(a.idx, []):
                    key = tuple(sorted((bond.a.name, bond.b.name)))
                    if key in seen_bonds:
                        continue
                    lines.append(bond.tikz_line())
                    seen_bonds.add(key)

        lines.append("\\end{molecule}")
        return "\n".join(lines)

    def write(self, path: Path) -> None:
        body = self.export_molecule_env()
        out = "\\begin{tikzpicture}\n" + body + "\n\\end{tikzpicture}\n"
        path.write_text(out, encoding="utf-8")


# ───────────────────── PyMOL hook ─────────────────────


@cmd.extend
def save_tikz(
    selection: str = "all",
    file_name: str | None = None,
    cell: int = 0,
) -> None:
    """
    PyMOL:
      save_tikz [selection], [file_name], [cell]

    cell:
      0 -> no unit cell
      1 -> include unit cell edges as \\draw[cell] ...
    """
    if not file_name:
        objs = cmd.get_object_list(selection) or cmd.get_object_list() or ["molecule"]
        file_name = objs[0]

    path = Path(file_name)
    if path.suffix.lower() != ".tex":
        path = path.with_suffix(".tex")

    mol = TikzMolecule2D(selection, cell=bool(cell))
    mol.load()
    mol.write(path)
