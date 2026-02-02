# save_as_tikz.py
# Export a PyMOL selection to TikZ using \molecule, \atom, \bond (2D only).
# Units: Angstroms (Å), with perspective scaling taken from the current PyMOL view.
#
# Output format:
#   \begin{tikzpicture}
#   \begin{molecule}
#     \coordinate (C1) at (...,...);
#     ...
#     \atom[r=...,<color>]{C1};
#     \bond[r=...,p=...,q=...]{C1}{H1};
#     ...
#   \end{molecule}
#   \end{tikzpicture}

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from pymol import cmd, stored

from braunmol.core import ELEMENTS

DEFAULT_COLOR = "gray"
DEPTH_TOL = 0.10  # camera-space z binning (Å)
CLAMP_FRAC = 0.95  # keep p+q <= frac * projected bond length
STATE = 1  # PyMOL state to read settings from (keep it simple)


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
        return cam[:2] * self.scale_at(cam)


# ───────────────────── settings collection ─────────────────────


@dataclass(frozen=True, slots=True)
class AtomSettings:
    sphere_scale: float
    stick_radius: float
    stick_h_scale: float


def _collect_atom_settings(selection: str) -> Dict[int, AtomSettings]:
    """
    Collect per-atom effective settings using iterate_state() and the setting-proxy `s`.
    Keys are PyMOL atom 'index' (same as iterate's `index` and model.atom[i].index).
    """
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
    """
    Try to read a per-bond setting (e.g. set_bond stick_radius, ...).
    If unsupported or unset, return None.

    Note: this is called per bond; fine for typical molecules.
    """
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


# ───────────────────── molecule primitives ─────────────────────


@dataclass(slots=True)
class Atom2D:
    idx: int  # PyMOL atom index
    sym: str
    xyz: np.ndarray  # model coords (Å)
    vdw: float  # VDW radius (Å)

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
        # effective sphere radius in projected Å
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
    stick_radius_override: Optional[float] = None  # Å

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
        """
        Stick radius in projected Å, including stick_h_scale for X-H bonds and
        perspective scaling (average of endpoints).
        """
        base = self.stick_radius_override
        if base is None:
            base = 0.5 * (self.a.stick_radius + self.b.stick_radius)

        # Apply stick_h_scale when hydrogen is involved (PyMOL behavior)
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


# ───────────────────── exporter ─────────────────────


class TikzMolecule2D:
    def __init__(self, selection: str = "all"):
        self.selection = selection
        self.atoms: List[Atom2D] = []
        self.bonds: List[Bond2D] = []

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

        # Stable naming: element-count in ascending PyMOL index order (C1, H1, H2, ...)
        counts: Dict[str, int] = {}
        for at in sorted(atoms, key=lambda x: x.idx):
            counts[at.sym] = counts.get(at.sym, 0) + 1
            at.name = f"{at.sym}{counts[at.sym]}"

        self.atoms = atoms

        idx_to_atom = {a.idx: a for a in self.atoms}
        self.bonds = []

        # Prefer model.bond so we don't depend on cmd.get_bonds index conventions.
        for b in getattr(model, "bond", []):
            i0, i1 = int(b.index[0]), int(b.index[1])  # indices into model.atom list
            a_idx = int(model.atom[i0].index)
            b_idx = int(model.atom[i1].index)

            a_atom = idx_to_atom.get(a_idx)
            b_atom = idx_to_atom.get(b_idx)
            if a_atom is None or b_atom is None:
                continue

            stick_override = _try_get_bond_setting(
                "stick_radius", self.selection, a_idx, b_idx
            )
            self.bonds.append(
                Bond2D(a=a_atom, b=b_atom, stick_radius_override=stick_override)
            )

    def export_molecule_env(self) -> str:
        view = View.from_pymol()
        for a in self.atoms:
            a.prepare(view)

        atoms_sorted = sorted(self.atoms, key=lambda a: float(a.cam[2]))

        # group by quantized camera z
        layers: List[List[Atom2D]] = []
        cur_key: Optional[int] = None
        for a in atoms_sorted:
            k = round(float(a.cam[2]) / DEPTH_TOL)
            if cur_key is None or k != cur_key:
                layers.append([a])
                cur_key = k
            else:
                layers[-1].append(a)

        lines: List[str] = []
        lines.append(f"\\begin{{molecule}}[angstrom={stored.tikz_angstrom_cm:.3f}cm]")

        for a in atoms_sorted:
            lines.append(a.coord_line())

        seen: set[Tuple[str, str]] = set()
        for layer in layers:
            for a in layer:
                lines.append(a.atom_line())

            for bond in self._bonds_touching(layer):
                key = tuple(sorted((bond.a.name, bond.b.name)))
                if key in seen:
                    continue
                lines.append(bond.tikz_line())
                seen.add(key)

        lines.append("\\end{molecule}")
        return "\n".join(lines)

    def _bonds_touching(self, layer: List[Atom2D]) -> List[Bond2D]:
        layer_ids = {a.idx for a in layer}
        return [b for b in self.bonds if b.a.idx in layer_ids or b.b.idx in layer_ids]

    def write(self, path: Path) -> None:
        body = self.export_molecule_env()
        out = "\\begin{tikzpicture}\n" + body + "\n\\end{tikzpicture}\n"
        path.write_text(out, encoding="utf-8")


# ───────────────────── PyMOL hook ─────────────────────


@cmd.extend
def save_tikz(selection: str = "all", file_name: str | None = None) -> None:
    """
    PyMOL:
      save_tikz [selection], [file_name]

    Examples:
      save_tikz all, methane
      save_tikz (organic), organic_mol.tex
    """
    if not file_name:
        objs = cmd.get_object_list(selection) or cmd.get_object_list() or ["molecule"]
        file_name = objs[0]

    path = Path(file_name)
    if path.suffix.lower() != ".tex":
        path = path.with_suffix(".tex")

    mol = TikzMolecule2D(selection)
    mol.load()
    mol.write(path)
