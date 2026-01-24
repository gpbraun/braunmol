# save_as_tikz.py
# Export the current PyMOL selection to TikZ (flat 2-D or true 3-D).

from __future__ import annotations

import importlib
import itertools
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from pymol import cmd

from braunmol.core import ELEMENTS

# ───────────────────── constants ──────────────────────
ANG_TO_CM = 0.7  # 1 Å → 0.7 cm on page
SPHERE_K = 0.2  # multiplies PyMOL VDW radius
DEPTH_TOL = 0.10  # z-tolerance for painte's layers
DEFAULT_COLOR = "gray"


# ─────────────── helper utilities ─────────────────────
def _pymol_view():
    """Return (rotation, origin_model, origin_camera, z_near, perspective?)."""
    view = np.array(cmd.get_view())
    rotation = view[:9].reshape(3, 3)
    origin_model = view[12:15]
    origin_camera = view[9:12]
    z_near = view[15]
    perspective = view[17] < 0
    return rotation, origin_model, origin_camera, z_near, perspective


def _view_angles(rotation):
    forward = -rotation[:, 2]
    az = np.degrees(np.arctan2(forward[1], forward[0]))
    el = np.degrees(np.arcsin(forward[2] / np.linalg.norm(forward)))
    return az, el


def _persp_vectors(rotation, z_near, length_cm=10.0):
    scale = z_near * length_cm / max(abs(z_near), 1e-6)
    return rotation[:, 0] * scale, rotation[:, 1] * scale, -rotation[:, 2] * scale


def _fmt_vec(vec):
    return f"({vec[0]:.3f},{vec[1]:.3f},{vec[2]:.3f})"


# ───────────────────── data objects ───────────────────
@dataclass(slots=True)
class TikzAtom:
    idx: int
    sym: str
    xyz: np.ndarray
    vdw: float

    cam: np.ndarray = field(init=False)
    proj: np.ndarray = field(init=False)
    px: float = field(init=False)
    col: str = field(init=False)

    def prepare(self, rotation, o_mod, o_cam, z_near, perspective):
        self.cam = (self.xyz - o_mod) @ rotation + o_cam
        x, y, z_val = self.cam
        scale = z_near / -z_val if perspective and abs(z_val) > 1e-6 else 1.0
        print(scale)

        self.proj = np.array([x, y]) * scale
        self.px = self.vdw * ANG_TO_CM * SPHERE_K * scale

        if self.sym == "C":
            self.col = "carbon"
        else:
            self.col = (
                ELEMENTS.get(self.sym, {}).get("color", DEFAULT_COLOR).replace("_", "-")
            )

    # TikZ helpers
    def coord_2d(self):
        return f"\\coordinate (A{self.idx}) at ({self.proj[0]:.3f},{self.proj[1]:.3f});"

    def coord_3d(self):
        x, y, z_val = self.xyz
        return f"\\coordinate (A{self.idx}) at ({x:.3f},{y:.3f},{z_val:.3f});"

    def coord_tpp(self):
        x, y, z_val = self.xyz
        return (
            f"\\coordinate (A{self.idx}) "
            f"at (tpp cs:x={x:.3f},y={y:.3f},z={z_val:.3f});"
        )

    def node(self):
        return (
            f"\\node (A{self.idx}) "
            f"[ball={{vdw={self.px:.3f},color={self.col}}}] "
            f"at (A{self.idx}) {{}};"
        )


@dataclass(slots=True)
class TikzBond:
    a: TikzAtom
    b: TikzAtom
    order: int

    def tikz(self):
        delta = self.b.cam - self.a.cam
        sin_t = np.hypot(delta[0], delta[1]) / (np.linalg.norm(delta) or 1.0)
        off_a = self.a.px * sin_t
        off_b = self.b.px * sin_t
        macro = "doublebond" if self.order == 2 else "bond"
        return (
            f"\\{macro}[{off_a:.3f}cm]{{A{self.a.idx}}}"
            f"[{off_b:.3f}cm]{{A{self.b.idx}}};"
        )


# ─────────────────────── molecule ─────────────────────
class TikzMolecule:
    def __init__(self, selection="all"):
        self.sel = selection
        self.atoms = []
        self.bonds = []
        self._adj = {}

    def load(self):
        model = cmd.get_model(self.sel)
        self.atoms = [
            TikzAtom(
                a.index,
                a.symbol,
                ANG_TO_CM * np.asarray(a.coord),
                getattr(a, "vdw", 1.5),
            )
            for a in model.atom
        ]
        idx_map = {a.idx: a for a in self.atoms}
        for i, j, order in cmd.get_bonds(self.sel):
            a = idx_map.get(i + 1)
            b = idx_map.get(j + 1)
            if a and b:
                bond = TikzBond(a, b, order)
                self.bonds.append(bond)
                self._adj.setdefault(a.idx, []).append(bond)
                self._adj.setdefault(b.idx, []).append(bond)

    # ------------- TikZ writer -------------
    def write(self, outfile, mode="2d"):
        if mode not in ("2d", "3d"):
            raise ValueError("mode must be '2d' or '3d'")

        r, o_mod, o_cam, z_near, perspective = _pymol_view()
        for a in self.atoms:
            a.prepare(r, o_mod, o_cam, z_near, perspective)

        atoms_sorted = sorted(self.atoms, key=lambda a: a.cam[2])
        layer_key = lambda a: round(a.cam[2] / DEPTH_TOL)
        layers = [list(g) for _, g in itertools.groupby(atoms_sorted, layer_key)]

        header, coord_fn = self._header_coord(mode, perspective, r, z_near)

        lines = []
        lines.extend(header)
        lines.append("%% coordinates")
        lines.extend(coord_fn(atom) for atom in atoms_sorted)
        lines.append("")

        seen = set()
        for layer in layers:
            lines.append("%% nodes")
            lines.extend(atom.node() for atom in layer)
            lines.append("%% bonds")
            for atom in layer:
                for bond in self._adj[atom.idx]:
                    key = tuple(sorted((bond.a.idx, bond.b.idx)))
                    if key not in seen:
                        lines.append(bond.tikz())
                        seen.add(key)
            # lines.append("")
        lines.append("\\end{tikzpicture}")

        Path(outfile).write_text("\n".join(lines), encoding="utf-8")

    # ----------- header & coord switcher -----------
    def _header_coord(self, mode, perspective, rotation, z_near):
        if mode == "2d":
            return [
                "\\begin{tikzpicture}[xscale=\\molscale,yscale=\\molscale]"
            ], TikzAtom.coord_2d

        az, el = _view_angles(rotation)
        header = [
            "\\begin{tikzpicture}[",
            f"  3d view={{{az:.1f}}}{{{el:.1f}}},",
        ]
        if perspective:
            p, q, rr = _persp_vectors(rotation, z_near)
            header.append(
                f"  perspective={{ p={{{_fmt_vec(p)}}}, "
                f"q={{{_fmt_vec(q)}}}, r={{{_fmt_vec(rr)}}} }}]"
            )
            coord = TikzAtom.coord_tpp
        else:
            header[-1] += "]"
            coord = TikzAtom.coord_3d
        return header, coord


# ─────────────────── PyMOL command hook ──────────────────
@cmd.extend
def save_as_tikz(file_name=None, mode="2d"):
    """
    save_as_tikz <selection> , <file.tex> , 0|1  (0 = 2-D, 1 = 3-D).
    """
    file_name = file_name or cmd.get_object_list()[0]
    file_path = Path(file_name).with_suffix(".tex")

    mol = TikzMolecule("all")
    mol.load()
    mol.write(file_path, mode)
