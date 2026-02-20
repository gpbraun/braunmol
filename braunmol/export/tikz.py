# braunmol/export/tikz.py
# Export a PyMOL selection to TikZ using \molecule, \atom, \bond, \cell (2D only).
# Coordinates are in projected Angstroms (Å). Perspective scaling is taken from cmd.get_view().

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
from pymol import cmd, stored

from braunmol.core import ELEMENTS
from braunmol.show.cell import _cell_registry

# =============================================================================
# Configuration knobs
# =============================================================================

STATE = 1
DEFAULT_COLOR = "gray"

# Ensure p+q <= CLAMP_FRAC * projected segment length (Å)
CLAMP_FRAC = 0.95

# Your TeX does:  r/.code = { \pgfmathsetlengthmacro\bond@r { 1.3*#1*\mol@unit} }
# We therefore pass r_param such that (1.3 * r_param) matches the projected stick radius.
LATEX_R_FACTOR = 1.3

# Visual calibration:
# - BOND_RADIUS_SCALE lets you tune thickness if PyMOL looks visually thicker/thinner vs TikZ.
# - CAP_CLEARANCE adds extra p/q to avoid rounded caps poking into spheres (set 0 for strict).
BOND_RADIUS_SCALE = 1.00
CAP_CLEARANCE = (
    0.00  # recommended starting point; increase slightly (0.2–0.6) if needed
)

# Cell tolerances (Å, camera space)
VERTEX_ATOM_TOL_A = 0.06
EDGE_ATOM_TOL_A = 0.05
CUT_MIN_SEP = 1e-3

# Mid-bond occlusion splitting (hide bond parts behind nearer spheres)
USE_BOND_OCCLUSION = True
OCCLUDER_PAD = 2e-3  # widen cut interval slightly to avoid hairline artifacts
OCCLUDER_INFLATE = 0.00  # DO NOT inflate by bond thickness for ball-and-stick look


# =============================================================================
# View / projection
# =============================================================================


@dataclass(frozen=True, slots=True)
class View:
    r"""
    Minimal view model derived from cmd.get_view().

    Notes:
      - get_view rotation is column-major; reshape with order="F".
      - In perspective: projected_xy = cam_xy * (z_near / -z_cam)
    """

    R: np.ndarray  # (3,3) model->cam rotation
    o_cam: np.ndarray  # (3,) camera translation vector
    o_model: np.ndarray  # (3,) model origin for rotation
    z_near: float
    perspective: bool

    @staticmethod
    def from_pymol() -> "View":
        v = np.asarray(cmd.get_view(), dtype=np.float64)
        R = v[:9].reshape((3, 3), order="F")
        o_cam = v[9:12].copy()
        o_model = v[12:15].copy()
        z_near = float(v[15])
        perspective = bool(v[17] < 0.0)
        return View(
            R=R, o_cam=o_cam, o_model=o_model, z_near=z_near, perspective=perspective
        )

    def cam_xyz(self, xyz: np.ndarray) -> np.ndarray:
        return self.R @ (xyz - self.o_model) + self.o_cam

    def scale_at_cam(self, cam: np.ndarray) -> float:
        if not self.perspective:
            return 1.0
        z = float(cam[2])
        if abs(z) < 1e-12:
            return 1.0
        return float(self.z_near / (-z))

    def project(self, xyz: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
        cam = self.cam_xyz(xyz)
        s = self.scale_at_cam(cam)
        return cam, cam[:2] * s, s

    def project_xy_from_cam(self, cam: np.ndarray) -> np.ndarray:
        return cam[:2] * self.scale_at_cam(cam)

    def cam_from_proj_param_on_segment(
        self,
        A_cam: np.ndarray,
        B_cam: np.ndarray,
        p0: np.ndarray,
        p1: np.ndarray,
        t2d: float,
    ) -> np.ndarray:
        """
        Given a param t2d along the *projected* segment p(t)=p0+t*(p1-p0),
        return a consistent point on the camera-space 3D segment A_cam+u*(B_cam-A_cam).

        In perspective, t2d != u; we solve u from x and y equations and average.
        """
        t2d = float(t2d)
        if not self.perspective:
            return A_cam + t2d * (B_cam - A_cam)

        p = p0 + t2d * (p1 - p0)
        near = float(self.z_near)

        x0, y0, z0 = map(float, A_cam)
        dx, dy, dz = map(float, (B_cam - A_cam))

        def solve_u(coord0: float, dcoord: float, proj_coord: float) -> Optional[float]:
            # proj = (coord0 + u*dcoord) * (near / -(z0 + u*dz))
            # => coord0 + u*dcoord = k*(z0 + u*dz), with k = -proj/near
            k = -float(proj_coord) / near
            denom = dcoord - k * dz
            if abs(denom) < 1e-14:
                return None
            return (k * z0 - coord0) / denom

        ux = solve_u(x0, dx, float(p[0]))
        uy = solve_u(y0, dy, float(p[1]))

        vals = [u for u in (ux, uy) if u is not None and np.isfinite(u)]
        if not vals:
            u = t2d
        elif len(vals) == 1:
            u = vals[0]
        else:
            u = 0.5 * (vals[0] + vals[1])

        u = 0.0 if u < 0.0 else (1.0 if u > 1.0 else float(u))
        return A_cam + u * (B_cam - A_cam)


def _begin_molecule_line() -> str:
    return r"\begin{molecule}"


# =============================================================================
# PyMOL settings + reps
# =============================================================================


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


def _collect_rep_indices(selection: str, rep: str) -> Set[int]:
    idxs: Set[int] = set()
    sel = f"({selection}) and rep {rep}"
    cmd.iterate_state(
        STATE,
        sel,
        "idxs.add(index)",
        space={"idxs": idxs},
        quiet=1,
    )
    return idxs


def _collect_bond_setting_map(
    name: str, selection: str
) -> Dict[Tuple[int, int], float]:
    """
    Build a (lo_idx, hi_idx) -> value map for a per-bond setting.

    cmd.get_bond commonly returns:
      [(model_name, [(idx1, idx2, value), ...]), ...]
    """
    get_bond = getattr(cmd, "get_bond", None)
    if get_bond is None:
        return {}

    out: Dict[Tuple[int, int], float] = {}
    try:
        res = get_bond(name, selection, selection, state=STATE, quiet=1)
    except Exception:
        return out

    for _model, vlist in res or []:
        for idx1, idx2, value in vlist or []:
            try:
                v = float(value)
            except Exception:
                continue
            a, b = int(idx1), int(idx2)
            lo, hi = (a, b) if a < b else (b, a)
            out[(lo, hi)] = v

    return out


def _tikz_color(sym: str) -> str:
    return str(ELEMENTS.get(sym, {}).get("color", DEFAULT_COLOR)).replace("_", "-")


# =============================================================================
# Shared geometry utilities
# =============================================================================


def _coord_line(name: str, xy: np.ndarray) -> str:
    return rf"  \coordinate ({name}) at ({float(xy[0]):.3f},{float(xy[1]):.3f});"


def _coord_alias(name: str, target: str) -> str:
    return rf"  \coordinate ({name}) at ({target});"


def _closest_param(p: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    ab = b - a
    denom = float(ab @ ab)
    if denom <= 1e-18:
        return 0.0
    t = float(((p - a) @ ab) / denom)
    return 0.0 if t < 0.0 else (1.0 if t > 1.0 else t)


def _dist_point_segment_3d(p: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    t = _closest_param(p, a, b)
    q = a + t * (b - a)
    return float(np.linalg.norm(p - q))


def _segment_circle_interval(
    p0: np.ndarray, p1: np.ndarray, c: np.ndarray, r: float
) -> Optional[Tuple[float, float]]:
    """
    Intersection interval [t1,t2] for segment p(t)=p0+t*(p1-p0) with disk |p-c|<=r.
    Returns None if no intersection.
    """
    d = p1 - p0
    f = p0 - c

    a = float(d @ d)
    if a <= 1e-18:
        return None

    b = 2.0 * float(f @ d)
    cc = float(f @ f) - float(r * r)

    disc = b * b - 4.0 * a * cc
    if disc <= 1e-14:
        return None

    sdisc = float(np.sqrt(disc))
    t1 = (-b - sdisc) / (2.0 * a)
    t2 = (-b + sdisc) / (2.0 * a)
    if t1 > t2:
        t1, t2 = t2, t1

    t1 = 0.0 if t1 < 0.0 else (1.0 if t1 > 1.0 else float(t1))
    t2 = 0.0 if t2 < 0.0 else (1.0 if t2 > 1.0 else float(t2))

    if (t2 - t1) <= 1e-10:
        return None
    return t1, t2


def _merge_intervals(iv: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    if not iv:
        return []
    iv.sort()
    out = [iv[0]]
    for a, b in iv[1:]:
        pa, pb = out[-1]
        if a <= pb:
            out[-1] = (pa, max(pb, b))
        else:
            out.append((a, b))
    return out


def _subtract_intervals(
    base: Tuple[float, float], cut: List[Tuple[float, float]]
) -> List[Tuple[float, float]]:
    a0, b0 = base
    if not cut:
        return [(a0, b0)]
    cut = _merge_intervals(cut)

    out: List[Tuple[float, float]] = []
    cur = a0
    for a, b in cut:
        if a > cur:
            out.append((cur, min(a, b0)))
        cur = max(cur, b)
        if cur >= b0:
            break
    if cur < b0:
        out.append((cur, b0))

    return [(a, b) for a, b in out if (b - a) > 1e-6]


# =============================================================================
# Nodes / atoms (2D projected)
# =============================================================================


@dataclass(frozen=True, slots=True)
class Node2D:
    r"""
    Endpoint for \bond or \cell.

    radius_proj: projected sphere radius in Å (0 for cuts/vertices or hidden spheres)
    radius_world: world sphere radius in Å (0 for cuts/vertices or hidden spheres)
    """

    name: str
    cam: np.ndarray
    proj: np.ndarray
    cam_scale: float
    radius_proj: float
    radius_world: float
    sym: str = ""


@dataclass(slots=True)
class Atom2D:
    idx: int
    sym: str
    xyz: np.ndarray
    vdw: float
    sphere_scale: float
    stick_radius: float
    stick_h_scale: float
    show_sphere: bool
    show_stick: bool

    name: str = ""
    color: str = ""

    cam: np.ndarray = field(default_factory=lambda: np.zeros(3))
    proj: np.ndarray = field(default_factory=lambda: np.zeros(2))
    cam_scale: float = 1.0

    _node_trim: Optional[Node2D] = field(default=None, init=False, repr=False)
    _node_notrim: Optional[Node2D] = field(default=None, init=False, repr=False)

    def prepare(self, view: View) -> None:
        self.color = _tikz_color(self.sym)
        self.cam, self.proj, self.cam_scale = view.project(self.xyz)
        self._node_trim = None
        self._node_notrim = None

    @property
    def radius_world(self) -> float:
        return float(self.vdw * self.sphere_scale)

    @property
    def radius_proj(self) -> float:
        return float(self.radius_world * self.cam_scale)

    @property
    def z_front(self) -> float:
        # Front-most surface point in camera space (+z toward viewer)
        return float(self.cam[2] + self.radius_world)

    def node(self, trim: bool) -> Node2D:
        if trim and self.show_sphere:
            if self._node_trim is None:
                self._node_trim = Node2D(
                    name=self.name,
                    cam=self.cam,
                    proj=self.proj,
                    cam_scale=self.cam_scale,
                    radius_proj=self.radius_proj,
                    radius_world=self.radius_world,
                    sym=self.sym,
                )
            return self._node_trim

        if self._node_notrim is None:
            self._node_notrim = Node2D(
                name=self.name,
                cam=self.cam,
                proj=self.proj,
                cam_scale=self.cam_scale,
                radius_proj=0.0,
                radius_world=0.0,
                sym=self.sym,
            )
        return self._node_notrim

    def coord_line(self) -> str:
        return _coord_line(self.name, self.proj)

    def atom_line(self) -> str:
        # Only called if show_sphere
        return rf"  \atom[r={self.radius_proj:.3f},{self.color}]{{{self.name}}};"


# =============================================================================
# 3D-correct p/q trimming (the key fix)
# =============================================================================


def _trim_len_from_sphere_3d(view: View, A: Node2D, B: Node2D) -> float:
    """
    Trimming length (projected Å) from endpoint A toward B based on the 3D
    exit point of the bond centerline from A's sphere.

    This fixes the over-trimming of tilted bonds that happens when using
    radius_proj directly.
    """
    ra = float(A.radius_world)
    if ra <= 0.0:
        return 0.0

    d3 = B.cam - A.cam
    L3 = float(np.linalg.norm(d3))
    if L3 <= 1e-12:
        return 0.0

    u = ra / L3
    if u >= 1.0:
        u = 1.0
    exit_cam = A.cam + u * d3
    exit_proj = view.project_xy_from_cam(exit_cam)

    d2 = B.proj - A.proj
    L2 = float(np.linalg.norm(d2))
    if L2 <= 1e-12:
        return 0.0

    u2 = d2 / L2
    p = float(np.dot(exit_proj - A.proj, u2))
    return 0.0 if p < 0.0 else p


def _pq_offsets_3d(
    view: View, A: Node2D, B: Node2D, r_draw: float
) -> Tuple[float, float]:
    """
    3D-correct p/q trimming in projected Å.

    r_draw is the actual drawn bond radius in projected units (after TeX factor),
    used only for optional cap clearance.
    """
    d2 = B.proj - A.proj
    L2 = float(np.linalg.norm(d2))
    if L2 <= 1e-12:
        return 0.0, 0.0

    p = _trim_len_from_sphere_3d(view, A, B)
    q = _trim_len_from_sphere_3d(view, B, A)

    cap = float(CAP_CLEARANCE) * float(r_draw)
    if A.radius_world > 0.0:
        p += cap
    if B.radius_world > 0.0:
        q += cap

    max_sum = CLAMP_FRAC * L2
    s = p + q
    if s > max_sum and s > 1e-12:
        k = max_sum / s
        p *= k
        q *= k

    return float(p), float(q)


# =============================================================================
# Cuts
# =============================================================================


class CutPool:
    def __init__(self, view: View):
        self.view = view
        self._n = 0
        self.nodes: Dict[str, Node2D] = {}

    def new(self, cam: np.ndarray) -> Node2D:
        self._n += 1
        name = f"cut{self._n}"
        proj = self.view.project_xy_from_cam(cam)
        s = self.view.scale_at_cam(cam)
        n = Node2D(
            name=name,
            cam=cam,
            proj=proj,
            cam_scale=s,
            radius_proj=0.0,
            radius_world=0.0,
            sym="",
        )
        self.nodes[name] = n
        return n

    def coord_lines(self, used: Set[str]) -> List[str]:
        return [
            _coord_line(name, self.nodes[name].proj)
            for name in sorted(used)
            if name in self.nodes
        ]


# =============================================================================
# Bonds (split by occlusion; p/q computed with 3D correction)
# =============================================================================


@dataclass(frozen=True, slots=True)
class BondSegment2D:
    a: Node2D
    b: Node2D
    r_param: float  # passed to \bond[r=...]
    r_draw: float  # drawn radius in projected units (LATEX_R_FACTOR * r_param)
    p: float
    q: float

    @property
    def z_key(self) -> float:
        # Use the front-most trimmed endpoint in camera space (better than midpoint)
        A = self.a.cam
        B = self.b.cam
        d = B - A
        n = float(np.linalg.norm(d))
        if n <= 1e-12:
            return float(A[2])

        u = d / n
        za = float(A[2] + float(self.a.radius_world) * float(u[2]))
        zb = float(B[2] - float(self.b.radius_world) * float(u[2]))
        return max(za, zb)

    def tikz_line(self) -> str:
        return rf"  \bond[r={self.r_param:.3f},p={self.p:.3f},q={self.q:.3f}]{{{self.a.name}}}{{{self.b.name}}};"


class BondBuilder:
    def __init__(
        self,
        view: View,
        atoms: Sequence[Atom2D],
        cut_pool: CutPool,
        stick_radius_map: Dict[Tuple[int, int], float],
    ):
        self.view = view
        self.atoms = list(atoms)
        self.cut_pool = cut_pool
        self.stick_radius_map = stick_radius_map

        self.idx_to_atom: Dict[int, Atom2D] = {a.idx: a for a in self.atoms}
        self.occluders: List[Atom2D] = [a for a in self.atoms if a.show_sphere]

    def _stick_radius_world(self, a: Atom2D, b: Atom2D) -> float:
        lo, hi = (a.idx, b.idx) if a.idx < b.idx else (b.idx, a.idx)
        base = self.stick_radius_map.get((lo, hi), None)
        if base is None:
            base = 0.5 * (a.stick_radius + b.stick_radius)
        base = float(abs(base))

        # Hydrogen scaling rule
        if a.sym == "H" and b.sym != "H":
            base *= float(a.stick_h_scale)
        elif b.sym == "H" and a.sym != "H":
            base *= float(b.stick_h_scale)
        elif a.sym == "H" and b.sym == "H":
            base *= 0.5 * (float(a.stick_h_scale) + float(b.stick_h_scale))

        return base

    def _bond_radii_projected(self, a: Atom2D, b: Atom2D) -> Tuple[float, float]:
        """
        Compute projected bond radius to match PyMOL cylinder radius under perspective.

        We use the average cam_scale. You can switch to max() if you prefer thicker
        near-camera appearance.
        """
        r_world = self._stick_radius_world(a, b)
        cam_scale = 0.5 * (float(a.cam_scale) + float(b.cam_scale))

        r_proj = float(r_world * cam_scale * float(BOND_RADIUS_SCALE))

        # TeX draws radius as (LATEX_R_FACTOR * r_param)
        r_param = float(r_proj / float(LATEX_R_FACTOR))
        r_draw = float(LATEX_R_FACTOR * r_param)
        return r_param, r_draw

    def _visible_intervals(
        self, a: Atom2D, b: Atom2D, r_draw: float
    ) -> List[Tuple[float, float]]:
        """
        Determine visible intervals along projected segment a->b by removing parts
        occluded by nearer spheres.

        IMPORTANT: For ball-and-stick, do not inflate the circle by bond thickness,
        or occlusion becomes too aggressive.
        """
        p0, p1 = a.proj, b.proj
        A_cam, B_cam = a.cam, b.cam

        cut: List[Tuple[float, float]] = []
        inflate = float(OCCLUDER_INFLATE) * float(r_draw)

        for sph in self.occluders:
            if sph.idx == a.idx or sph.idx == b.idx:
                continue

            R = float(sph.radius_proj + inflate)
            if R <= 1e-9:
                continue

            iv = _segment_circle_interval(p0, p1, sph.proj, R)
            if iv is None:
                continue

            t1, t2 = iv
            t1 = max(0.0, t1 - float(OCCLUDER_PAD))
            t2 = min(1.0, t2 + float(OCCLUDER_PAD))
            if (t2 - t1) <= 1e-6:
                continue

            tm = 0.5 * (t1 + t2)
            cam_m = self.view.cam_from_proj_param_on_segment(A_cam, B_cam, p0, p1, tm)
            z_line = float(cam_m[2])

            # occluder is in front if its front surface is closer (higher z)
            if float(sph.z_front) <= z_line:
                continue

            cut.append((t1, t2))

        return _subtract_intervals((0.0, 1.0), cut)

    def build(self, bond_pairs: Iterable[Tuple[int, int]]) -> List[BondSegment2D]:
        segs: List[BondSegment2D] = []

        for ia, ib in bond_pairs:
            a_atom = self.idx_to_atom.get(ia)
            b_atom = self.idx_to_atom.get(ib)
            if a_atom is None or b_atom is None:
                continue

            # Draw bond only if sticks shown for at least one endpoint
            if not (a_atom.show_stick or b_atom.show_stick):
                continue

            r_param, r_draw = self._bond_radii_projected(a_atom, b_atom)

            if USE_BOND_OCCLUSION:
                vis = self._visible_intervals(a_atom, b_atom, r_draw)
            else:
                vis = [(0.0, 1.0)]

            p0, p1 = a_atom.proj, b_atom.proj
            A_cam, B_cam = a_atom.cam, b_atom.cam

            def node_at(t: float) -> Node2D:
                if t <= 1e-12:
                    return a_atom.node(trim=True)
                if t >= 1.0 - 1e-12:
                    return b_atom.node(trim=True)
                cam_t = self.view.cam_from_proj_param_on_segment(
                    A_cam, B_cam, p0, p1, t
                )
                return self.cut_pool.new(cam_t)

            for t0, t1 in vis:
                n0 = node_at(t0)
                n1 = node_at(t1)
                if float(np.linalg.norm(n1.proj - n0.proj)) <= 1e-9:
                    continue

                # 3D-correct p/q on this *visible* segment
                p, q = _pq_offsets_3d(self.view, n0, n1, r_draw=r_draw)

                segs.append(
                    BondSegment2D(a=n0, b=n1, r_param=r_param, r_draw=r_draw, p=p, q=q)
                )

        return segs


# =============================================================================
# Cell exporting (optional)
# =============================================================================


@dataclass(frozen=True, slots=True)
class CellSegment:
    a: Node2D
    b: Node2D
    p: float
    q: float

    @property
    def z_key(self) -> float:
        return float(0.5 * (self.a.cam[2] + self.b.cam[2]))

    def tikz_line(self) -> str:
        return rf"  \cell[p={self.p:.3f},q={self.q:.3f}]{{{self.a.name}}}{{{self.b.name}}};"


@dataclass(frozen=True, slots=True)
class CellVertex:
    coord_name: str
    endpoint: Node2D  # atom node if occupied else the vertex node itself


class CellExporter:
    r"""
    Exports a registered unit cell as:
      - vertex coordinates: cell1, cell2, ...
      - \cell segments between endpoints (atoms, vertices, cuts)
    """

    def __init__(self, view: View, atoms: Sequence[Atom2D], cut_pool: CutPool):
        self.view = view
        self.cut_pool = cut_pool

        self.atoms_sphere = [a for a in atoms if a.show_sphere]
        self.atom_nodes = [a.node(trim=True) for a in self.atoms_sphere]

        self.vertices: List[CellVertex] = []
        self.segments: List[CellSegment] = []

    def build(self, obj_name: str) -> None:
        cell_info = _cell_registry(cmd).get(obj_name)
        if cell_info is None:
            return

        lattice = getattr(cell_info, "lattice", None)
        repeats = getattr(cell_info, "repeats", None)
        if lattice is None or repeats is None:
            return

        repeats = tuple(int(x) for x in repeats)

        fracs = np.indices((2, 2, 2)).reshape(3, -1).T
        corners0 = lattice.get_cartesian_coords(fracs)

        cube_edges = [
            (i, j)
            for i, j in combinations(range(8), 2)
            if int(np.sum(np.abs(fracs[i] - fracs[j]))) == 1
        ]

        idxs = np.indices(repeats).reshape(3, -1).T
        offsets = idxs.dot(lattice.matrix)

        v_key_to_id: Dict[Tuple[int, int, int], int] = {}
        v_xyz: List[np.ndarray] = []

        def key(xyz: np.ndarray) -> Tuple[int, int, int]:
            return (
                int(round(float(xyz[0]) * 1e6)),
                int(round(float(xyz[1]) * 1e6)),
                int(round(float(xyz[2]) * 1e6)),
            )

        def get_vid(xyz: np.ndarray) -> int:
            k = key(xyz)
            vid = v_key_to_id.get(k)
            if vid is None:
                vid = len(v_xyz)
                v_key_to_id[k] = vid
                v_xyz.append(np.asarray(xyz, dtype=float))
            return vid

        edge_pairs: Set[Tuple[int, int]] = set()
        for off in offsets:
            corners = corners0 + off
            vids = [get_vid(corners[i]) for i in range(8)]
            for i, j in cube_edges:
                a, b = vids[i], vids[j]
                lo, hi = (a, b) if a < b else (b, a)
                edge_pairs.add((lo, hi))

        endpoints: List[Node2D] = []
        self.vertices = []
        for i, xyz in enumerate(v_xyz, start=1):
            cam, proj, s = self.view.project(xyz)
            v_node = Node2D(
                name=f"cell{i}",
                cam=cam,
                proj=proj,
                cam_scale=s,
                radius_proj=0.0,
                radius_world=0.0,
                sym="",
            )

            occ = self._match_atom(v_node)
            endpoint = occ if occ is not None else v_node
            self.vertices.append(CellVertex(coord_name=v_node.name, endpoint=endpoint))
            endpoints.append(endpoint)

        self.segments = []
        for lo, hi in sorted(edge_pairs):
            self._split_edge(endpoints[lo], endpoints[hi])

    def coord_lines(self) -> List[str]:
        out: List[str] = []
        for v in self.vertices:
            if v.endpoint.name != v.coord_name:
                out.append(_coord_alias(v.coord_name, v.endpoint.name))
            else:
                out.append(_coord_line(v.coord_name, v.endpoint.proj))
        return out

    def _match_atom(self, v: Node2D) -> Optional[Node2D]:
        best: Optional[Node2D] = None
        best_d = float("inf")
        for a in self.atom_nodes:
            d = float(np.linalg.norm(a.cam - v.cam))
            if d < best_d:
                best_d = d
                best = a
        if best is not None and best_d <= VERTEX_ATOM_TOL_A:
            return best
        return None

    def _split_edge(self, A: Node2D, B: Node2D) -> None:
        A_cam, B_cam = A.cam, B.cam
        splits: List[Tuple[float, Node2D]] = [(0.0, A), (1.0, B)]

        # atoms lying on edge (3D) — only visible spheres
        for at in self.atoms_sphere:
            if _dist_point_segment_3d(at.cam, A_cam, B_cam) <= EDGE_ATOM_TOL_A:
                t = _closest_param(at.cam, A_cam, B_cam)
                if CUT_MIN_SEP < t < 1.0 - CUT_MIN_SEP:
                    splits.append((t, at.node(trim=True)))

        splits.sort(key=lambda x: x[0])

        compact: List[Tuple[float, Node2D]] = []
        for t, n in splits:
            if not compact or abs(t - compact[-1][0]) > CUT_MIN_SEP:
                compact.append((t, n))
            else:
                if n.radius_proj > compact[-1][1].radius_proj:
                    compact[-1] = (t, n)

        for (_, n1), (_, n2) in zip(compact, compact[1:]):
            if float(np.linalg.norm(n2.proj - n1.proj)) <= 1e-9:
                continue

            # Cell uses simple trimming (no cap clearance), but keep 3D-exit correction for consistency
            p = _trim_len_from_sphere_3d(self.view, n1, n2)
            q = _trim_len_from_sphere_3d(self.view, n2, n1)

            d2 = n2.proj - n1.proj
            L2 = float(np.linalg.norm(d2))
            if L2 > 1e-12:
                max_sum = CLAMP_FRAC * L2
                s = p + q
                if s > max_sum and s > 1e-12:
                    k = max_sum / s
                    p *= k
                    q *= k

            self.segments.append(CellSegment(a=n1, b=n2, p=float(p), q=float(q)))


# =============================================================================
# Main exporter
# =============================================================================


class TikzMolecule2D:
    def __init__(self, selection: str = "all"):
        self.selection = selection
        self.atoms: List[Atom2D] = []
        self.bond_pairs: List[Tuple[int, int]] = []
        self.obj_name: Optional[str] = None

    def load(self) -> None:
        model = cmd.get_model(self.selection, state=STATE)

        aset = _collect_atom_settings(self.selection)
        sphere_idxs = _collect_rep_indices(self.selection, "spheres")
        stick_idxs = _collect_rep_indices(self.selection, "sticks")

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
                    show_sphere=(idx in sphere_idxs),
                    show_stick=(idx in stick_idxs),
                )
            )

        # Stable naming: element-count in ascending atom index
        counts: Dict[str, int] = {}
        for at in sorted(atoms, key=lambda x: x.idx):
            counts[at.sym] = counts.get(at.sym, 0) + 1
            at.name = f"{at.sym}{counts[at.sym]}"

        self.atoms = atoms

        # Bond pairs from model.bond, de-duped
        pairs: Set[Tuple[int, int]] = set()
        for b in getattr(model, "bond", []) or []:
            i0, i1 = int(b.index[0]), int(b.index[1])  # indices into model.atom list
            a_idx = int(model.atom[i0].index)
            b_idx = int(model.atom[i1].index)
            if a_idx == b_idx:
                continue
            lo, hi = (a_idx, b_idx) if a_idx < b_idx else (b_idx, a_idx)
            pairs.add((lo, hi))
        self.bond_pairs = sorted(pairs)

        objs = cmd.get_object_list(self.selection) or cmd.get_object_list() or []
        self.obj_name = objs[0] if objs else None

    def export_molecule_env(self) -> str:
        view = View.from_pymol()
        for a in self.atoms:
            a.prepare(view)

        # sort far->near by camera z (painter order)
        atoms_sorted = sorted(self.atoms, key=lambda a: float(a.cam[2]))

        # shared cut pool
        cut_pool = CutPool(view)

        # bonds (optional)
        any_sticks = any(a.show_stick for a in atoms_sorted)
        stick_radius_map = _collect_bond_setting_map("stick_radius", self.selection)

        bond_segments: List[BondSegment2D] = []
        if any_sticks:
            bond_builder = BondBuilder(view, atoms_sorted, cut_pool, stick_radius_map)
            bond_segments = bond_builder.build(self.bond_pairs)

        # cell (optional)
        cell_coords: List[str] = []
        cell_segments: List[CellSegment] = []
        try:
            if self.obj_name and _cell_registry(cmd).get(self.obj_name) is not None:
                cell = CellExporter(view, atoms_sorted, cut_pool)
                cell.build(self.obj_name)
                cell_coords = cell.coord_lines()
                cell_segments = cell.segments
        except Exception:
            pass

        # Which cuts are referenced?
        used_cuts: Set[str] = set()
        for s in bond_segments:
            if s.a.name.startswith("cut"):
                used_cuts.add(s.a.name)
            if s.b.name.startswith("cut"):
                used_cuts.add(s.b.name)
        for s in cell_segments:
            if s.a.name.startswith("cut"):
                used_cuts.add(s.a.name)
            if s.b.name.startswith("cut"):
                used_cuts.add(s.b.name)

        # primitives for painter sort: (z_key, kind, line)
        # kind ensures consistent ordering for equal depth: cell(0) -> bonds(1) -> atoms(2)
        prims: List[Tuple[float, int, str]] = []
        prims.extend((s.z_key, 0, s.tikz_line()) for s in cell_segments)
        prims.extend((s.z_key, 1, s.tikz_line()) for s in bond_segments)
        prims.extend(
            (a.z_front, 2, a.atom_line()) for a in atoms_sorted if a.show_sphere
        )

        prims.sort(key=lambda x: (x[0], x[1]))

        # emit
        lines: List[str] = [_begin_molecule_line()]

        # atom coordinates always (bonds/cell aliasing depend on them)
        for a in atoms_sorted:
            lines.append(a.coord_line())

        # cell vertex coordinates (may alias atoms)
        lines.extend(cell_coords)

        # cut coordinates used by bonds/cell
        lines.extend(cut_pool.coord_lines(used=used_cuts))

        # primitives in depth order
        for _, __, line in prims:
            lines.append(line)

        lines.append(r"\end{molecule}")
        return "\n".join(lines)

    def write(self, path: Path) -> None:
        body = self.export_molecule_env()
        path.write_text(
            "\\begin{tikzpicture}\n" + body + "\n\\end{tikzpicture}\n",
            encoding="utf-8",
        )


# =============================================================================
# PyMOL command hook (wrap exceptions to avoid PyMOL parser crashes on Python 3.12)
# =============================================================================


@cmd.extend
def save_tikz(selection: str = "all", file_name: str | None = None) -> None:
    r"""
    PyMOL:
      save_tikz [selection], [file_name]

    Rules:
      - \atom only for atoms with 'rep spheres'
      - \bond only if 'rep sticks' is shown for at least one endpoint
      - cell is automatic if registered
    """
    try:
        if not file_name:
            objs = (
                cmd.get_object_list(selection) or cmd.get_object_list() or ["molecule"]
            )
            file_name = objs[0]

        path = Path(file_name)
        if path.suffix.lower() != ".tex":
            path = path.with_suffix(".tex")

        mol = TikzMolecule2D(selection)
        mol.load()
        mol.write(path)

        print(f"save_tikz: wrote {path}")
    except Exception as e:
        import traceback

        print(f"save_tikz: ERROR: {e}")
        traceback.print_exc()
