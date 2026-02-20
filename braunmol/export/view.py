"""
view.py

Gabriel Braun, 2026
"""

import numpy as np
from pymol import cmd, stored

SIZE_CM = 10.0


def get_centroid(
    selection: str,
    state: int = 1,
    _self=cmd,
) -> np.ndarray:
    """
    Retorna: centro geométrico de uma seleção.
    """
    m = _self.get_model(selection, state=state)

    coords = np.array([a.coord for a in m.atom], dtype=float)
    centroid = coords.mean(axis=0)
    return centroid


@cmd.extend
def set_view(
    selection=None,
    state=1,
    cm_per_a=0.75,
    buffer=10.0,
    fov=20.0,
    _self=cmd,
):
    """
    Define: posicionamento determinístico da câmera.

    >>> PyMOL> set_view
    """
    selection = selection or "all"

    ortho = _self.get_setting_int("orthoscopic")
    _self.set("orthoscopic", 0)

    _self.center(selection, state=state)
    _self.origin(selection, state=state)

    _self.set("field_of_view", float(fov))
    fov_rad = np.deg2rad(float(fov))

    # Fixed camera-to-center distance (Å)
    zdist = (SIZE_CM / (2.0 * float(cm_per_a))) / np.tan(fov_rad / 2.0)

    v = np.asarray(_self.get_view(), dtype=float)
    # Deterministic camera distance.
    v[9] = 0.0
    v[10] = 0.0
    v[11] = -zdist
    # Set slab size.
    v[15] = zdist - float(buffer)
    v[16] = zdist + float(buffer)

    _self.set_view(v.tolist())

    stored.tikz_angstrom_cm = float(cm_per_a)

    _self.set("orthoscopic", ortho)


@cmd.extend
def align(
    plane: str,
    selection_1: str,
    selection_2: str,
    selection_3: str,
    state: int = 1,
    _self=cmd,
):
    """
    plane: xy, yx, xz, zx, yz, zy

    Uses the centroids of sele1/2/3. Rotates the CAMERA (view rotation only) so that:
      - the three centroids lie in the given camera plane
      - the line (sele1 -> sele2) aligns with the FIRST axis of the plane
      - sele3 fixes the in-plane roll so it lies on the +SECOND axis side

    >>> PyMOL> align
    """
    p = plane.strip().lower()
    a, b = p[0], p[1]
    c = ({"x", "y", "z"} - {a, b}).pop()

    # centroids
    p1 = get_centroid(selection_1, state, _self)
    p2 = get_centroid(selection_2, state, _self)
    p3 = get_centroid(selection_3, state, _self)

    # unit helpers (inline)
    v = p2 - p1
    v /= np.linalg.norm(v)

    n = np.cross(p2 - p1, p3 - p1)
    n /= np.linalg.norm(n)

    w = np.cross(n, v)
    w /= np.linalg.norm(w)

    # orient sign so that (p3 - p1) has positive projection on +b axis
    if np.dot(p3 - p1, w) < 0.0:
        n = -n
        w = -w

    # build camera axes in model coordinates per requested plane ordering
    # a-axis -> v, b-axis -> w, c-axis -> n
    axis = {a: v, b: w, c: n}
    ex, ey, ez = axis["x"], axis["y"], axis["z"]

    # two plausible model->cam rotations; choose best for v->+a and n->+c
    M1 = np.vstack([ex, ey, ez])
    M2 = M1.T

    idx = {"x": 0, "y": 1, "z": 2}
    ea = np.zeros(3)
    ea[idx[a]] = 1.0
    ec = np.zeros(3)
    ec[idx[c]] = 1.0

    s1 = (M1 @ v).dot(ea) + (M1 @ n).dot(ec)
    s2 = (M2 @ v).dot(ea) + (M2 @ n).dot(ec)
    M = M1 if s1 >= s2 else M2

    view = list(_self.get_view())
    view[0:9] = M.reshape(9, order="F").tolist()
    _self.set_view(view)
