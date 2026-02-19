"""
png.py

Gabriel Braun, 2026
"""

import numpy as np
from pymol import cmd, stored

SIZE_CM = 10.0


@cmd.extend
def set_view(
    selection: str | None = None,
    state: int = 1,
    cm_per_a: float = 0.6,
    buffer: float = 5.0,
    _self=cmd,
):
    """
    Define distância fixa da câmera para gerar PNG e TikZ com escala consistente.

    Stores only:
      stored.tikz_angstrom_cm  -> value to use in \\begin{molecule}[angstrom=...cm]
    """
    selection = selection or "all"

    ortho_prev = _self.get_setting_int("orthoscopic")
    _self.set("orthoscopic", 0)

    _self.orient(selection, state=state)
    _self.center(selection, state=state)
    _self.origin(selection, state=state)

    fov_deg = 20
    fov_rad = np.radians(fov_deg)

    zdist = (SIZE_CM / (2.0 * cm_per_a)) / np.tan(fov_rad / 2.0)

    view = list(_self.get_view())
    view[11] = -zdist
    _self.set_view(view)
    _self.set("field_of_view", fov_deg)
    _self.clip("atoms", buffer, selection, state)

    # Compute after clip(), since it can change the front plane distance (view[15]).
    v2 = np.asarray(_self.get_view(), dtype=float)
    front = float(v2[15])

    stored.tikz_angstrom_cm = SIZE_CM / (2.0 * front * np.tan(fov_rad / 2.0))

    _self.set("orthoscopic", ortho_prev)


@cmd.extend
def align_bond(axis: str, sel1: str, sel2: str, state: int = 1, quiet: int = 1):
    """
    Rotate ONLY the camera (view) so that the vector centroid(sel1)->centroid(sel2)
    aligns with a VIEW axis (+/- x,y,z). No translation: camera distance/center unchanged.

    PyMOL convention (per set_view docs):
      - view[0:9] is a column-major (Fortran order) 3x3 rotation matrix
      - it maps model/world -> camera axes
    """
    ax = axis.strip().lower()
    target = {
        "x": np.array([1.0, 0.0, 0.0]),
        "y": np.array([0.0, 1.0, 0.0]),
        "z": np.array([0.0, 0.0, 1.0]),
        "+x": np.array([1.0, 0.0, 0.0]),
        "+y": np.array([0.0, 1.0, 0.0]),
        "+z": np.array([0.0, 0.0, 1.0]),
        "-x": np.array([-1.0, 0.0, 0.0]),
        "-y": np.array([0.0, -1.0, 0.0]),
        "-z": np.array([0.0, 0.0, -1.0]),
    }.get(ax)
    if target is None:
        raise ValueError("axis must be one of: x,y,z,+x,+y,+z,-x,-y,-z")

    c1 = cmd.get_coords(sel1, state=state)
    c2 = cmd.get_coords(sel2, state=state)
    if c1 is None or len(c1) == 0:
        raise ValueError(f"Selection '{sel1}' has no atoms (state={state}).")
    if c2 is None or len(c2) == 0:
        raise ValueError(f"Selection '{sel2}' has no atoms (state={state}).")

    p1 = c1.mean(axis=0)
    p2 = c2.mean(axis=0)

    v = p2 - p1
    vn = np.linalg.norm(v)
    if vn < 1e-12:
        raise ValueError("The two centroids are identical (zero-length vector).")
    v /= vn  # unit in world/model coordinates

    view = np.asarray(cmd.get_view(), dtype=float)

    # Correct PyMOL convention: column-major, world/model -> camera
    R = view[:9].reshape(3, 3, order="F")

    # Bond vector in camera coordinates
    w = R @ v
    wn = np.linalg.norm(w)
    if wn < 1e-12:
        raise ValueError("Numerical issue: camera-space vector length ~ 0.")
    w /= wn

    # Construct minimal camera-space rotation S such that S @ w = target
    dot = float(np.clip(w @ target, -1.0, 1.0))
    cr = np.cross(w, target)
    s = np.linalg.norm(cr)  # = sin(theta) for unit vectors
    I = np.eye(3, dtype=float)

    if s < 1e-12:
        if dot > 0.0:
            S = I
        else:
            # 180° around any axis perpendicular to w: S = -I + 2uu^T
            ref = (
                np.array([1.0, 0.0, 0.0])
                if abs(w[0]) < 0.9
                else np.array([0.0, 1.0, 0.0])
            )
            u = np.cross(w, ref)
            u /= np.linalg.norm(u)
            S = -I + 2.0 * np.outer(u, u)
    else:
        u = cr / s  # unit axis
        ux, uy, uz = u
        K = np.array([[0.0, -uz, uy], [uz, 0.0, -ux], [-uy, ux, 0.0]], dtype=float)
        # Rodrigues: cos=dot, sin=s, since theta = atan2(s, dot)
        S = dot * I + (1.0 - dot) * np.outer(u, u) + s * K

    # Rotate the camera: left-multiply in camera space. No translation update.
    Rnew = S @ R
    view[:9] = Rnew.reshape(-1, order="F")
    cmd.set_view(view.tolist())

    if not int(quiet):
        # verify
        view2 = np.asarray(cmd.get_view(), dtype=float)
        R2 = view2[:9].reshape(3, 3, order="F")
        w2 = R2 @ v
        w2 /= np.linalg.norm(w2)
        score = float(np.clip(w2 @ target, -1.0, 1.0))
        err = float(np.degrees(np.arccos(score)))
        print(f"align_bond: axis={ax}, score={score:.6f}, err≈{err:.4f}°")
