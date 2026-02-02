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
    cm_per_a: float = 0.8,
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

    fov_deg = float(_self.get_setting_float("field_of_view"))
    fov_rad = np.radians(fov_deg)

    zdist = (SIZE_CM / (2.0 * cm_per_a)) / np.tan(fov_rad / 2.0)

    view = list(_self.get_view())
    view[11] = -float(zdist)
    _self.set_view(view)
    _self.clip("atoms", float(buffer), selection, state)

    # Compute after clip(), since it can change the front plane distance (view[15]).
    v2 = np.asarray(_self.get_view(), dtype=float)
    front = float(v2[15])

    stored.tikz_angstrom_cm = SIZE_CM / (2.0 * front * np.tan(fov_rad / 2.0))

    _self.set("orthoscopic", ortho_prev)
