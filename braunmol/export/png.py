"""
png.py

Gabriel Braun, 2026
"""

import subprocess
from pathlib import Path

import numpy as np
from pymol import cmd

SIZE_CM = 10.0


@cmd.extend
def set_view(
    selection: str | None = None,
    state: int = 1,
    cm_per_A: float = 1.0,
    buffer: float = 5.0,
    _self=cmd,
):
    """
    Define distância fixa da câmera.

    Geração sistemática de imagens com tamanho consistente.

    >>> PyMOL> set_camera [ selection [, state [, cm_per_A [, buffer ]]]]
    """
    selection = selection or "all"

    ortho_prev = _self.get_setting_int("orthoscopic")
    _self.set("orthoscopic", 0)

    # orientação
    _self.orient(selection, state=state)
    _self.center(selection, state=state)
    _self.origin(selection, state=state)

    # calcula zdist
    fov_deg = _self.get_setting_float("field_of_view")
    fov_rad = np.radians(fov_deg)

    # muda o view
    zdist = (SIZE_CM / (2.0 * cm_per_A)) / np.tan(fov_rad / 2.0)

    view = list(_self.get_view())
    view[11] = -float(zdist)
    _self.set_view(view)
    _self.clip("atoms", float(buffer), selection, state)

    # restora orthoscopic
    _self.set("orthoscopic", ortho_prev)


@cmd.extend
def save_png(
    file_name: str | None = None,
    dpi: int = 300,
    ray: bool = True,
    crop: bool = True,
    _self=cmd,
):
    """
    Salva a imagem.

    >>> PyMOL> save_png [ file_name [, dpi [, ray [, crop ]]]]
    """
    file_name = file_name or _self.get_object_list()[0]
    file_path = Path(file_name).with_suffix(".png")

    _self.png(
        str(file_path),
        width=f"{SIZE_CM}cm",
        height=f"{SIZE_CM}cm",
        dpi=dpi,
        ray=ray,
    )

    if int(crop):
        subprocess.run(["convert", str(file_path), "-trim", "+repage", str(file_path)])

    file_size = round(file_path.stat().st_size / 1024)

    print(f"Arquivo salvo em {file_path.resolve()} ({file_size} kB)")
