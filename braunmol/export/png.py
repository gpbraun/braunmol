"""
png.py

Gabriel Braun, 2026
"""

import subprocess
from pathlib import Path

from pymol import cmd

from .view import SIZE_CM


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
