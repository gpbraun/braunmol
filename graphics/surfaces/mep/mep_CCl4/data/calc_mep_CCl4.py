from pathlib import Path

import psi4

XYZFILE = Path("geom.xyz")

mol = psi4.geometry(XYZFILE.read_text())

psi4.set_options(
    {
        "basis": "cc-pvdz",
        "df_basis_scf": "cc-pvdz-jkfit",
        "cubeprop_tasks": ["esp", "density"],
    }
)

energy, wfn = psi4.energy("scf", return_wfn=True)

psi4.cubeprop(wfn)

print("ESP and density cube files written successfully.")
