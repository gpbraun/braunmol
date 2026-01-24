import psi4

MOLECULE = "imidazole"

mol = psi4.geometry(f"pubchem:{MOLECULE}")
mol.save_xyz_file("geom.xyz", 1)

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
