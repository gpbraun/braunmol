run ../../pymol/pymol_config.py

bg_color black
set orthoscopic, 1

load_cell KH8O4F.cif, reps=1-1-1
hide everything

show_ionic sticks=1
# add_vdw
add_cell

turn y, -20
turn x, 15

zoom complete=1
# move z, -10
clip slab, 20

view home, store

# teste zns.png