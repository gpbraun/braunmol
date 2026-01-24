################################
# tio2.pml
# 
# Gabriel Braun, 2026
################################

set orthoscopic, 1

load_cell tio2.cif
hide everything

show_ionic sticks=1
add_cell

turn x, 90
turn y, -20
turn x, 15

set_camera
save_png
