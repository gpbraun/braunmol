################################
# tio2.pml
# 
# Gabriel Braun, 2026
################################

set orthoscopic, 1

load_cell tio2.cif
hide everything

show_ion sticks=1
add_cell

set_view
turn x, 90
turn y, -20
turn x, 15

save_png
