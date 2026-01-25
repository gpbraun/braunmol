################################
# nacl.pml
# 
# Gabriel Braun, 2026
################################

load_cell nacl.cif
hide everything

show_ion sticks=1
add_cell

set orthoscopic, 1
set_view
turn y, 90
turn y, -20
turn x, 15

save_png
