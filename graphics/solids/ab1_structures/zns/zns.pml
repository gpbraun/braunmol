################################
# zns.pml
# 
# Gabriel Braun, 2026
################################

set orthoscopic, 1

load_cell zns.cif
hide everything

show_ionic sticks=1
add_cell

set_view
turn y, 90
turn y, -20
turn x, 15

save_png
