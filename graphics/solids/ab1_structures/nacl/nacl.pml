################################
#  nacl.pml
#
#  Gabriel Braun, 2026
################################

load_cell nacl.cif
hide everything

show_ion
add_cell
# add_vdw

set orthoscopic, 0
set_view
# turn x, 45
turn y, 90
turn y, -20
turn x, 15

save_png
save_tikz
