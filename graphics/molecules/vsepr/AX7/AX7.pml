################################
# AX7.pml
# 
# Gabriel Braun, 2026
################################

load AX7.mol
hide everything

color tw_amber_500,  elem A
color tw_indigo_500, elem X
alter (name A), vdw=1.8
alter (name X), vdw=1.6

show_bas

orient
turn z, 90
turn x, 90
turn x, 25
set_view

save_png
save_tikz
