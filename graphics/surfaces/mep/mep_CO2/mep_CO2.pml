################################
# mep_CCl4.pml
# 
# Gabriel Braun, 2026
################################

load data/geom.xyz, CCl4
load data/Da.cube,  __den_cube
load data/ESP.cube, __mep_cube
hide everything

show_bas
mep_surface __den_cube, __mep_cube, low=0.05, high=0.35

set_view
select a1, (CCl4 and id 1)
select a2, (CCl4 and id 2)
align_bond y, a1, a2
turn y, 90
turn x, 20

save_png mep_CCl4.png
