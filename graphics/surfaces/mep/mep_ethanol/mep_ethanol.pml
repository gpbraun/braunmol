################################
# mep_ethanol.pml
# 
# Gabriel Braun, 2026
################################

load data/geom.xyz, ethanol
load data/Da.cube,  __den_cube
load data/ESP.cube, __mep_cube
hide everything

show_bas
mep_surface __den_cube, __mep_cube

set_view
turn x, 180
turn z, -35
turn y, 30

save_png mep_ethanol.png
