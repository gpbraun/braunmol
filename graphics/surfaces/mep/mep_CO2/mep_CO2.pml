################################
# mep_CO2.pml
# 
# Gabriel Braun, 2026
################################

load data/geom.xyz, CO2
load data/Da.cube,  __den_cube
load data/ESP.cube, __mep_cube
hide everything

show_bas
mep_surface __den_cube, __mep_cube, low=-0.20, high=0.60

set_view

save_png mep_CO2.png
