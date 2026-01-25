################################
# phenol.pml
# 
# Gabriel Braun, 2025
################################

load data/geom.xyz, imidazole
load data/Da.cube,  __den_cube
load data/ESP.cube, __mep_cube
hide everything

show_bas
mep_surface __den_cube, __mep_cube

set_view

save_png mep_imidazole.png, crop=1
