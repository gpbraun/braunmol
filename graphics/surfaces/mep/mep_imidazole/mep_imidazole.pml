################################
# phenol.pml
# 
# Gabriel Braun, 2025
################################

load data/geom.xyz, imidazole
load data/Da.cube,  __den_cube
load data/ESP.cube, __mep_cube

hide everything

ball_and_stick
mep_surface __den_cube, __mep_cube

save_png mep_imidazole.png, crop=1
