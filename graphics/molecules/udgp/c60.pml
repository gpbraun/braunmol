################################
# udgp.pml
# 
# Gabriel Braun, 2025
################################

run ../../_pymol/pymol_config.py

load c60.xyz
hide everything

ball_and_stick

zoom c60
turn x, 10

save_png
save_as_tikz
