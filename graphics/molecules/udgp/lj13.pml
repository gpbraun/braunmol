################################
# udgp.pml
# 
# Gabriel Braun, 2025
################################

run ../../_pymol/pymol_config.py

load lj13.xyz
hide everything

ball_and_stick
turn y, 40

save_png
save_as_tikz
