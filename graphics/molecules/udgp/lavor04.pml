################################
# udgp.pml
# 
# Gabriel Braun, 2025
################################

run ../../_pymol/pymol_config.py

load lavor04.xyz
hide everything

ball_and_stick

turn z, 70
turn x, 30
turn y, 180
zoom lavor04

save_png
save_as_tikz
