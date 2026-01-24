################################
# SF6.pml
# 
# Gabriel Braun, 2025
################################

load SF6.sdf
hide everything

set_zdist_cm
ball_and_stick

turn y, 45
turn x, 25
# move z, 2.6

save_png crop=0

# save_as_tikz
