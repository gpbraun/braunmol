################################
# SF6.pml
# 
# Gabriel Braun, 2025
################################

load caffeine.sdf
hide everything

ball_and_stick

# turn y, 45
# turn x, 25
# move z, 2.6

save_png crop=1

set orthoscopic, 1

save_as_tikz
