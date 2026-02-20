################################
# AX2.pml
# 
# Gabriel Braun, 2026
################################

load AX2.mol
hide everything

color tw_amber_500,  elem A
color tw_indigo_500, elem X
alter (name A), vdw=1.8
alter (name X), vdw=1.6

show_bas

orient

select A, (AX2 and id 1)
set_view A

save_png
save_tikz
