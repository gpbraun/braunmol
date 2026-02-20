################################
# AX2E.pml
# 
# Gabriel Braun, 2026
################################

load AX2E.mol
hide everything

color tw_amber_500,  elem A
color tw_indigo_500, elem X
alter (name A), vdw=1.8
alter (name X), vdw=1.6

show_bas

select A,  (AX2E and id 1)
select X1, (AX2E and id 2)
select X2, (AX2E and id 3)
align xy, X1, X2, A

set_view A

save_png
save_tikz
