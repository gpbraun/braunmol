################################
# AX4.pml
# 
# Gabriel Braun, 2026
################################

load AX4.mol
hide everything

color tw_amber_500,  elem A
color tw_indigo_500, elem X
alter (name A), vdw=1.8
alter (name X), vdw=1.6

show_bas

select A,  (AX4 and id 1)
select X1, (AX4 and id 2)
select X2, (AX4 and id 3)
align yz, A, X1, X2
turn x, 25

set_view A

save_png
save_tikz
