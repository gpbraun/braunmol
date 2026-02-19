################################
# CCl4.pml
# 
# Gabriel Braun, 2025
################################

load CCl4.mol
hide everything

show_bas

set_view
select a1, (CCl4 and id 1)
select a2, (CCl4 and id 2)
align_bond y, a1, a2
turn y, 90
turn x, 20

save_png
save_tikz
