################################
# CH4.pml
# 
# Gabriel Braun, 2025
################################

load CH4.sdf
hide everything

show_bas

select C,  (CH4 and id 2)
select H1, (CH4 and id 1)
select H2, (CH4 and id 3)
align yz, C, H1, H2
turn x, 25
set_view

save_png
save_tikz
