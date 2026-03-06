################################
# NH3.pml
# 
# Gabriel Braun, 2025
################################

load NH3.mol
hide everything

show_bas

select N,  (NH3 and id 1)
select H1, (NH3 and id 2)
select H2, (NH3 and id 3)
select H3, (NH3 and id 4)
align xz, H1, H2, H3
turn z, 180
turn x, 25
set_view

save_png 
save_tikz
