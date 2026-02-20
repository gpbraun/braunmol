################################
# NO2.pml
# 
# Gabriel Braun, 2025
################################

load NO2.mol
hide everything

show_bas

orient
select N,  (NO2 and id 1)
select O1, (NO2 and id 2)
select O2, (NO2 and id 3)
align xy, O1, O2, N
set_view N

save_png
save_tikz
