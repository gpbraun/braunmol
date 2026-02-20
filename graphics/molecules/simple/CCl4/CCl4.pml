################################
# CCl4.pml
# 
# Gabriel Braun, 2025
################################

load CCl4.mol
hide everything

show_bas

select C,   (CCl4 and id 1)
select Cl1, (CCl4 and id 2)
select Cl2, (CCl4 and id 3)
align yz, C, Cl1, Cl2
turn x, 25
set_view

save_png
save_tikz
