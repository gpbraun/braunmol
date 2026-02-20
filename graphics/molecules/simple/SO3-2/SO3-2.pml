################################
# SO3-2.pml
# 
# Gabriel Braun, 2025
################################

load SO3-2.mol
hide everything

show_bas

select S,  (SO3-2 and id 1)
select O1, (SO3-2 and id 2)
select O2, (SO3-2 and id 3)
select O3, (SO3-2 and id 4)
align xz, O1, O2, O3
turn z, 180
turn x, 25
set_view

save_png
save_tikz
