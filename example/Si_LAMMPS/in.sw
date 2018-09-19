units           metal
atom_style      atomic
boundary        p p p

read_data       tmp.lammps

pair_style      sw
pair_coeff 	* * Si.sw Si

dump            1 all custom 1 FORCE id fx fy fz
dump            2 all custom 1 COORD id xu yu zu
dump_modify     1 format float "%20.15f"
dump_modify     2 format float "%20.15f"
run             0

