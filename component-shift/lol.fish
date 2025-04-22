#!/usr/bin/env fish

set base $HOME/Documents/lammps/2024_dec_paper/results

set angles 0 10 20 30 40 50 60 65 70 75 80
set table shift_components.txt

for angle in $angles
    set dir $base/0K_8keV_angle$angle
    echo $dir
    echo -n '' >$dir/$table
    for run in (seq 1 100)
        ../target/release/component-shift --depth 0.0 \
            $dir/run_$run/dump.initial \
            $dir/run_$run/dump.final_no_cluster >>$dir/$table
    end
    cp $dir/$table $base/shift_components/{$angle}_$table
end
