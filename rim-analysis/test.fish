#!/usr/bin/env fish

for i in (seq 1 50)
    set DIR ~/Documents/lammps/results/0K_8keV_angle0/run_$i
    echo $i
    RUST_LOG=info ../target/release/rim-analysis $DIR/dump.initial $DIR/dump.final_no_cluster $DIR -c 3 2>&1
end
