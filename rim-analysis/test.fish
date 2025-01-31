#!/usr/bin/env fish

for i in (seq 1 50)
    set DIR ~/Documents/lammps/2024_dec_paper/results/0K_8keV_angle0/run_$i
    RUST_LOG=info ../target/release/rim-analysis $DIR/dump.initial $DIR/dump.final_no_cluster $DIR -c 3
end
