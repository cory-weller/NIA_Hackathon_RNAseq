#!/usr/bin/env bash

hard_code() {
    sbatch STAR_align.sh Foliate 1
    sbatch STAR_align.sh Foliate 2
    sbatch STAR_align.sh Foliate 3
    sbatch STAR_align.sh Foliate 4
    sbatch STAR_align.sh Foliate 5
    sbatch STAR_align.sh Organoid 1
    sbatch STAR_align.sh Organoid 2
    sbatch STAR_align.sh Organoid 3
    sbatch STAR_align.sh Organoid 4
    sbatch STAR_align.sh Organoid 5
}

for_loop() {
    for treatment in Foliate Organoid; do
        for sample in 1 2 3 4 5; do
            sbatch STAR_align.sh ${treatment}_${sample}
        done
    done
}

gnu_parallel() {
    parallel -j 1 sbatch STAR_align.sh ::: Foliate Organoid ::: $(seq 1 5)
}