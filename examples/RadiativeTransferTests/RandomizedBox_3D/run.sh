#!/bin/bash

# make run.sh fail if a subcommand fails
set -e
set -o pipefail

if [ ! -f 'randomized-sine.hdf5' ]; then
    echo "Generating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
# ../../swift \
gdb -ex run --args ../../swift \
    --hydro \
    --threads=13 \
    --verbose=0  \
    --radiation \
    --self-gravity \
    --stars \
    --feedback \
    --steps=3000 \
    -e \
    ./randomized-rt.yml 2>&1

# echo "running sanity checks"
# python3 ../UniformBox_3D/rt_sanity_checks.py | tee sanity_check.log
