#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

if [ ! -f ./ionization_balance_test.hdf5 ]; then
    echo "creating ICs"
    python3 makeIC.py
fi

# Run SWIFT with RT
../../swift \
    --hydro \
    --threads=1 \
    --verbose=0  \
    --radiation \
    --stars \
    --feedback \
    --external-gravity \
    ./ion_balance.yml 2>&1 | tee output.log

