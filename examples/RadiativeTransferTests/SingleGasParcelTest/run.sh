#!/bin/bash

# make run.sh fail if a subcommand fails
set -e

if [ ! -f ./single_gas_parcel_test.hdf5 ]; then
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
    ./single_gas_parcel.yml 2>&1 | tee output.log

