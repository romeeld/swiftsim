#!/bin/bash
set -o xtrace

# Resolution
N_label=n70

# Create the initial particle planets
python3 make_init_cond.py

# Download the equation of state tables if not already present
if [ ! -e ../EoSTables/CD21_HHe.txt ]
then
    cd ../EoSTables
    ./get_eos_tables.sh
    cd -
fi

# Run SWIFT settling simulations
../../../swift --hydro --self-gravity --threads=28 demo_target_"$N_label".yml \
    2>&1 | tee output_"$N_label"_t.txt

# Plot the settled particles
python3 plot_snapshots.py
python3 plot_profiles.py
