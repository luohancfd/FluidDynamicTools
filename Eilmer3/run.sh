#!/bin/bash

# 1. Prepare gas file
gasfile o2-2T.inp o2-2T.lua
# We want harmonic oscillators instead of truncated harmonic
# oscillators for this literature comparison.
#sed -ie 's/truncated harmonic/harmonic/g' o2-2T.lua
sed -ie 's/truncated harmonic/anharmonic/g' o2-2T.lua
mv o2-2T.lua o2-2T-HO.lua

# 2. Run simulation
poshax3.x o2shock.cfg

python poshax2tec.py 4400
# 3. Plot results
#gnuplot profiles.gplot

# 4. Convert to png files
#convert -density 600x600 -quality 90 number_density_profiles.eps -scale 750 number_density_profiles.png
#convert -density 600x600 -quality 90 temperature_profiles.eps -scale 750 temperature_profiles.png 

# 5. Compute average errors with Panesi solution
#python compute_errors.py 
