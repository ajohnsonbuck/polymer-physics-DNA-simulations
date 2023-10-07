# polymer-physics-DNA-simulations
Tools for simulation of DNA polymer dynamics, calculating local effective concentrations and rms end-to-end distances


# Extended Description

## swing_arm_simulation_pollack_multi_trials.m:

Performs polymer physics simulations of engineered DNA swinging arm constructs to determine local effective 
concentration of functional components as described in J. Fu et al. Nature Nanotechnology 9, pages 531–536 (2014).

For flexible polymer simulations (ssDNA), a freely-rotating chain model with electrostatics is used, after 
Meisburger et al. (Lois Pollack) Biopolymers, 2013, “Polyelectrolyte properties of single stranded DNA measured 
using SAXS and single molecule FRET: Beyond the wormlike chain model. For rigid polymer simulations (dsDNA), 
a rigid rod model is used.

Requires functions **swing_arm_simulation_pollack_finitewidth.m** and **FRC_pollack_ssDNA_finitewidth.m** to
be in the same directory, or on the MATLAB path.

## fjc_rmsd

Standalone function that estimates root-mean-squared end-to-end distance of a freely jointed chain polymer. See script header for instructions on use.

## wlc_rmsd

Standalone function that estimates root-mean-squared end-to-end distance of a worm-like chain polymer. See script header for instructions on use.

# License
This project is licensed under the BSD 3-Clause License.

# Author

Alex Johnson-Buck, 2014, The University of Michigan
