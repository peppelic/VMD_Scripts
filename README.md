# VMD_Scripts
Contains several useful VMD scripts to analyze trajectories, mainly obtained from NAMD simulations.

**CheckRingPiercing.tcl**\
Check if there is a ring piercing in the current frame loaded in VMD.
A PSF file is needed to run this script, as the bond information are 
extracted from there. Be aware that the script is based solely on
bond distances, i.e. if a bond is longer than 2 Angstrom it will be
considered as a ring piercing bond. Some bonds may natually be longer
than this cutoff.

**StaticDielectricConstant.tcl**\
Determine the dielectric constant of a solvent from a trajectory.
Run this script from the TK console with the trajectory already
loaded in VMD. The system should contain only the solvent and be 
equilibrated, otherwise the calculation may not give an accurate 
value for the dielectric constant.
