# Mask file for LHC and HL-LHC

In this repository you can retrieve and contribute to improve the MADX masks for LHC and HL-LHC.

This work is based on the work of many colleagues  who shared with the community the efforts spent for enriching this simulation framework.

The generic mask file has two main parts:
 1. the definition of the configuration parameters. Their value can be assigned explicitly or masked by a placeholder (to be used for SixTrack scans).
 2. the call of the madx files 
    - to load the sequence and optics without beam-beam, to define the beam crossing angle and separation, the status of the experimental magnet...
    - to level the luminosity and install/configure the BB lenses...
    - to load the beam to track and install the magnetic errors, power the octupoles, match the tunes and the chromaticities...
    - to make the final twiss and prepare the input files for SixTrack.

The separation of the configuration parameters of the mask with the MADX code aims 
- to improve the trasparency and readibility for the users,
- to define better interfaces for the maintenance of the MADX code.

## Contributors:
<!-- use two spaces for new line -->
R. De Maria  
S. Fartoukh  
M. Giovannozzi  
G. Iadarola  
Y. Papaphilippou  
D. Pellegrini  
G. Sterbini  
F. Van Der Veken  

