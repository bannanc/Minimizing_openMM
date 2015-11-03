# Minimizing_openMM
Script to minimize using openmm intended as a way around the GROMACS minimizer

This is a rough draft for a script that would take in GROMACS .gro and .top files and minimize the system. It would return a .gro file of the minimized configuration. It uses openmm (from openmoltools) and ParmEd both are available here on gitHub. 

If this is a fix we continue to use I will make it more general so it can be customized from the command line. Since it uses ParmEd it has the potential to be made more general than just GROMACS input files, but right now that is the only immediate use I have for it.

