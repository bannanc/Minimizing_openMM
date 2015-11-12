import os
import sys
# Import for openmm
import simtk.openmm as mm
from simtk.openmm import app
# Import for parmed
import parmed
from parmed import unit as u

"""
Written to work around the GROMACS minimizer that is getting consistently stuck. It uses methods from ParmEd and OpenMM to read in GROMACS topology and coordinate files. Do a minimization and then writes out a coordinate file for the minimized system. 

Modules needed:
    openMM - available at http://github.com/pandegroup/openmm
    ParmEd - available at http://github.com/ParmEd/ParmEd
"""

def run_minimize(top, system, integrator, tolerance, maxIt):

    """
    Runs a openmm minimization on an input system and topology. 
    input:
        top = ParmEd or openmm object with properties topology and positions
        system = ParmEd or openmm system created from top
        integrator = openmm integrator with assigned parameters
        tolerance = float, energy tolerance in kJ/mol
        maxIt = integer, max number of iterations in the minimization, if zero it will continue until it reaches the tolerance specified.

    output:
        returns openmm simulation in the conditions after it has been run. 
    """
    # Build simulation
    print "\tAssigning simulation parameters..."
    simulation = app.Simulation(top.topology, system, integrator)
    print "\tSetting positions in simulation..."
    simulation.context.setPositions(top.positions)
    print "\tMinimizing system..." 
    simulation.minimizeEnergy(tolerance = tolerance * u.kilojoule_per_mole, maxIterations = maxIt)
    
    return simulation 

def load_and_minimize(topFile, groFile, maxIt,
        output = "output.gro",
        EnTol = 10.0,
        temperature = 300.0, 
        friction = 1.0, 
        timestep = 2.0,
        cutOff = 1.2,
        IntTol = 0.00001):

    """
    Performs a minimization using openmm on GROMACS input files and returns the ParmEd system with the input topology and the minimized coordinates. 

    input:
        topFile = string, GROMACS .top topology file
        groFile = string, GROMACS .gro coordinate file
        maxIt = integer, maximum number of steps the minimization should take. If zero it will continue to minimize until it reaches the specified tolerance
        output = ("output.gro") string, GROMACS .gro coordinate file
        EnTol = (10.0) float, energy tolerance in kJ/mol 
        temperature = (300.) float, temperature in kelvin
        friction = (1.0) float, friction in inverse picoseconds
        timestep = (2.0) float, size of each step in femtoseconds
        cutOff = (1.2) float, nonbondedCutoff in nanometers
        IntTol = (0.00001) float, fraction of a distance within which constraints are maintained in the openMM integrator. 

    output:
        returns top, openmm/parmEd system with coordinates after minimization
        outputs GROMACS coordinate file with final positions
    """
    print "Checking status of input files..."
    if not os.path.isfile(topFile):
        topologyFileError = Exception("Topology file not found please check the path was entered correctly")
        raise topologyFileError
    if not os.path.isfile(groFile):
        coordinateFileError = Exception("Coordinate file was not found, please check the path was entered correctly")
        raise coordinateFileError
    
    print "Files found! openMM setup and minimization: "

    # Load .gro and .top files into parmed
    print "\tLoading %s and %s into ParmEd..." % (topFile, groFile)
    try:
        top = parmed.load_file(topFile)
    except:
        ParmEdTopologyError = Exception("ParmEd could not load your topology file, please check the ParmEd documentation at https://parmed.github.io/ParmEd for more information")
        raise ParmEdTopologyError
    try:
        gro = parmed.load_file(groFile)
    except:
        ParmEdCoordinateError = Exception("ParmEd could not load your coordinate file, please check the ParmEd documentation at https://parmed.github.io/ParmEd for more information")
        raise ParmEdCoordinateError

    print "\tAssign coordinates to the topology..."
    top.box = gro.box
    top.positions = gro.positions

    # Make openmm system
    print "\tCreating system in openMM..."
    system = top.createSystem(
        nonbondedMethod = app.PME, 
        nonbondedCutoff = cutOff * u.nanometer, 
        constraints = app.HBonds)
    print "\tCreating integrator..."
    integrator = mm.LangevinIntegrator(
        temperature * u.kelvin, 
        friction / u.picosecond, 
        timestep * u.femtosecond)
    integrator.setConstraintTolerance(IntTol)

    print "\tAssigning platform:", platform
    simulation = run_minimize(top, system, integrator, EnTol, maxIt)
    
    # Save Final positions
    print "\tSaving final positions of OpenMM Minimization to %s ... " % output
    top.positions = simulation.context.getState(getPositions=True).getPositions()
    try:
        parmed.gromacs.GromacsGroFile.write(top, output)
    except:
        print "ParmEd could not could not write the file coordinates to %s, this module writes out to a GROMACS .gro coordinate file, please check the ParmEd documentation at https://parmed.github.io/ParmEd for more information" % output

    print "Returning the miminimized system as ParmEd system: ",top
    return top

# From the command line
if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage = "This takes a GROMACS topology and coordinate file and runs an openmm minimization. It outputs a new GROMACS coordinate file with the resulting minimization. It can optionally output a PDB file with the trajectory of the minimization.")

    parser.add_option('-t', '--top',
            help = "REQUIRED! GROMACS topology file, must end in .top or an error will be raised.",
            dest = 'top')

    parser.add_option('-g','--gro',
            help = "REQUIRED! GROMACS coordinate file, must end in .gro or an error will be raised.",
            dest = 'gro')

    parser.add_option('-n','--maxIt',
            help = "REQUIRED! integer, max number of iterations in the minimization, if zero it will continue until it reaches the tolerance specified.",
            type = "int",
            dest = 'maxIt')

    parser.add_option('-o','--output',
            default = 'output.gro',
            help = "output file, GROMACS coordinate file for minimization results",
            dest = 'output')

    parser.add_option('-l', '--EnTol',
            default = 10.0,
            help = "float, energy tolerance in kJ/mol",
            type = "float",
            dest = 'EnTol')

    parser.add_option('-k', '--temp',
            default = 300.0,
            help = "float, temperature in kelvin",
            type = "float",
            dest = 'temperature')

    parser.add_option('-f', '--frict',
            default = 0.3,
            help = "float, friction in inverse picoseconds",
            type = "float",
            dest = 'friction')

    parser.add_option('-s','--timestep',
            default = 2.0,
            help = "float, time for each step in femtoseconds for integrator",
            type = "float",
            dest = 'timestep')

    parser.add_option('-c', '--cutoff',
            default = 1.2,
            help = "float, distance in nanometers for nonbonded interaction cutoff",
            type = "float",
            dest = 'cutOff')

    parser.add_option('-d','--IntTol',
            default = 0.00001,
            help = "float, fraction of a distance within which constraints are maintained in the openMM integrator",
            type = "float",
            dest = 'IntTol')

    parser.add_option('-p','--platform',
            default = 'CPU',
            help = "string, openMM name for the platform you are using",
            dest = 'platform')

    (opt, args) = parser.parse_args()

    if opt.top == None:
        parser.error("ERROR: No topology file was provided")

    if opt.gro == None:
        parser.error("ERROR: No coordinate file was provided")

    if opt.maxIt == None:
        parser.error("ERROR: Max iterations was not provided and is required")

    load_and_minimize(opt.top, opt.gro, opt.maxIt, opt.output, opt.EnTol, opt.temperature, opt.friction, opt.timestep, opt.cutOff, opt.IntTol, opt.platform)
    
