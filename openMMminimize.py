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

def run_minimize(top, system, integrator, tolerance, steps, writeSteps = 0, pdbFile = "output.pdb"):
    """
    Runs a openmm minimization on an input system and topology. 
    input:
        top = ParmEd or openmm object with properties topology and positions
        system = ParmEd or openmm system created from top
        integrator = openmm integrator with assigned parameters
        tolerance = float, energy tolerance in kJ/mol
        steps = integer, number of steps the minimization should take
        writeSteps = integer (0), frequency to write coordinates to pdb file
        pdbFile = string ('output.pdb'), file to save steps
    
    output:
        returns openmm simulation in the conditions after it has been run. Outputs pdb file of trajectory if writeSteps > 0.
    """
    # Build simulation
    print "Assigning simulation parameters..."
    simulation = app.Simulation(top.topology, system, integrator)
    print "Setting positions in simulation..."
    simulation.context.setPositions(top.positions)
    print "Assigning minimization to simulation ..."
    simulation.minimizeEnergy(tolerance = tolerance)
    if writeSteps > 0:
        print "Assigning conditions for writing trajectory..."
        simulation.reporters.append(app.PDBReporter(pdbFile, writeSteps))

    # Running minimization
    print "OpenMM Minimization is running..."
    simulation.step(steps)
    return simulation 

def load_and_minimize(topFile, groFile, steps,
        output = "output.gro",
        temperature = 300.0, 
        friction = 0.3, 
        timestep = 2.0,
        cutOff = 1.2,
        writeSteps = 0,
        tolerance = 2.0,
        pdbFile = "output.pdb"):

    """
    Performs a minimization using openmm on GROMACS input files and returns the ParmEd system with the input topology and the minimized coordinates. 

    input:
        topFile = string, GROMACS .top topology file
        groFile = string, GROMACS .gro coordinate file
        steps = integer, number of steps the minimization should take
        output = ("output.gro") string, GROMACS .gro coordinate file
        temperature = (300.) float, temperature in kelvin
        friction = (0.3) float, friction in inverse picoseconds
        timestep = (2.0) float, size of each step in femtoseconds
        cutOff = (1.2) float, nonbondedCutoff in nanometers
        writeSteps = (0) integer, if non-zero write coordinates to file
        tolerance = (2.0) float, energy tolerance in kJ/mol 
        pdbFile = ("output.pdb") string, where to writeSteps

    output:
        returns top, openmm/parmEd system with coordinates after minimization
        outputs GROMACS coordinate file with final positions
        optionally outputs pdb file with trajectory of minimization
    """
    print "Beginning minimization with openMM: "

    # Load .gro and .top files into parmed
    print "Loading %s and %s into ParmEd..." % (topFile, groFile)
    top = parmed.load_file(topFile)
    gro = parmed.load_file(groFile)
    print "Assign gro box to the topology..."
    top.box = gro.box
    top.positions = gro.positions

    # Make openmm system
    print "Creating system in openMM..."
    system = top.createSystem(
        nonbondedMethod = app.PME, 
        nonbondedCutoff = cutOff * u.nanometer, 
        constraints = app.HBonds)
    print "Creating integrator..."
    integrator = mm.LangevinIntegrator(
        temperature * u.kelvin, 
        friction / u.picosecond, 
        timestep * u.femtosecond)
    
    simulation = run_minimize(top, system, integrator, tolerance, steps, writeSteps, pdbFile)
    
    # Save Final positions
    print "Saving final positions of OpenMM Minimization to %s ... " % output
    top.positions = simulation.context.getState(getPositions=True).getPositions()
    parmed.gromacs.GromacsGroFile.write(top, output)

    print "OpenMM minimization complete"
    return top

# From the command line
if __name__ == '__main__':
    from optparse import OptionParser
    import os
    parser = OptionParser(usage = "This takes a GROMACS topology and coordinate file and runs an openmm minimization. It outputs a new GROMACS coordinate file with the resulting minimization. It can optionally output a PDB file with the trajectory of the minimization.")

    parser.add_option('-t', '--top',
            help = "REQUIRED! GROMACS topology file, must end in .top or an error will be raised.",
            dest = 'top')

    parser.add_option('-g','--gro',
            help = "REQUIRED! GROMACS coordinate file, must end in .gro or an error will be raised.",
            dest = 'gro')

    parser.add_option('-n','--steps',
            help = "REQUIRED! integer, the number of steps the minimization should take",
            type = "int",
            dest = 'steps')

    parser.add_option('-o','--output',
            default = 'output.gro',
            help = "output file, GROMACS coordinate file for minimization results",
            dest = 'output')

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
            help = "float, time for each step in femtoseconds",
            type = "float",
            dest = 'timestep')

    parser.add_option('-c', '--cutoff',
            default = 1.2,
            help = "float, distance in nanometers for nonbonded interaction cutoff",
            type = "float",
            dest = 'cutOff')

    parser.add_option('-w', '--write',
            default = 0,
            help = "integer, if non-zero frequency for writing trajectory to pdb",
            type = 'int',
            dest = 'writeSteps')

    parser.add_option('-l', '--tolerance',
            default = 2.0,
            help = "float, energy tolerance in kJ/mol",
            type = "float",
            dest = 'tolerance')

    parser.add_option('-p', '--pdb',
            help = "file name for pdb, if write is non-zero this is required!",
            dest = 'pdbFile')

    (opt, args) = parser.parse_args()

    if opt.top == None:
        parser.error("ERROR: No topology file was provided")
    elif not os.path.isfile(opt.top):
        parser.error("ERROR: topology file not found.")
    elif opt.top.split('.')[1] != 'top':
        parser.error("ERROR: topology file must be a GROMACS .top file, extension incorrect")

    if opt.gro == None:
        parser.error("ERROR: No coordinate file was provided")
    elif not os.path.isfile(opt.gro):
        parser.error("ERROR: coordinate file not found.")
    elif opt.gro.split('.')[1] != 'gro':
        parser.error("ERROR: coordinate file must be a GROMACS .gro file, extension incorrect")

    if opt.steps == None:
        parser.error("ERROR: Number of steps was not provided and is required")

    if opt.pdbFile == None and opt.writeSteps > 0:
        parser.error("ERROR: if you want to write the trajectory with frequency write, you must provide a pdb file name")

    if opt.pdbFile != None and opt.writeSteps == 0:
        print "WARNING: you provided a pdb file name, but the frequency for writing the trajectory was zero, no pdb file will be created."

    load_and_minimize(opt.top, opt.gro, opt.steps, opt.output, opt.temperature, opt.friction, opt.timestep, opt.cutOff, opt.writeSteps, opt.tolerance, opt.pdbFile)
    
