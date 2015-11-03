# Import for openmm
import simtk.openmm as mm
from simtk.openmm import app
# Import for parmed
import parmed
from parmed import unit as u

print "Beginning minimization with openMM: "

# Load .gro and .top files into parmed
groFiles = 'morph'
print "Loading %s.top and %s.gro into ParmEd..." % (groFiles, groFiles)
top = parmed.load_file('%s.top' % groFiles)
gro = parmed.load_file('%s.gro' % groFiles)
top.box = gro.box

# Make openmm system
temp = 300 * u.kelvin
frict = 0.3 / u.picosecond
timestep = 2.0 * u.femtosecond
cutOff = 1.2 * u.nanometer
print "Creating system in openMM..."
system = top.createSystem(nonbondedMethod = app.PME, nonbondedCutoff = cutOff, constraints = app.HBonds)
integrator = mm.LangevinIntegrator(temp, frict, timestep)

# Build simulation
steps = 2
writeSteps = 1
tolerance = 2.0
print "Assigning simulation parameters..."
simulation = app.Simulation(top.topology, system, integrator)
simulation.context.setPositions(gro.positions)
simulation.minimizeEnergy(tolerance = tolerance)
simulation.reporters.append(app.PDBReporter('openMin.pdb',writeSteps))

# Running minimization
print "OpenMM Minimization is running..."
simulation.step(steps)

# Save Final positions
out = 'min0'
output = '%s.gro' % out
print "Saving final positions of OpenMM Minimization to %s ... " % output
top.positions = simulation.context.getState(getPositions=True).getPositions()
parmed.gromacs.GromacsGroFile.write(top, output)

print "OpenMM minimization complete"
