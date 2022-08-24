import mdtraj as md
import nglview as nv
from openmm import app
import openmm as mm
from simtk import unit
from sys import stdout

# Load the PDB file into an object
pdb = app.PDBFile('trpcage.pdb')

# Load the force field file into an object
forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')

# Create system object using information in the force field:
# forcefield: contains parameters of interactions
# topology: lists of atoms, residues, and bonds
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True,
    ewaldErrorTolerance=0.0005)

# Create a Langevin integrator for temperature control
integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds,
    2.0*unit.femtoseconds)

# Add a Monte Carlo barostat to the system for pressure control
system.addForce(mm.MonteCarloBarostat(1*unit.atmospheres, 300*unit.kelvin, 25))

# Use the CPU platform
platform = mm.Platform.getPlatformByName('CPU')

## Add forces between sytem setup and simulation

frc = mm.CustomBondForce("0.5*k*(r-r0)^2")
frc.addGlobalParameter("k", 10000)
frc.addGlobalParameter("r0", 2)

frc.addBond(56, 200)
system.addForce(frc)

# Create a Simulation object using definitions from above
simulation = app.Simulation(pdb.topology, system, integrator, platform)

# Set positions in the Simulation object
simulation.context.setPositions(pdb.positions)

# Minimize the energy of the system
print('Minimizing...')
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
simulation.minimizeEnergy(maxIterations=20, tolerance=100)
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

# Initialize the random velocities of the system from a Maxwell-Boltzmann distribution
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)

# Create .dcd traj file
simulation.reporters.append(app.DCDReporter('trajectory-spr.dcd', 100))

# Report information as simulation runs
simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True,
    potentialEnergy=True, temperature=True, density=True, progress=True,
    remainingTime=True, speed=True, totalSteps=10000, separator='\t'))

# Run simulation
print('Running Production...')
simulation.step(10000)
print('Done!')

traj = md.load('trajectory-spr.dcd', top='trpcage.pdb')
viewsim = nv.show_mdtraj(traj)
viewsim

import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline

traj = md.load('trajectory-spr.dcd', top='trpcage.pdb')

bond_indices = [56, 200]
Calength = md.compute_distances(traj, [bond_indices])


plt.plot(10*Calength, color='Tomato')
# note above that we have multiplied NClength by 10 to convert from nm to Å
plt.title(r'C$\alpha$-C$\alpha$ bond length')
plt.xlabel('Time (fs)')
plt.ylabel(r'Bond length ($\AA$)')
plt.show()

#
