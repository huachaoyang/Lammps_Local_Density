#usr/bin/env python


import sys, os
sys.path.append('./')

import numpy as np

import sim

LammpsExec = "./lammpsLD/src/lmp_ZIN_test"
doPairWrite = False     # flag to run pair_localdensity.cpp with just pair-write (no MD run)
compare_values = True   # flag to compare energies by re-evaluating them with sim


#parameter for number of molecules
NMol = 2


#setpoint temperature
TempSet = 1.0


np.random.seed(12345)


#define an atom type; give them names, masses, and charge
AtomTypeA = sim.chem.AtomType("A", Mass = 1.0, Charge = 0., Color = (1,0,0))

#define a molecule type; give it a name and list of component atom types
MolType = sim.chem.MolType("M", [AtomTypeA])

#define the world in terms of a list of molecular species and dimensionality 
World = sim.chem.World([MolType], Dim = 3, Units = sim.units.DimensionlessUnits)

#make a system that exists in this world
SysName = "testmdlocaldensity"
Sys = sim.system.System(World, Name = SysName)

#add instances of the molecule type to the system
for i in range(NMol):
    Sys += MolType.New()


Cut = 10.00
print "Cutoff is %f" % Cut

#set the system box length sizes
BoxL = 2*Cut + 3.0
Sys.BoxL = BoxL


## Adding a Localdensity and a Soft-sphere potential
FilterAA = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeA], Ordered = True)
P1 = sim.potential.LocalDensity2(Sys, Cut = Cut, InnerCut = 0.8*Cut,
                                 Knots = [0., 1., -2., 3., 4., -5.],
                                 RhoMin = 0.01, RhoMax = 2.0, 
                                 Label = "LD", Filter = FilterAA)

P2 = sim.potential.SoftSphere(Sys, Sigma = 0.01, Epsilon = -100.0, Label = "SS",
                              Cut = Cut, Shift = True, Filter = sim.atomselect.Pairs)

Sys.ForceField.extend([P1])


Sys.TempSet = TempSet
Int = Sys.Int
Int.Method = Sys.Int.Methods.VVIntegrate
Int.Method.Thermostat = Int.Method.ThermostatLangevin


Sys.Measures.PEnergy.SetupHist(0, 1000, 10)


Sys.Load()

#set initial positions and velocities
sim.system.init.positions.CubicLattice(Sys, Random = 0.1)
sim.system.init.velocities.Canonical(Sys, Temp = TempSet)


##### Log LD potential for a set of positions 'r' using sim 

print 'Writing sim computed local density potential to file...'
r_min = 2.0
r_max = BoxL + 2.
delta = (r_max - r_min)/500
f = open('LD_potential_sim.txt', 'w')
f.write("%11s %11s" % ("# rij", "ene_locdens"))
for rij in np.arange(r_min, r_max, delta):
    f.write("%11.5f %11.5f %11.5f" % (rij, P1.Val(P1.RhoVal(rij)), P1.DVal(P1.RhoVal(rij))))
    f.write('\n')
f.close()

 
#test to export to lammps
try:
	print "Creating and running Lammps script..."
	if os.path.isfile('parselog.txt'): os.remove('parselog.txt')
	if os.path.isfile('forcelog.txt'): os.remove('forcelog.txt')	
	sim.export.lammps.LammpsExec = LammpsExec
	Traj, TrajFile = sim.export.lammps.MakeLammpsTraj(Sys, Prefix = "LD", NStepsEquil = 10000, 
                                                 NStepsProd = 40000, WriteFreq = 10,
                                                 DelTempFiles = False)
except sim.export.lammps.LammpsError:	
	print sim.export.lammps.LammpsError
	compare_values = False
	

##### LOG A COMPARISON OF THE LAMMPS DUMP FILE WITH ENERGIES AND FORCES COMPUTED BY SIM FROM THE LAMMPS TRAJ 
if compare_values:
	print 'Comparing Lammps energies with sim energies...'
	Sys.Flags.Alg.CalcTerms = True  
	f = open('compare.txt', 'w ')
	f.write("%11s "*6  % ("rij", "rhoi&rhoj", "ene_softsph", "ene_locdens", "ene_total", "lammps_ene"))
	f.write('\n')                                           
	for Pos in Traj:
	    # Minimum imaging co-ordinates
	    rij = sim.geom.Length(sim.geom.Minimage(Pos[0] - Pos[1], Sys.BoxL)) 
	    if len(Sys.ForceField)==2:
			SSEne = P2.Val(rij)
	    else:
			SSEne = 0.0
	    rhoi = P1.RhoVal(rij)
	    LDEne = P1.Val(rhoi)
	    #Sys.Arrays.Pos = Pos
	    #Sys.ForceField.Eval()
	    f.write ("%11.4e "*6 % (rij, rhoi, SSEne, LDEne*2, SSEne+2*LDEne, Traj.FrameData["PEnergy"]))
	    f.write('\n')
	#TrajEne, SysEne, TrajLogWeight = sim.srel.base.GetEneTraj(Traj, Sys,
        #                                                 ErrorNoEne = True, ErrorDiffEne = True)



######## BUILD ANOTHER LAMMPS SCRIPT JUST TO RUN THE PAIR_WRITE COMMAND THIS TIME AND LOG THE DATA 
if doPairWrite:
	if not os.path.isfile('LDlammps.in'):
		raise ValueError('Please run the rest of this script to generate the Lammps input script, first')
	
	print 'Re-running with pair_write...'
	if os.path.isfile('forcefield.txt'): os.remove('forcefield.txt')
	s = open('LDlammps.in', 'r').readlines()
	start = 0
	stop = [s.index(line) for line in s if line.startswith('special_bonds')][0] - 1
	cmdstring = 'pair_write 1 1 80 r %g %g LD_potential_lammps.txt LD' % (r_min, r_max)
	cmdstring += '\npair_modify shift yes'
	
	s_new = s[start:stop]
	s_new.append(cmdstring)
	f = open('LDlammps_new.in', 'w')
	[f.write(this) for this in s_new]
	f.close()
	os.system('%s -in LDlammps_new.in -log LDlammps.log' % LammpsExec)





