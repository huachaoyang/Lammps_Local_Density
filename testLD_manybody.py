#usr/bin/env python


import sys, os
sys.path.append('./')

import numpy as np
import sim
import matplotlib.pyplot as plt

LammpsExec = "/home/cask0/home/tsanyal/tanmoy_lammps/lammps-15May15/src/lmp_ZIN"

runLammps = True
compare_values = True   # flag to compare energies by re-evaluating them with sim
plot = True

Traj = ''


if runLammps:
	#parameter for number of molecules
	NMol = 30

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

	Cut = 2.5
	print "Cutoff is %f" % Cut

	#set the system box length sizes
	BoxL = 2*Cut + 3.0
	Sys.BoxL = BoxL

	## Adding a Localdensity and a Soft-sphere potential
	FilterAA = sim.atomselect.PolyFilter([AtomTypeA, AtomTypeA], Ordered = True)
	P1 = sim.potential.LocalDensity2(Sys, Cut = Cut, InnerCut = 0.8*Cut,
                                 Knots = [0., -1., -2., -3., -4., -5.],
                                 RhoMin = 0., RhoMax = 2.0, 
                                 Label = "LD", Filter = FilterAA)

	P2 = sim.potential.SoftSphere(Sys, Sigma = 1.0, Epsilon = 1.0, Label = "SS",
                              Cut = Cut, Shift = True, Filter = sim.atomselect.Pairs)

	Sys.ForceField.extend([P1,P2])

	Sys.TempSet = TempSet
	Int = Sys.Int
	Int.Method = Sys.Int.Methods.VVIntegrate
	Int.Method.Thermostat = Int.Method.ThermostatLangevin

	Sys.Measures.PEnergy.SetupHist(0, 1000, 10)
	Sys.Load()

	#set initial positions and velocities
	sim.system.init.positions.CubicLattice(Sys, Random = 0.1)
	sim.system.init.velocities.Canonical(Sys, Temp = TempSet)

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
	



if compare_values:
	if not Traj:	
		Trajfile = raw_input('Enter trajfile name:' )
		Traj = sim.traj.lammps.Lammps(Trajfile)
	
	print 'Comparing Lammps energies with sim energies...'
	Sys.Flags.Alg.CalcTerms = True  
	f = open('compare.txt', 'w ')
	f.write("%11s "*4  % ("ene_softsph", "ene_locdens", "ene_total", "lammps_ene"))
	f.write('\n')

	NFrames = len(Traj)                           
	
	pb = sim.utility.ProgressBar(Text = 'Processing frame by frame...', Steps = NFrames)
	for frame in range(NFrames):
		SSEne = 0.
		LDEne = 0.
		Pos = Traj[frame]

	    # calculate soft-sphere potential
		if len(Sys.ForceField) == 2:
			for i in range(len(Pos-1)):
				for j in range(i+1,len(Pos)):
					rij = sim.geom.Length(sim.geom.Minimage(Pos[j] - Pos[i], Sys.BoxL)) 
					SSEne += P2.Val(rij)

	    # calculate the localdensity potential
		for i in range(len(Pos)):
			rhoi = 0.
			for j in range(len(Pos)):
				if not i==j:	
					rij = sim.geom.Length(sim.geom.Minimage(Pos[j] - Pos[i], Sys.BoxL)) 
					rhoi += P1.RhoVal(rij)
			LDEne += P1.Val(rhoi)
		Sys.Arrays.Pos = Pos
		Sys.ForceField.Eval()
		
		f.write ("%11.6e "*4 % (SSEne, LDEne, SSEne+LDEne, Traj.FrameData["PEnergy"]))
		f.write('\n')
		pb.Update(frame)
	#TrajEne, SysEne, TrajLogWeight = sim.srel.base.GetEneTraj(Traj, Sys,
    #                                                          ErrorNoEne = True, ErrorDiffEne = True)



if plot:
		if not os.path.isfile('compare.txt'):
				raise IOError('Data file missing')

		WriteFreq = 10
		data = np.loadtxt('compare.txt', skiprows = 1)
		PE_sim = data[:,2]
		PE_lammps = data[:,3]
		

		plt.figure()
		ax = plt.subplot('121')
		ax.plot(PE_sim, linestyle = 'None', marker = 'o', markersize = 2, markerfacecolor = 'None', color = 'blue', label = 'PE_sim')
		ax.plot(PE_lammps, linewidth = 1, color = 'red', label = 'PE_lammps')
		ax.set_ylabel('Total Potential Energy')
		ax.set_xlabel('Frame X %g' % WriteFreq)
		ax.legend()

		ax = plt.subplot('122')
		ax.plot((1 - PE_lammps/PE_sim)*100, linewidth = 2, color = 'blue', label = '(PE_sim - PE_lammps)/PE_sim')
		ax.set_xlabel('Frame X %g' % WriteFreq)
		ax.set_ylabel('Percent difference')
		plt.legend()
		plt.show()
