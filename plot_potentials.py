#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt


compare = True      # flag to plot the comparison plot between Lammps dump file and values computed by sim using the Lammps traj
doPairWrite = True  # flag to plot comparison between sim and Lammps pair_write
doParseLog = True   # flag to plot data parsed from input LD file supplied by user. Used to verify if the parsing is correct
plot_force = True   # flag to plot the force (derivative of the LD potential)

def plot1():
	"""
	Plot the soft sphere and LD energies computed in sim, using the Lammps trajectory.
	Also compares the total energies computed by sim from the Lammps traj and the Lammps dump file
	"""
	if compare and os.path.isfile('compare.txt'):
		data = np.loadtxt('compare.txt', skiprows = 1)
		data = data[np.argsort(data[:,0])]	
		r = data[:,0]
		e_ss = data[:,2]
		e_ld = data[:,3]
		e_tot = data[:,4]
		e_lmp = data[:,5]
	
		plt.figure(figsize = (10,10))
		sp = plt.subplot(221)
		sp.set_xlabel('r (LJ units)', fontsize = 20)
		sp.set_ylabel('Energy (kT)', fontsize = 20)
		sp.plot(r, e_ss, linewidth = 3)
		sp.set_title('Ene_SS with sim using Lammps traj')
		
		sp = plt.subplot(222)
		sp.set_xlabel('r (LJ units)', fontsize = 20)
		sp.set_ylabel('Energy (kT)', fontsize = 20)
		sp.plot(r,e_ld, linewidth = 3)
		sp.set_title('Ene_LD with sim using Lammps traj')
		
		sp = plt.subplot(223)
		sp.set_xlabel('r (LJ units)', fontsize = 20)
		sp.set_ylabel('Energy (kT)', fontsize = 20)
		sp.plot(r,e_tot, linewidth = 1, linestyle = 'solid', color = 'red', label = 'Total PE sim')
		sp.plot(r,e_lmp, linewidth = 1, marker = 'o', markersize = 4, markerfacecolor= 'None', linestyle = 'None', color = 'blue', label = 'Total PE Lammps')
		sp.legend(loc = 4)

		#plt.tight_layout

def plot2():
	"""
	Compare the LD potential and the force as drawn with sim
	and with pair_write function of Lammps
	"""
	if doPairWrite and os.path.isfile('LD_potential_sim.txt') and os.path.isfile('LD_potential_lammps.txt'):
		data1 = np.loadtxt('LD_potential_sim.txt'); 
		data2 = np.loadtxt('LD_potential_lammps.txt', skiprows = 5)
		plt.figure(figsize = (10,10))
		
		sp = plt.subplot(211)
		sp.set_xlabel('r(LJ units)', fontsize = 20); sp.set_ylabel(r'$U_{LD} (r) \mathrm{kT}$', fontsize = 20)
		sp.plot(data1[:,0], data1[:,1], linestyle = 'solid', linewidth = 3, color = 'green', label = 'Energy - Sim')
		if plot_force:	sp.plot(data1[:,0], data1[:,2], linestyle = 'solid', linewidth = 3, color = 'red', label = 'Force - Sim')
		sp.legend(loc = 4, prop = {'size' : 12})
		
		sp = plt.subplot(212)
		sp.set_xlabel('r(LJ units)', fontsize = 20); sp.set_ylabel(r'$U_{LD} (r) \mathrm{kT}$', fontsize = 20)
		sp.plot(data2[:,1], data2[:,2], linewidth = 3, color = 'green', label = 'Energy - Lammps (pair_write)')
		if plot_force:	sp.plot(data2[:,1], data2[:,3], linewidth = 3, color = 'red', label = 'Force-Lammps(pair_write)')
		sp.legend(loc = 4, prop = {'size': 12})
	
		#plt.tight_layout





def plot3():
	"""
	-----------------Test of the Lammps interpolation routine------------------------
	"""
	if doParseLog and os.path.isfile('parselog.txt'):
		
		lines = open('parselog.txt', 'r').readlines()
	
		
		### EXTRACTING SPLINE COEFFICIENTS CALCULATED INSIDE THE LAMMPS ROUTINE - interpolate()
		start = [lines.index(line) for line in lines if line.startswith('>>> LD POTENTIAL 1')][0] + 1
		stop = start + 500
		rho = np.zeros([500], np.float64)
		frho = np.zeros([500], np.float64)
		for i, line in enumerate(lines[start:stop]):
			rho[i] = float(line.split()[0])
			frho[i] = float(line.split()[1])
		
		start = [lines.index(line) for line in lines if line.startswith('FRHO SPLINE COEFFICIENTS')][0] + 2
		stop = start + 500
		frho_spline = np.zeros([500,7], np.float64)
		for i, line in enumerate(lines[start:stop]):
			for j in range(7):
				frho_spline[i][j] = float(line.split()[j])


		
		### REPLICATING THE LAMMPS ROUTINE interpolate()
		nrho = 500
		delta = rho[1]-rho[0]
		rho_min = 0.01
		rho_max = 2.0		

		spline = np.zeros([nrho, 7], np.float64)
		for m in range(nrho): 
			spline[m,6] = frho[m]
		
		spline[0, 5] = spline[1, 6] - spline[0, 6]
	 	spline[1, 5] = 0.5 * (spline[2,6]-spline[0,6])
  		spline[nrho-2,5] = 0.5 * (spline[nrho-1, 6]-spline[nrho-3,6])
  		spline[nrho-1,5] = spline[nrho-1,6] - spline[nrho-2,6]
		
		for m in range(2,nrho-3):
    			spline[m,5] = ((spline[m-2,6]-spline[m+2,6]) + 8.0*(spline[m+1,6]-spline[m-1,6])) / 12.0
		
		for m in range(nrho-2):
    			spline[m,4] = 3.0*(spline[m+1,6]-spline[m,6]) -2.0*spline[m,5] - spline[m+1,5];
    			spline[m,3] = spline[m,5] + spline[m+1,5] -2.0*(spline[m+1,6]-spline[m,6]);
 
		spline[nrho-1,4] = 0.0
  		spline[nrho-1,3] = 0.0

  		for m in range(nrho):
    			spline[m,2] = spline[m,5]/delta;
	    		spline[m,1] = 2.0*spline[m,4]/delta;
    			spline[m,0] = 3.0*spline[m,3]/delta;
	
		


		### LD POTENTIAL CALCULATION USING A LAMMPS-LIKE CUBIC SPLINE INTERPOLATION
		uLD1 = np.zeros([nrho], np.float64); uLD2 = np.zeros([nrho], np.float64)
		dFdrho1 = np.zeros([nrho], np.float64); dFdrho2 = np.zeros([nrho], np.float64)
		p = 0.; m = 0;
		for i in range(nrho):
			p = (rho[i] - rho_min)/delta + 1
			m = max(1,min(m,nrho-1))
			p -= m
			p = min(p, 1.0)
			coeff1 = frho_spline[i]
			coeff2 = spline[i]
			uLD1[i] = ((coeff1[3]*p + coeff1[4])*p + coeff1[5])*p + coeff1[6]
			uLD2[i] = ((coeff2[3]*p + coeff2[4])*p + coeff2[5])*p + coeff2[6]
			dFdrho1[i] = (coeff1[0]*p + coeff1[1])*p + coeff1[2];
			dFdrho2[i] = (coeff2[0]*p + coeff2[1])*p + coeff2[2];
			
		
		#### REPLICATING THE LAMMPS PAIR WRITE FUNCTION 
		r = np.linspace(2.2, 4.0, 400)
		LD = np.zeros([400], np.float64)
		ene_LD = np.zeros([400], np.float64)
		force_LD = np.zeros([400], np.float64)
	
		
		def get_phi(r, rsq, R1, R2, option):
			lowercutsq =  R1*R1
 			uppercutsq =  R2*R2
 			uppercutfourth = uppercutsq*uppercutsq
 			uppercutsixth   = uppercutfourth*uppercutsq
 			cut_ratio = lowercutsq/uppercutsq
 			denom = (1-cut_ratio)**3
 			c0 = (1 - 3*cut_ratio) / denom 
 			c2 =  (1/uppercutsq) * (6*cut_ratio) / denom
 			c4 = -(1/uppercutfourth) * (3 + 3*cut_ratio) / denom
 			c6 = (1/uppercutsixth) * 2.00 / denom

  			if (this_r < R1): 
				phi = 1.0
				dphidr = 0.
  			if (this_r > R2): 
				phi = 0.0
				dphidr = 0
			if (this_r >= R1 and this_r <= R2): 
				phi = c0 + rsq * (c2 + rsq * (c4 + c6*rsq));
  				dphidr = (1./this_r) * (rsq * (2*c2 + rsq * (4*c4 + 6*c6*rsq)));


			if not option: return phi
			else:	       return dphidr	



		R2 = 3.0; R1 = 0.8 * R2
		phi = 0.; dphidr = 0.
		dFdrho = 0.0
		for i, this_r in enumerate(r):
			rsq = this_r*this_r
			phi = get_phi(this_r, rsq, R1, R2, 0)
			dphidr = get_phi(this_r, rsq, R1, R2, 1)
			
			LD[i] += phi
			
			
			if (LD[i] <= rho_min):
				p = 0.0; m = 1
				coeff = frho_spline[i]
				dFdrho = coeff[2]
				ene_LD[i] = coeff[6] + dFdrho*(LD[i] - rho_min) 
		
			elif (LD[i] >= rho_max):
				p = 1.0; m = nrho - 1
				coeff = frho_spline[i];
				dFdrho = coeff[0] + coeff[1]  + coeff[2]
				ene_LD[i] = (coeff[3] + coeff[4] + coeff[5] + coeff[6]) + dFdrho*(LD[i] - rho_max) ;
		
			else:
				p = (LD[i] - rho_min) / delta + 1.0
				m = int(p)
				m = max(1,min(m, nrho-1))
				p -= m;
				p = min(p, 1.0);
				coeff = frho_spline[i];
				dFdrho = (coeff[0]*p + coeff[1])*p + coeff[2]
      				ene_LD[i] =  ((coeff[3]*p + coeff[4])*p + coeff[5])*p + coeff[6]
		
			force_LD[i] =  -dFdrho * dphidr / this_r;

			 
 
	

		plt.figure(figsize = (10,10))
		sp = plt.subplot(221)
		sp.set_title('Verify Spline Interpolation')
		sp.plot(rho, frho, linestyle = 'solid', linewidth = 3, color = 'red', label = 'F(rho_LD) supplied')
		sp.plot(rho, uLD1, linestyle = 'solid', linewidth = 3, color = 'blue', label = 'F(rho_LD) from Lammps')
		sp.plot(rho, uLD2, linestyle = 'solid', linewidth = 3, color = 'green', label = 'F(rho_LD) from this code')
		sp.set_xlabel('rho_LD'); sp.set_ylabel('F(rho_LD)')		
		sp.legend()

		
		sp = plt.subplot(222)
		sp.set_title('Verify dF/drho')
		sp.plot(rho, dFdrho1, linestyle = 'solid', linewidth = 3, color = 'blue', label = 'dF/drho from Lammps')
		sp.plot(rho, dFdrho2, linestyle = 'solid', linewidth = 3, color = 'green', label = 'dF/drho from this code')	
		sp.set_xlabel('rho_LD'); sp.set_ylabel('dF/drho')			
		sp.legend()
		
		sp = plt.subplot(223)
		sp.set_title('Verify energy and force')
		sp.plot(r, ene_LD, linestyle = 'solid', linewidth = 3, color = 'green', label = 'ene_LD from this code')
		if plot_force: sp.plot(r, force_LD, linestyle = 'solid', linewidth = 3, color = 'red', label = 'force_LD from this code')
		sp.set_xlabel('r (LJ unit)'); sp.set_ylabel('ene_LD (lJ unit)')			
		sp.legend()

		sp = plt.subplot(224)
		sp.set_title('Verify Local Density')
		sp.plot(r, LD, linestyle = 'solid', linewidth = 3, color = 'green', label = 'Local Density from this code')
		sp.set_xlabel('r (LJ unit)'); sp.set_ylabel('Local-Density')			
		sp.legend()
		
		#plt.tight_layout
	



## MAIN ## 

plot1()
plot2()
plot3()

plt.show()





