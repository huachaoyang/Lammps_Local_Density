#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(facecolor = 'w', edgecolor = 'w', figsize = (5,5))
files = ['compare_2.txt', 'compare_30.txt']
titles = ['N = 2', 'N = 30']

for i, this_file in enumerate(files):
	ax = plt.subplot(1,2,i+1)
	data = np.loadtxt(this_file, skiprows = 1)
	e_lmp = data[:,-1]
	e_sim = data[:,-2]
	ax.plot(e_lmp, color = 'red', linestyle = 'None', marker = 'o', markeredgecolor = 'red', markersize = 6, label = 'Lammps')
	ax.plot(e_sim, color = 'blue', linestyle = 'solid', linewidth = 2, label = 'group-code')
	ax.set_xlabel('MD Steps', fontsize = 'medium')
	if i == 0:
		ax.legend(loc = 'best', prop = {'size': 12})
		ax.set_ylabel('Non dimensional energy', fontsize = 'medium')
	ax.set_title(titles[i])
	
plt.subplots_adjust(wspace = 0.2, hspace = 0.3, bottom = 0.12, left = 0.2)
plt.show()
