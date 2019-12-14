"""
Created on Tue Dec 10 16:44:36 2019

@author: Vedang
"""

#==============================================================================
# LIBRARIES
#==============================================================================

# imports libraries
from math import pi
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import time

#==============================================================================
# FUNCTIONS
#==============================================================================

# computes Euclidean distance considering periodic boundaries
def periodic_distance(pos_1, pos_2):
    dx = pos_1[0] - pos_2[0]
    if (abs(dx) > box_length*0.5):
       dx = box_length - dx
    dy = pos_1[1] - pos_2[1]
    if (abs(dy) > box_length*0.5):
       dy = box_length - dy
    dz = pos_1[2] - pos_2[2]
    if (abs(dz) > box_length*0.5):
       dz = box_length - dz
    distance = np.sqrt(dx**2 + dy**2 + dz**2)
    return distance

#==============================================================================
# INITIALIZATION
#==============================================================================

# starts stopwatch for execution time
start_time = time.time()

# enter known values here
n_particles = 500  # number of particles
rho = 0.95  # density

# sets number of bins
n_bins = 100

# calculates box length
box_length = pow(n_particles/rho, 1/3)

# calculates bin size
bin_size = box_length/(2*n_bins)

# initializes sample count 
sample_count = 0

# initializes g(r) list for all bins
g = [0] * (n_bins)

#==============================================================================
# SAMPLING
#==============================================================================

# reads all coordinate files in folder
rho_file_path = '/Users/Vedang/Desktop/Nanosystems/supporting_materials/rho_0.95'
list_dir = os.listdir(rho_file_path)
list_dir = [f.lower() for f in list_dir]   # converts to lower case
for coords_name in sorted(list_dir):
    if coords_name.startswith('data.coords.all.'):
        raw_data = pd.read_csv(os.path.join(rho_file_path, coords_name), header=None, skiprows = 9, usecols=[3,4,5], sep='\s+')
        coords = np.asarray(raw_data.values.tolist())
        sample_count = sample_count + 1
        for i in range(0, n_particles-1):
            for j in range(i+1, n_particles):
                absolute_distance = periodic_distance(coords[i], coords[j])
                if absolute_distance < (box_length/2):
                    int_absolute_distance = int(absolute_distance/bin_size)  # rounding down to integer value
                    g[int_absolute_distance] = g[int_absolute_distance] + 2

#==============================================================================
# RDF CALCULATION
#==============================================================================

# initializes centering value to plot bins
r = [0] * (n_bins)

# calculates g(r)
for i in range(0, n_bins):
    r[i] = bin_size*(i+1+0.5)  # centering value
    bin_volume = (((i+1)**3) - (i**3)) * (bin_size**3)
    n_ideal = (4/3) * pi * bin_volume * rho  # number of ideal gas particles in bin volume
    g[i] = g[i]/(sample_count * n_particles * n_ideal)  
            
# plots g(r)
plt.figure()
plt.title('My Code (rho = 0.7)')
plt.plot(r, g)
plt.xlabel('r')
plt.ylabel('g(r)')
plt.legend()
plt.grid(True, alpha = 0.5)
plt.show()

# prints execution time
print("\n--- Execution Time: %s seconds ---" % (time.time() - start_time))

# compares output
comp_data = pd.read_csv('/Users/Vedang/Desktop/Nanosystems/supporting_materials/rho_0.95/gr.out', header=None, sep='\s+')
comp_r = comp_data.iloc[:, 0]
comp_gr = comp_data.iloc[:, 1]

# plots g(r)
plt.figure()
plt.title('Class Code (rho = 0.7)')
plt.plot(comp_r, comp_gr)
plt.xlabel('r')
plt.ylabel('g(r)')
plt.legend()
plt.grid(True, alpha = 0.5)
plt.show()
