#####################################################################################################
##          libraries
#####################################################################################################

from pyevtk.hl import gridToVTK
from scipy.linalg import norm

from scipy import signal

import numpy as np

from numpy import linalg as la

import csv
import os

import time as clock
start_time = clock.time()

#####################################################################################################
##          Script
#####################################################################################################

# size of the snapshot
MM = 381600
# number of snapshots
N = 5

#Definition of the grid for writing a vtk file with this solution
nx, ny, nz = 106, 60, 60
lx, ly, lz = 0.07, 0.04, 0.04
dx, dy, dz = lx/(nx-1), ly/(ny-1), lz/(nz-1)

ncells = nx * ny * nz

# Coordinates
X = np.arange(0, lx+dx, dx, dtype='float64')
Y = np.arange(0, ly+dx, dy, dtype='float64')
Z = np.arange(0, lz+dx, dz, dtype='float64')

x = np.zeros((nx, ny, nz))
y = np.zeros((nx, ny, nz))
z = np.zeros((nx, ny, nz))

Ux_Mode = np.zeros((nx, ny, nz))
Uy_Mode = np.zeros((nx, ny, nz))
Uz_Mode = np.zeros((nx, ny, nz))

#Store the POD modes for each scalar
POD_mode_Ux = np.zeros([MM,N],'float')
POD_mode_Uy = np.zeros([MM,N],'float')
POD_mode_Uz = np.zeros([MM,N],'float')

print '\nReading POD modes from txt files...\n'

for i in range(0,N,1):
  POD_mode_Ux[:,i] = np.loadtxt('output/sortedModes/U/mode_Ux_%s.dat' % i, skiprows=0)
  POD_mode_Uy[:,i] = np.loadtxt('output/sortedModes/U/mode_Uy_%s.dat' % i, skiprows=0)
  POD_mode_Uz[:,i] = np.loadtxt('output/sortedModes/U/mode_Uz_%s.dat' % i, skiprows=0)
  print 'Finished reading POD mode %.0f of %.0f' % (i+1, N)

print '\nWriting the modes in VTK files...\n'

for idx in range(N):
   for k in range(nz):
      for i in range(nx):
         for j in range(ny):
            x[i, j, k] = X[i]
            y[i, j, k] = Y[j]
            z[i, j, k] = Z[k]

            Ux_Mode[i, j, k] = POD_mode_Ux[i*ny+j+k*ny*nx, idx]
            Uy_Mode[i, j, k] = POD_mode_Uy[i*ny+j+k*ny*nx, idx]
            Uz_Mode[i, j, k] = POD_mode_Uz[i*ny+j+k*ny*nx, idx]
   U_Mode = (Ux_Mode, Uy_Mode, Uz_Mode)
   gridToVTK("output/VTK/mode_U_%s" %(idx), x, y, z, pointData={"U": U_Mode})
   print 'Finished writing POD mode %.0f of %.0f' % (idx+1, N)

print '\nDone\n'