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
##          inputs
#####################################################################################################

# size of the snapshot
MM = 381600
# number of snapshots
N = 5
# number of modes to write
nModes = N
# number of the first snapshot
First_Snapshot=1


#####################################################################################################
##          POD analysis
#####################################################################################################

#####1 Read the data and store them in matrix data_time#####
#####################################################################################################
data_Ux = np.zeros(MM,'float')
data_Uy = np.zeros(MM,'float')
data_Uz = np.zeros(MM,'float')
data_time_U = np.zeros([3*MM,N],'float') 
data_time_Ux = np.zeros([MM,N],'float') 
data_time_Uy = np.zeros([MM,N],'float') 
data_time_Uz = np.zeros([MM,N],'float') 

chronos = np.zeros([N,N],'float')

#Tables for projection
data1 = np.zeros([3*MM],'float') #For projection
data2 = np.zeros([3*MM],'float') #For projection
projection_matrix = np.zeros([N,N],'float')

time = np.loadtxt('input/times.txt', skiprows=0)
data_time = np.zeros([N], 'float')

#Store the POD modes for each scalar
POD_mode_Ux = np.zeros([MM,N],'float')
POD_mode_Uy = np.zeros([MM,N],'float')
POD_mode_Uz = np.zeros([MM,N],'float')

print 'Creating database...\n'
for  t in range(N) :
    loop_time = clock.time()
    data_time[t] = time[t]
    Snapshot = First_Snapshot+t
    data_U = np.loadtxt('input/U/U_%s.dat' % Snapshot, skiprows=0)
    data_Ux = data_U[:, 0]
    data_Uy = data_U[:, 1]
    data_Uz = data_U[:, 2]
    for i in range(MM) :
        data_time_Ux[i,t] = data_Ux[i] 
        data_time_Uy[i,t] = data_Uy[i] 
        data_time_Uz[i,t] = data_Uz[i] 
        data_time_U[i+0*MM,t] = data_Ux[i]
        data_time_U[i+1*MM,t] = data_Uy[i]
        data_time_U[i+2*MM,t] = data_Uz[i]
    print 'Reading input data %.0f of %.0f... executed in %.2f s' % (t, N, clock.time()-loop_time)

#####Proceed to POD analysis for the Q criterion
#####################################################################################################
Loop_time = clock.time()
print '\nComputing the correlation matrix for U... \n'
for t1 in range(0,N,1) :
    for t2 in range(0,N,1) :
        data1 = data_time_U[:,t1]
        data2 = data_time_U[:,t2]
        projection_matrix[t1,t2]=np.dot(data1,np.transpose(data2))

projection_matrix = (1.0/N)*projection_matrix

print 'Calculating eigenvalues...\n'

[A,B]=np.linalg.eig(projection_matrix)

#order the eigenvalues the last one is the largest here
C=np.msort(np.real(A))
X1=A.argsort(axis=0)
C[N-1] = 0.0
sumC = sum(C)

write_A = open('chronos/A.txt', 'w')
for t1 in range(0,N,1) :
   loop_time = clock.time()
   t3 = X1[N-(t1+1)]
   write_A.write('%4.8f \n' % A[t3])
   if t1 <= nModes:
      write_chronos = open('./chronos/chronos.%s.txt' % t1, 'w')
   for t2 in range(0,N,1) :
         if t1 <= nModes:
            write_chronos.write('%4.8f \n' % (np.sqrt(A[t3]*N) * B[t2, t3]))
         POD_mode_Ux[:,t1] = POD_mode_Ux[:,t1] + (1/(A[t3]*N))* (np.sqrt(A[t3]*N)) * B[t2,t3] * data_time_Ux[:,t2]
         POD_mode_Uy[:,t1] = POD_mode_Uy[:,t1] + (1/(A[t3]*N))* (np.sqrt(A[t3]*N)) * B[t2,t3] * data_time_Uy[:,t2]
         POD_mode_Uz[:,t1] = POD_mode_Uz[:,t1] + (1/(A[t3]*N))* (np.sqrt(A[t3]*N)) * B[t2,t3] * data_time_Uz[:,t2]
   print 'Computing POD mode %.0f of %.0f... executed in %.2f s' % (t1, N, clock.time()-loop_time)

print '\nWriting the modes in txt files...\n'

write_POD_time = open('mode/time_POD.txt', 'w')

for t1 in range(0,nModes,1) :
   write_POD_time.write('%5e \n' % data_time[t1])
   write_POD_mode_Ux = open('./mode/mode_Ux.%s.txt' % (t1), 'w')
   write_POD_mode_Uy = open('./mode/mode_Uy.%s.txt' % (t1), 'w')
   write_POD_mode_Uz = open('./mode/mode_Uz.%s.txt' % (t1), 'w')
   for i in range(0,MM,1) :
      write_POD_mode_Ux.write('%5e \n' % POD_mode_Ux[i,(t1)])
      write_POD_mode_Uy.write('%5e \n' % POD_mode_Uy[i,(t1)])
      write_POD_mode_Uz.write('%5e \n' % POD_mode_Uz[i,(t1)])

#####################################################################################################
##          Export the grid data in a vtk format readable by paraview
#####################################################################################################

print 'Writing the modes in VTK files...\n'

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

for idx in range(nModes):
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
   gridToVTK("VTK/mode_U_%s" %(idx), x, y, z, pointData={"U": U_Mode})