#####################################################################################################
##          libraries
#####################################################################################################

#from pylab import *
from pyevtk.hl import gridToVTK
from scipy.linalg import norm

from scipy import signal

import numpy as np
#import matplotlib as mpl
#import matplotlib.pyplot as plt

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
N = 2804  # 1243
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
reconstructed_Ux = np.zeros([MM,N],'float') 
reconstructed_Uy = np.zeros([MM,N],'float') 
reconstructed_Uz = np.zeros([MM,N],'float') 
mean_Ux = np.zeros(MM,'float')
mean_Uy = np.zeros(MM,'float')
mean_Uz = np.zeros(MM,'float')

chronos = np.zeros([N,N],'float')

#Tables for projection
data1 = np.zeros([3*MM],'float') #For projection
data2 = np.zeros([3*MM],'float') #For projection
projection_matrix = np.zeros([N,N],'float')

time = np.loadtxt('../../../data/longest3M/times.txt', skiprows=0)
data_time = np.zeros([N], 'float')

#Store the POD modes for each scalar
POD_mode_Ux = np.zeros([MM,N],'float')
POD_mode_Uy = np.zeros([MM,N],'float')
POD_mode_Uz = np.zeros([MM,N],'float')

print ('Create database')
for  t in range(N) :
    loop_time = clock.time()
    data_time[t] = time[t]
    Snapshot = First_Snapshot+t
    data_U = np.loadtxt('../../../data/longest3M/sortedField/U/U_%s.dat' % Snapshot, skiprows=0)
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
    print('reading data ', t, ' over ', N, ' executed in ', '%.2f' % (clock.time()-loop_time), ' s')

print(data_time_Ux)

#####Proceed to POD analysis for the Q criterion
#####################################################################################################
Loop_time = clock.time()
print ('Compute the correlation matrix for Ux')
for t1 in range(0,N,1) :
    for t2 in range(0,N,1) :
        data1 = data_time_U[:,t1]
        data2 = data_time_U[:,t2]
        projection_matrix[t1,t2]=np.dot(data1,np.transpose(data2))

projection_matrix = (1.0/N)*projection_matrix

print (projection_matrix)

[A,B]=np.linalg.eig(projection_matrix)

print('EigenValues')
#order the eigenvalues the last one is the largest here
C=np.msort(np.real(A))
X1=A.argsort(axis=0)
C[N-1] = 0.0
sumC = sum(C)
print (C/sumC)
print (A[X1])

write_A = open('./chronos/A.txt', 'w')
for t1 in range(0,N,1) :
   loop_time = clock.time()
   t3 = X1[N-(t1+1)]
   write_A.write('%4.8f ' % A[t3])
   write_A.write('\n')
   write_B = open('./chronos/B.%s.txt' % t1, 'w')
   for t2 in range(0,N,1) :
         write_B.write('%4.8f ' % (np.sqrt(A[t3]*N) * B[t2,t3]))
         write_B.write('\n')
         POD_mode_Ux[:,t1] = POD_mode_Ux[:,t1] + (1/(A[t3]*N))* (np.sqrt(A[t3]*N)) * B[t2,t3] * data_time_Ux[:,t2]
         POD_mode_Uy[:,t1] = POD_mode_Uy[:,t1] + (1/(A[t3]*N))* (np.sqrt(A[t3]*N)) * B[t2,t3] * data_time_Uy[:,t2]
         POD_mode_Uz[:,t1] = POD_mode_Uz[:,t1] + (1/(A[t3]*N))* (np.sqrt(A[t3]*N)) * B[t2,t3] * data_time_Uz[:,t2]
   print('POD mode ', t1, ' over ', N, ' executed in ', '%.2f' % (clock.time()-loop_time), ' s')

mean_Ux[:] = POD_mode_Ux[:,0]
mean_Uy[:] = POD_mode_Uy[:,0]
mean_Uz[:] = POD_mode_Uz[:,0]

for t1 in range(0,N,1) :
   loop_time = clock.time()
   t3 = X1[N-(t1+1)]
   for t2 in range(0,N,1) :
      chronos[t1,t2] = (np.sqrt(A[t3]*N) * B[t2,t3])
   print('chronos ', t1, ' over ', N, ' executed in ', '%.2f' % (clock.time()-loop_time), ' s')

for i in range(0,MM,1) :
   for t1 in range(0,N,1) : # KEEP N-1 !!! 
      reconstructed_Ux[i,0] = reconstructed_Ux[i,0] + POD_mode_Ux[i,t1]*chronos[t1,0] 
      reconstructed_Uy[i,0] = reconstructed_Uy[i,0] + POD_mode_Uy[i,t1]*chronos[t1,0] 
      reconstructed_Uz[i,0] = reconstructed_Uz[i,0] + POD_mode_Uz[i,t1]*chronos[t1,0] 

print('Write the txt files')

write_POD_time = open('./mode/time_POD.txt', 'w')
for t1 in range(0,N,1) :
   write_POD_time.write('%5e ' % data_time[t1])
   write_POD_mode_Ux = open('./mode/mode_Ux.%s.txt' % (t1), 'w')
   write_POD_mode_Uy = open('./mode/mode_Uy.%s.txt' % (t1), 'w')
   write_POD_mode_Uz = open('./mode/mode_Uz.%s.txt' % (t1), 'w')
   write_chronos = open('./chronos/chronos.%s.txt' % (t1), 'w')
   for i in range(0,MM,1) :
      write_POD_mode_Ux.write('%5e ' % POD_mode_Ux[i,(t1)])
      write_POD_mode_Ux.write('\n')
      write_POD_mode_Uy.write('%5e ' % POD_mode_Uy[i,(t1)])
      write_POD_mode_Uy.write('\n')
      write_POD_mode_Uz.write('%5e ' % POD_mode_Uz[i,(t1)])
      write_POD_mode_Uz.write('\n')
   for t2 in range(0,N,1) :
      write_chronos.write('%4.8f ' % chronos[t1,t2])
      write_chronos.write('\n')


write_reconstructed_Ux = open('./mode/reconstructed_Ux.txt', 'w')
for i in range(MM):
   write_reconstructed_Ux.write('%5e %5e' % (reconstructed_Ux[i,0], data_time_Ux[i,0]))
   write_reconstructed_Ux.write('\n')


write_mean_Ux = open('./mode/mean_Ux.txt', 'w')
for i in range(MM):
   write_mean_Ux.write('%5e ' % mean_Ux[i])
   write_mean_Ux.write('\n')

#####################################################################################################
##          Export the grid data in a vtk format readable by paraview
#####################################################################################################


print('Write the VTK files')

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

nModes = 50

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








