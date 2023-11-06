### PREAMBLE ###
# Philip Daniel Blocher
# Quantum Optics Group
# Department of Physics and Astronomy
# Aarhus University, Denmark
# Email: philipblocher@gmail.com or pblocher@phys.au.dk
# Last updated: May 2019

### PROGRAM DESCRIPTION ###
# This program plots the data obtained from the program 'LoschmidtEcho.py' in order to visualize
# the Loschmidt echo in a chaotic spin 7/2 system.

### IMPORTED LIBRARIES AND GLOBAL VARIABLES ###
import math, cmath
import random
import numpy as np
import os, time
from scipy import linalg
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
pi = math.pi


### FIGURE PARAMETERS ###
# If saveFigure = True, Loschmidt echo figures are saved in the folder address '../figures/'.
saveFigure = True
saveLocation = os.getcwd() + '/LE figures/'


# If saveFigure = true and plotTrajectory = True, the Loschmdt echo trajectory is saved in the folder address '../LE trajectory figures/'
saveLocation2 = os.getcwd() + '/LE trajectory figures/'

### IMPORT DATA ###
# Data is stored in 'LEdata.npy':
data = np.load('LEdata.npy',allow_pickle=True)

# Unpack the important variables:
LE = data[0]
supportargs = data[1]

# Time-related variables:
f = supportargs[0]
tStart = supportargs[2]
tEnd = supportargs[3]
tPeriod = (1/f)
NPeriods = int((tEnd-tStart)*f)

# Search-related variables:
phiLin = supportargs[4]
thetaLin = supportargs[5]
phiN = len(phiLin)
thetaN = len(thetaLin)

# Print some info to the user:
sysp = supportargs[-1]
I = sysp[0]

if supportargs[1] == 'rotating':
	print('Loschmidt echo in rotating frame for parameters (I, alpha, beta, gamma, f) = ({I}, {alpha}, {beta}, {gamma}, {f}).'.format(I=sysp[0],alpha=sysp[1],beta=sysp[2],gamma=sysp[3],f=sysp[4]))
elif supportargs[1] == 'lab':
	print('Loschmidt echo in lab frame for parameters (I, gamma_n, A, B0, B1, f) = ({I}, {gn}, {A}, {B0}, {B1}, {f}).'.format(I=sysp[0],gn=sysp[1],A=sysp[2],B0=sysp[3],B1=sysp[4],f=sysp[5]))
else:
	raise Exception('Data did not originate from either lab or rotating frame!!')


### PLOTTING LOOP ###
# We here plot the Loschmidt echo for increasing time evolution periods t.
print('- - - - - - - - - -')
print('Now plotting Loschmidt echo values for a range of time durations t:')
print('- - - - - - - - - -')

# We first determine the min and max values for all Loschmidt echos:
z_min, z_max = np.min(np.real(LE[:,:,0])), np.max(np.real(LE[:,:,0]))

for ti in range(NPeriods+1):
	z_min = min(z_min, np.min(np.real(LE[:,:,ti])))
	z_max = max(z_max, np.max(np.real(LE[:,:,ti])))

# Now plot a pcolor figure with the LE data for each time step:
X = np.array([x/pi for x in phiLin])
Y = np.array([y/pi for y in thetaLin])
X, Y = np.meshgrid(X,Y)
statusPrint = 0

# Keep track of elapsed time:
tic = time.time()

for ti in range(NPeriods+1):
	t = tStart + ti*(1/f)

	fig = plt.figure(figsize=(24, 12))
	ax = plt.axes()

	c = plt.pcolormesh(X.T, Y.T, (LE[:,:,ti]), cmap='magma', vmin=z_min, vmax=z_max)
	ax.set_ylabel(r'$\theta/\pi$')
	ax.set_xlabel(r'$\phi/\pi$')
	ax.set_title('I = {I}, t = {t}'.format(I=I,t=round(t*100)/100))
	ax.axis([X.min(), X.max(), Y.min(), Y.max()])
	plt.colorbar(c, ax=ax)

	if saveFigure == True:
		if statusPrint == 0:
			if not os.path.exists(saveLocation):
				os.mkdir(saveLocation)
				print('\tDirectory ' + saveLocation + ' created.')
			else:
				print('\tDirectory ' + saveLocation + ' already exists.')

			print('\tPrinting figure to ' + saveLocation + ' !\n') 

			statusPrint = 1

		if ti%5 == 0 and ti>0:
			toc = time.time() - tic
			print('\tPrinting time step {ti} out of {NPeriods} after elapsed time {toc:.1f} seconds'.format(ti=ti,NPeriods=NPeriods+1,toc=toc))
			statusPrint2 = 1

		fig.savefig(saveLocation + 'TimeStep{n}.png'.format(n=ti+1))
	
	plt.close()

# We also want to do an average over all Loschmidt echos to see which initial states |phi,theta> are in particular disturbed:
LE_avr = np.real(LE[:,:,0])

for ti in range(1,NPeriods+1):
	LE_avr = LE_avr + np.real(LE[:,:,ti])

LE_avr = LE_avr / (NPeriods+1)

# We then calculate the extrema of LE_avr:
z_min, z_max = np.min(LE_avr), np.max(LE_avr)

for ti in range(NPeriods+1):
	z_min = min(z_min, np.min(LE_avr))
	z_max = max(z_max, np.max(LE_avr))

fig = plt.figure(figsize=(24, 12))
ax = plt.axes()

c = plt.pcolormesh(X.T, Y.T, LE_avr, cmap='magma', vmin=z_min, vmax=z_max)
ax.set_ylabel(r'$\theta/\pi$')
ax.set_xlabel(r'$\phi/\pi$')
ax.set_title('I = {I}, average Loschmidt echo'.format(I=I))
ax.axis([X.min(), X.max(), Y.min(), Y.max()])
plt.colorbar(c, ax=ax)

if saveFigure == True:
	fig.savefig(saveLocation + 'TimeStepAverage.png')


tic = time.time()
print('- - - - - - - - - -')
print('Now plotting Loschmidt echo trajectory:')
print('- - - - - - - - - -')

# We now want to plot a single trajectory as defined in 'LoschmidtEcho.py':
# Data is stored in 'LEtrajectory.npy':
data = np.load('LEtrajectory.npy',allow_pickle=True)

# Unpack the data (we already unpacked the system parameters previously):
LEtrajectory = data[0]
supportargs = data[1]

phi0 = supportargs[0]
theta0 = supportargs[1]

# We first determine the min and max values for all Loschmidt echos:
z_min, z_max = np.min(np.real(LEtrajectory[:,:,0])), np.max(np.real(LEtrajectory[:,:,0]))

for ti in range(NPeriods+1):
	z_min = min(z_min, np.min(np.real(LEtrajectory[:,:,ti])))
	z_max = max(z_max, np.max(np.real(LEtrajectory[:,:,ti])))

statusPrint = 0

# Keep track of elapsed time:
tic = time.time()

for ti in range(2*NPeriods+2):
	t = tStart + ti*(1/f)

	fig = plt.figure(figsize=(24, 12))
	ax = plt.axes()

	c = plt.pcolormesh(X.T, Y.T, np.real(LEtrajectory[:,:,ti]), cmap='magma', vmin=z_min, vmax=z_max)
	ax.set_ylabel(r'$\theta/\pi$')
	ax.set_xlabel(r'$\phi/\pi$')
	ax.set_title('I = {I}, t = {t}'.format(I=I,t=round(t*100)/100))
	ax.axis([X.min(), X.max(), Y.min(), Y.max()])
	plt.colorbar(c, ax=ax)

	if saveFigure == True:
		if statusPrint == 0:
			if not os.path.exists(saveLocation2):
				os.mkdir(saveLocation2)
				print('\tDirectory ' + saveLocation2 + ' created.')
			else:
				print('\tDirectory ' + saveLocation2 + ' already exists.')

			print('\tPrinting figure to ' + saveLocation2 + ' !\n') 

			statusPrint = 1

		if ti%5 == 0 and ti>0:
			toc = time.time() - tic
			print('\tPrinting time step {ti} out of {NPeriods} after elapsed time {toc:.1f} seconds'.format(ti=ti,NPeriods=2*NPeriods+2,toc=toc))
			statusPrint2 = 1

		fig.savefig(saveLocation2 + 'TimeStep{n}.png'.format(n=ti+1))
	
	plt.close()