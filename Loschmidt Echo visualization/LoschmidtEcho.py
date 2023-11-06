### PREAMBLE ###
# Philip Daniel Blocher
# Quantum Optics Group
# Department of Physics and Astronomy
# Aarhus University, Denmark
# Email: philipblocher@gmail.com or pblocher@phys.au.dk
# Last updated: May 2019

### PROGRAM DESCRIPTION ###
# This program simulates a Loschmidt echo carried out for the chaotic spin-7/2
# system used by the Morello group. Different perturbations W can be chosen, however
# these perturbations should remain small and not reflect symmetries of the Hamiltonian.
# For now this simulation is carried out in a closed quantum system.

### IMPORTED LIBRARIES AND GLOBAL VARIABLES ###
import math, cmath
import random, time
import numpy as np
from scipy import linalg
from CreateSpinOperators import *
from CreateFloquetOperator import *
pi = math.pi

### HAMILTONIAN FUNCTION ###
def ChaosHamiltonian(t, args):

	# Load parameters from variable args:
	f = args[0]
	gamma_n = args[1]
	A = args[2]
	Q = args[3]
	B0 = args[4]
	B1 = args[5]
	Ix = args[6]
	Iy = args[7]
	Iz = args[8]

	# Create Hamiltonian
	H = (gamma_n*B0-A)*Iz + Q*Ix@Ix + gamma_n*B1*math.cos(2*pi*f*t)*Iy
	return H


def ChaosHamiltonianRWA(t, args):

	# Load parameters from variable args:
	f = args[0]
	alpha = args[1]
	beta = args[2]
	gamma = args[3]
	Ix = args[4]
	Iy = args[5]
	Iz = args[6]
	I = args[7]

	# Create the RWA Hamiltonian:
	H = alpha*Iz + beta*(Ix@Ix)/I + gamma*math.cos(2*pi*f*t)*Iy
	return H


# # # MAIN CODE STARTS HERE # # #

### GENERAL SYSTEM PARAMETERS ###
# Nuclear spin:
I = 7/2
N = int(2*I+1)


### SELECT WHETHER TO USE LAB FRAME OR ROTATING FRAME ###
# Select either 'lab' or 'rotating' frame:
statusFrame = 'lab'


### LAB FRAME SPECIFIC PARAMETERS ###
# Nuclear gyromagnetic ratio (in MHz/T):
gamma_n = 5.55

# Hyperfine interaction strength in lab frame (in MHz):
A = 0*101.52

# Quadrupole term (in MHz):
Q = 0.8

# Magnetic fields in lab frame (in T):
B0 = 0.5
B1 = 0.01

# AC B-field frequency in MHz:
f_lab = 5


### ROTATING FRAME SPECIFIC PARAMETERS ###
# Linear term alpha = 0.5*gamma_n*B_1I can be tuned from 0 up to 30 kHz:
alpha = 1

# Quadratic term beta = 0.5*Q is locked at 40 kHz:
beta = 1.5*alpha

# Drive term gamma = 0.5*gamma_n*B_1Q can be tuned from 0 up to 30 kHz:
gamma = .05*alpha

# AC B-field frequency in the rotating frame (in MHz):
f_rot = 1.5*alpha


### CREATE SPIN OPERATORS ###
[Ix, Iy, Iz, Ip, Im] = CreateSpinOperators(I)


### PERTURBATION OPERATOR ###
phiW = pi/40
W = linalg.expm(-1j*phiW*(Iy))


### ASSIGN FREQUENCY DEPENDING ON CHOICE OF FRAME ###
if statusFrame == 'rotating':
	f = f_rot
elif statusFrame == 'lab':
	f = f_lab
else:
	raise Exception("Please enter either \'lab\' or \'rotating\'.")
 

### SIMULATION PARAMETERS ###
# Time evolution starting point:
tStart = 0

# Time evolution ending point (followed by W and then time reversal):
tEnd = 60*(1/f)

# Number of Floquet periods in the chosen time evolution:
NPeriods = int((tEnd-tStart)*f)


### INITIALIZATION OF FLOQUET OPERATOR ###
# Number of segments in approximation of Floquet operator:
NFloquet = 1000

# Create the Floquet object:
FOp = None
U = None

if statusFrame == 'lab':
	args = [f,gamma_n,A,Q,B0,B1,Ix,Iy,Iz]
	FOp = FloquetOperator(ChaosHamiltonian, tStart, f, NFloquet, args)
elif statusFrame == 'rotating':
	args = [f,alpha,beta,gamma,Ix,Iy,Iz,I]
	FOp = FloquetOperator(ChaosHamiltonianRWA, tStart, f, NFloquet, args)
else:
	raise Exception("Please enter either \'lab\' or \'rotating\'.")

# Initialize time evolution operator U and Ud:
U = FOp.GetTimeEvolutionOperator(tStart)
Ud = U.T.conj()


### AZIMUTHAL AND POLAR ANGLES ###
# Azimuthal angle phi:
phiStart = -pi
phiEnd = phiStart + 2*pi
phiN = 80

# Polar angle theta:
thetaStart = 0
thetaEnd = thetaStart + pi
thetaN = 40

# Corresponding linspaces:
phiLin = [phiStart + (phiEnd-phiStart)/(phiN - 1)*i for i in range(phiN)]
thetaLin = [thetaStart + (thetaEnd-thetaStart)/(thetaN-1)*i for i in range(thetaN)]


### LOSCHMIDT ECHO SIMULATION LOOP ###
# We will store the Loschmidt echo in the variable LE:
LE = np.zeros((phiN,thetaN,NPeriods+1))

tic = time.time()
print('- - - - - - - - - -')
print('Now starting the Loschmidt echo calculation:')
print('- - - - - - - - - -')

for ti in range(NPeriods+1):
	t = tStart + ti*(1/f)
	tPrior = tStart + (ti-1)*(1/f)

	if t > tStart:
		# Update the time evolution operator if t > tStart:
		U = FOp.GetTimeEvolutionOperator(t,tPrior,U)
		Ud = U.T.conj()

		# Also broadcast info pertaining to the simulation progress:
		if ti%5 == 0:
			toc = time.time() - tic
			if toc > 120:
				print('\tTime step t = {ti}/{tTotal} reached after {tocM:.0f} minutes and {toc:.1f} seconds.'.format(ti=ti, tTotal = NPeriods+1, tocM=int(toc/60), toc = toc - 60*(int(toc/60))))
			elif toc > 60:
				print('\tTime step t = {ti}/{tTotal} reached after {tocM:.0f} minute and {toc:.1f} seconds.'.format(ti=ti, tTotal = NPeriods+1, tocM=int(toc/60), toc = toc - 60*(int(toc/60))))
			else:
				print('\tTime step t = {ti}/{tTotal} reached after {toc:.1f} seconds.'.format(ti=ti, tTotal = NPeriods+1, toc=toc))


	# We sample initial states according to the phiLen and thetaLen distributions:
	Psi_ini = np.zeros((N,1))
	Psi_ini[0] = 1

	for i in range(phiN):
		phi = phiLin[i]
		for j in range(thetaN):
			theta = thetaLin[j]

			# Create the spin-coherent state |phi,theta>:
			RotMat = linalg.expm(-1j*theta*(Ix*math.sin(phi) - Iy*math.cos(phi)))
			Psi0 = RotMat@Psi_ini

			# Calculate the Loschmidt echo:
			value = abs((Psi0.T.conj())@Ud@W@U@Psi0)**2
			LE[i,j,ti] = value

#Save the LE data to a file:
if statusFrame == 'lab':
	systemparameters = (I,gamma_n,A,Q,B0,B1,f)
elif statusFrame == 'rotating':
	systemparameters = (I,alpha,beta,gamma,f)
else:
	raise Exception("Please enter either \'lab\' or \'rotating\'.")
supportargs = (f,statusFrame,tStart,tEnd,phiLin,thetaLin,systemparameters)
np.save('LEdata',(LE,supportargs))

# Display some info to let the user know when the Loschmidt echo calculation ended:
toc = time.time()-tic
print('- - - - - - - - - -')
if toc > 120:
	print('Loschmidt echo calculation finished after {m:.0f} minutes and {s:.1f} second(s).'.format(m=int(toc/60),s=toc-60*int(toc/60)))
elif toc > 60:
	print('Loschmidt echo calculation finished after {m:.0f} minute and {s:.1f} second(s).'.format(m=int(toc/60),s=toc-60*int(toc/60)))
else:
	print('Loschmidt echo calculation finished after {s:.1f} second(s).'.format(s=toc-60*int(toc/60)))
print('- - - - - - - - - -')


### LOSCHMIDT ECHO SINGLE TRAJECTORY LOOP ###
# We now want to follow a single trajectory during the Loschmidt echo!
# We follow the state starting at |phi0,theta0>:
phi0 = 0*pi
theta0 = 0.4*pi
Psi_ini = np.zeros((N,1))
Psi_ini[0] = 1
Psi0 = linalg.expm(-1j*theta0*(Ix*math.sin(phi0)-Iy*math.cos(phi0)))@Psi_ini

# We will store the projection of the resulting state onto spin-coherent states in the variable LEtrajectory:
LEtrajectory = np.zeros((phiN,thetaN,2*NPeriods+2), dtype=np.complex_)

# Initialize time evolution operator U and a storage for Us:
U = FOp.GetTimeEvolutionOperator(tStart)
U_storage = [U]

tic = time.time()
print(' ')
print('- - - - - - - - - -')
print('Now calculating a single Loschmidt echo trajectory:')
print('- - - - - - - - - -')

# Forward time evolution:
for ti in range(NPeriods+1):
	t = tStart + ti*(1/f)
	tPrior = tStart + (ti-1)*(1/f)

	if t > tStart:
		# Update the time evolution operator if t > tStart:
		U = FOp.GetTimeEvolutionOperator(t,tPrior,U)
		U_storage.append(U)

		# Also broadcast info pertaining to the simulation progress:
		if ti%5 == 0:
			toc = time.time() - tic
			if toc > 60:
				print('\tTime step t = {ti}/{tTotal} reached after {tocM:.0f} minutes and {toc:.1f} seconds.'.format(ti=ti, tTotal = 2*NPeriods+2, tocM=int(toc/60), toc = toc - 60*(int(toc/60))))
			elif toc > 60:
				print('\tTime step t = {ti}/{tTotal} reached after {tocM:.0f} minute and {toc:.1f} seconds.'.format(ti=ti, tTotal = 2*NPeriods+2, tocM=int(toc/60), toc = toc - 60*(int(toc/60))))
			else:
				print('\tTime step t = {ti}/{tTotal} reached after {toc:.1f} seconds.'.format(ti=ti, tTotal = 2*NPeriods+2, toc=toc))

	# Do time evolution of the initial state:
	Psi_t = U@Psi0

	# We now do the Husimi-map of Psi_t:
	for i in range(phiN):
		phi = phiLin[i]

		for j in range(thetaN):
			theta = thetaLin[j]

			# Create the spin-coherent state |phi,theta>:
			RotMat = linalg.expm(-1j*theta*(Ix*math.sin(phi) - Iy*math.cos(phi)))
			Psi_SCS = RotMat@Psi_ini

			# Calculate the overlap |<phi,theta|psi_t>|^2:
			value = abs((Psi_SCS.T.conj())@Psi_t)**2
			LEtrajectory[i,j,ti] = value


# Apply the perturbation:
Psi0 = W@Psi_t

# We only stored the Us going forward and now need to convert them to Us going backward in time.
Ud_storage = [np.identity(N, dtype=np.complex_)]
for i in range(1,len(U_storage)):
	U = U_storage[-1]@(U_storage[-1-i].T.conj())
	Ud_storage.append(U.T.conj())

# Backward time evolution:
for ti in range(NPeriods+1):
	t = tEnd - ti*(1/f)

	Ud = Ud_storage[ti]

	# Also broadcast info pertaining to the simulation progress:
	if (ti+1)%5 == 0 and t < tEnd:
		toc = time.time() - tic
		if toc > 60:
			print('\tTime step t = {ti}/{tTotal} reached after {tocM:.0f} minutes and {toc:.1f} seconds.'.format(ti=ti+NPeriods+1, tTotal = 2*NPeriods+2, tocM=int(toc/60), toc = toc - 60*(int(toc/60))))
		elif toc > 60:
			print('\tTime step t = {ti}/{tTotal} reached after {tocM:.0f} minute and {toc:.1f} seconds.'.format(ti=ti+NPeriods+1, tTotal = 2*NPeriods+2, tocM=int(toc/60), toc = toc - 60*(int(toc/60))))
		else:
			print('\tTime step t = {ti}/{tTotal} reached after {toc:.1f} seconds.'.format(ti=ti+NPeriods+1, tTotal = 2*NPeriods+2, toc=toc))

	# Do backward time evolution of the initial state:
	Psi_t = Ud@Psi0

	# We now do the Husimi-map of Psi_t:
	for i in range(phiN):
		phi = phiLin[i]

		for j in range(thetaN):
			theta = thetaLin[j]

			# Create the spin-coherent state |phi,theta>:
			RotMat = linalg.expm(-1j*theta*(Ix*math.sin(phi) - Iy*math.cos(phi)))
			Psi_SCS = RotMat@Psi_ini

			# Calculate the overlap |<phi,theta|psi_t>|^2:
			value = abs((Psi_SCS.T.conj())@Psi_t)**2
			LEtrajectory[i,j,NPeriods+1+ti] = value

#Save the LEtrajectory data to a file:
np.save('LEtrajectory',(LEtrajectory,supportargs))

# Display some info to let the user know when the program ended.
toc = time.time()-tic
print('- - - - - - - - - -')
if toc > 120:
	print('Loschmidt echo trajectory calculation finished after {m:.0f} minutes and {s:.1f} second(s).'.format(m=int(toc/60),s=toc-60*int(toc/60)))
elif toc > 60:
	print('Loschmidt echo trajectory calculation finished after {m:.0f} minute and {s:.1f} second(s).'.format(m=int(toc/60),s=toc-60*int(toc/60)))
else:
	print('Loschmidt echo trajectory calculation finished after {s:.1f} second(s).'.format(s=toc-60*int(toc/60)))
print('- - - - - - - - - -')

