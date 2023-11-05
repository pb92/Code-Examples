### PREAMBLE ###
# Philip Daniel Blocher
# Center for Quantum Information and Control
# Department of Physics and Astronomy, University of New Mexico
# Email: blocher@unm.edu (current email available on Google Scholar)


### PROGRAM DESCRIPTION ###
# This program creates the matrix representation of the transverse field Ising
# model (TFIM) in the Pauli operator basis.


### IMPORTED LIBRARIES AND GLOBAL VARIABLES ###
import math, cmath
import random, time
import numpy as np
from scipy import linalg
from scipy import special
import matplotlib.pyplot as plt
import itertools

pi = math.pi


### FUNCTION DEFINITIONS ###

# Hermitian conjugate of matrix:
def herm(A):
    return np.conjugate(A.T)

# Single-site Pauli algebra:
def PauliSymProdSingleSite(A,B):
    
    # Check for errors:
    if A  != 'I' and A != 'X' and A != 'Y' and A != 'Z':
        print('ERROR: PauliSymProdSingleSite received invalid input A')
        return None
    
    if B  != 'I' and B != 'X' and B != 'Y' and B != 'Z':
        print('ERROR: PauliSymProdSingleSite received invalid input B')
        return None
    
    # Return the result of AB according to the Pauli operator algebra:
    if A == 'I':
        return (1, B)
    elif B == 'I':
        return (1, A)
    elif A == 'X':
        if B == 'X':
            return (1, 'I')
        elif B == 'Y':
            return (1j, 'Z')
        else:
            return (-1j, 'Y')
    elif A == 'Y':
        if B == 'X':
            return (-1j, 'Z')
        elif B == 'Y':
            return (1, 'I')
        else:
            return (1j, 'X')
    else:
        if B == 'X':
            return (1j, 'Y')
        elif B == 'Y':
            return (-1j, 'X')
        else:
            return (1, 'I')
        

# Pauli algebra for all sites:
def PauliSymProdFull(Avec,Bvec):
    
    # Check for errors in the length of the Pauli strings:    
    N = len(Avec)
    if len(Bvec) != N:
        print('ERROR: PauliSymProdFull was given vectors of different lenghts')
        return None
    
    PauliString = []
    PauliStringCoefficient = 1
    for index in range(N):
        
        result = PauliSymProdSingleSite(Avec[index],Bvec[index])
        PauliString.append(result[1])
        PauliStringCoefficient *= result[0]
        
    return (PauliStringCoefficient, PauliString)


# Function that returns the weight of a given Pauli operator:
def PauliWeight(Avec):
    weight = 0
    
    N = len(Avec)
    
    for index in range(N):
        
        if Avec[index] != 'I':
            weight += 1;
        
    return weight


# Function that creates the transverse field Ising model as a list of Pauli
# strings with appropriate coefficients:
def createSFIMHamiltonian(N,J,Bperp,Bpara):
    Hamiltonian = []
    
    # Nearest-neighbor interaction terms:
    if J != 0:
        for siteIndex in range(N-1):
            PauliString = []
            
            for index in range(siteIndex):
                PauliString.append('I')
                
            PauliString.append('Z')
            PauliString.append('Z')
            
            for index in range(siteIndex+2,N):
                PauliString.append('I')
                
            Hamiltonian.append((-J,PauliString))
    
    # Transverse field terms:
    if Bperp != 0:
        for siteIndex in range(N):
            PauliString = []
            
            for index in range(siteIndex):
                PauliString.append('I')
                
            PauliString.append('X')
            
            for index in range(siteIndex+1,N):
                PauliString.append('I')
                
            Hamiltonian.append((-Bperp,PauliString))
        
    # Parallel field terms:
    if Bpara != 0:
        for siteIndex in range(N):
            PauliString = []
            
            for index in range(siteIndex):
                PauliString.append('I')
                
            PauliString.append('Z')
            
            for index in range(siteIndex+1,N):
                PauliString.append('I')
                
            Hamiltonian.append((-Bpara,PauliString))
        
    return Hamiltonian


# Function that returns whether two Pauli strings are identical:
def PauliStringIsIdentical(Avec,Bvec):
    # Check whether the two Pauli strings are of the same length:
    N = len(Avec)
    if len(Bvec) != N:
        return False
    
    # Check whether the two Pauli strings are of the same weight:
    if PauliWeight(Avec) != PauliWeight(Bvec):
        return False
    
    # Given that the two Pauli strings are of the same length and weight, we 
    # now compare the two Pauli strings element-wise:
    for siteIndex in range(N):
        if Avec[siteIndex] != Bvec[siteIndex]:
            return False
    
    # As all elements turned out to be identical, we will return True:
    return True






# # # MAIN CODE STARTS HERE # # #

### SYSTEM PARAMETERS ###
# Number of sites:
N = 6

# Number of on-site degrees of freedom (for qubits d = 2):
d = 2

# Operator Hilbert space size:
D = d**(2*N)

# Nearest-neighbor spin interaction for the slanted field Ising model:
J = 1

# Transverse B-field amplitude:
Bperp = 1;

# Parallel B-field amplitude:
Bpara = 0;

# Depolarization rate:
Gamma = 0;

# Pauli weight cutoff (controls the precisiozn of the Pauli string truncation):
weightCutoff = 1


### Pauli matrices ###
Sx = np.matrix([[0,1],[1,0]])
Sy = np.matrix([[0,-1j],[1j,0]])
Sz = np.matrix([[1,0],[0,-1]])
Id = np.eye(2)


### TRUNCATED PAULI HILBERT SPACE ###
# Create list with Pauli strings (not ordered by weights at first):
PauliList = []

# Create index vector:
IndexVector = np.zeros(N)

# Run through each possible Pauli and append it to PauliList if its weight is
# below the Pauli weight cutoff:
for index in range(D):
    PauliString = []
    
    for siteIndex in range(N):
        if IndexVector[siteIndex] == 0:
            PauliString.append('I')
        elif IndexVector[siteIndex] == 1:
            PauliString.append('X')
        elif IndexVector[siteIndex] == 2:
            PauliString.append('Y')
        elif IndexVector[siteIndex] == 3:
            PauliString.append('Z')
    
    if PauliWeight(PauliString) <= weightCutoff:
        PauliList.append(PauliString)
    
    IndexVector[N-1] += 1
    
    for index2 in range(N-1,-1,-1):
        if IndexVector[index2] > 3 and index != D:
            IndexVector[index2] = 0
            IndexVector[index2-1] += 1
            
# Now create a list of Pauli strings sorted by weights:
PauliListSorted = []
PauliListWeightIndex = np.zeros((N+2,1))

for weightIndex in range(weightCutoff+1):
    
    if weightIndex > 0:
        PauliListWeightIndex[weightIndex+1] = PauliListWeightIndex[weightIndex]
    
    for index in range(len(PauliList)):
        
        weight = PauliWeight(PauliList[index])
        
        if weight == weightIndex:
            PauliListSorted.append(PauliList[index])
            PauliListWeightIndex[weightIndex+1] += 1 
    
    
### MATRIX REPRESENTATION OF THE LINDBLAD ME ###
# Load the MR of the LME appropriate for the chosen truncation level:
fileName = "N" + str(N) + "_C" + str(weightCutoff) + "_TFIMcritical_GammaUD0_GammaDU0_GammaZ0.npy"

M = np.load(fileName)


### CREATE INITIAL STATE ###
# We want to create an initial product state in the (possibly truncated) Pauli
# operator basis. We will do (1+Y)^N to avoid eigenstates of H. To do so, we 
# need to assign the appropriate probability to all terms that contain only Zs.

rhoInitial = np.zeros((len(PauliListSorted),1))

for PauliIndex in range(len(PauliListSorted)):
    
    PauliOperator = PauliListSorted[PauliIndex]
    
    statusPassed = True
    
    for SiteIndex in range(N):
        
        if PauliOperator[SiteIndex] != 'I' and PauliOperator[SiteIndex] != 'Y':
            statusPassed = False
            continue
        
    if statusPassed == True:
        rhoInitial[PauliIndex] = 1/math.sqrt(D)
        
        

### TIME EVOLUTION ###
# Time evolution step parameters (in units of 1/J):
tStart = 0
tEnd = 5
dt = 1e-3
tN = int((tEnd-tStart)/dt)    
    
# Number of times where we store the state:
NStore = 100

# Time evolution storage variables:
rhoStorage = np.zeros((len(PauliListSorted),NStore+1),dtype = "complex_")
rhoStorage[:,0] = rhoInitial[:,0]

# Create the time-evolved density matrix:    
rho = rhoInitial

# Now evolve rho forward in time:    
for tIndex in range(tN):
    
    # Calculate drho for a single time step:
    drho = M@rho*dt
    
    # Assign rho the updated value:    
    rho = rho + drho
    
    # Store the state according to the number of store times defined above:
    if (tIndex+1) % int(tN/NStore) == 0:
        rhoStorage[:,int((tIndex+1)/int(tN/NStore))] = rho[:,0]
        

# Save rhoStorage for later comparison:
fileName = "N" + str(N) + "_C" + str(weightCutoff) + "_TFIMcritical_GammaUD0_GammaDU0_GammaZ0_rhoStorage.npy"
np.save(fileName,rhoStorage)
        
    






