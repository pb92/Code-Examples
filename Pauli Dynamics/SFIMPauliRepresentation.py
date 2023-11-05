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
N = 4

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

# Flip up->down rate:
GammaUD = 0;

# Flip down->up rate:
GammaDU = 0;

# Dephasing rate:
GammaZ = 0;

# Pauli weight cutoff (controls the precision of the Pauli string truncation):
weightCutoff = N


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



# Create matrix with site combinations for every weight for future reference:
NumberOfConfigurations = 0;
for w in range(N):
    NumberOfConfigurations += math.comb(N,w+1)
    
SiteConfigurationMatrix = np.zeros((NumberOfConfigurations,N))

index = 0

for w in range(N):
    indexVector = np.zeros((1,N))
    
    
    for i in range(2**N-1):
        
        indexVector[0,N-1] = indexVector[0,N-1]+1
        
        for j in range(N-1,0,-1):
            if indexVector[0,j] > 1:
                indexVector[0,j] = 0
                indexVector[0,j-1] = indexVector[0,j-1]+1
                
        if np.sum(indexVector) == w:
            index += 1
            
            for j in range(N):
                SiteConfigurationMatrix[index,j] = indexVector[0,j]
                


    
### MATRIX REPRESENTATION OF THE LINDBLAD ME ###
# We now need to find the MR of the LME using the weight-sorted Pauli list. We
# start by finding the MR of the system Hamiltonian in the truncated Pauli 
# Hilbert space:    
HamiltonianList = createSFIMHamiltonian(N,J,Bperp,Bpara)
    
# Create a storage variable for the MR of the Hamiltonian:
H = np.zeros((len(PauliListSorted),len(PauliListSorted)),dtype = "complex_")
Htest = np.zeros((len(PauliListSorted),len(PauliListSorted)),dtype = "complex_")
    
# Run over every single Pauli in the truncated space:
for j in range(len(PauliListSorted)):
    
    # Pauli string:
    PauliIn = PauliListSorted[j]    
    
    # Run over all terms in the Hamiltonian:
    for h in range(len(HamiltonianList)):
        
        # Read out the hth term in the Hamiltonian:        
        HamiltonianTerm = HamiltonianList[h]
        
        # Pauli string of the term:
        PauliH = HamiltonianTerm[1]
        
        # Apply the hth term of the Hamiltonian to the jth Pauli string:        
        PauliHIn = PauliSymProdFull(PauliH,PauliIn)
        PauliHInCoeff = PauliHIn[0]
        PauliHInString = PauliHIn[1]
        
        # Calculate the resulting weight of applying the hth term to the jth
        # Pauli string:
        weightHIn = PauliWeight(PauliHInString)
        
        # If the resulting Pauli weight is greater than the truncation cutoff,
        # we ignore it:
        if weightHIn > weightCutoff:
            continue
        
        # If it is smaller than the truncation cutoff, we find the resulting
        # Pauli string and store the coefficient in the MR of H:
        else:
            # Run over all Pauli strings that have the same weight as the
            # resulting Pauli string:
            for i in range(int(PauliListWeightIndex[weightHIn]),int(PauliListWeightIndex[weightHIn+1])):
                
                PauliOut = PauliListSorted[i]
                
                if PauliStringIsIdentical(PauliOut,PauliHInString):
                    H[i,j] += HamiltonianTerm[0]*PauliHInCoeff                    
                    continue
                
                
# With the MR of the Hamiltonian found, we can construct the MR of the unitary
# part of the Lindblad master equation in the following way:
M = -1j*(H - H.T)


# The decoherence is now taken into account by adding the jump operators in the
# Lindblad master equation style. The dephasing noise is the easiest to add due
# to the jump operator being one of the Pauli operators and being diagonal in
# the basis:
if GammaZ != 0:

    for i in range(len(PauliListSorted)):
        
        Pauli = PauliListSorted[i]
        
        coefficient = 0
        
        for siteIndex in range(N):
            
            if Pauli[siteIndex] != 'Z' and Pauli[siteIndex] != 'I':
                coefficient -= 2
                
        coefficient *= GammaZ
        
        M[i,i] += coefficient


# The dissipative effects require more work:
if GammaUD != 0 or GammaDU != 0:
    
    # First we add the diagonal terms that appear as soon as GammaUD or GammaDU
    # are non-zero:
    for i in range(len(PauliListSorted)):
        
        # Pauli string:
        PauliIn = PauliListSorted[i]  
        
        coefficient = 0.0
        
        for siteIndex in range(N):
            
            if PauliIn[siteIndex] == 'X':
                coefficient -= 0.5
                
            if PauliIn[siteIndex] == 'Y':
                coefficient -= 0.5            
            
            if PauliIn[siteIndex] == 'Z':
                coefficient -= 1
                
        coefficient *= (GammaUD+GammaDU)
        
        M[i,i] += coefficient
    
    
    # If GammaUD isn't equal to GammaDU, there will be non-diagonal terms that
    # we need to add. These have to be handled case-by-case:
    if GammaUD != GammaDU:
        
        for j in range(len(PauliListSorted)):
            
            # Pauli string:
            PauliIn = PauliListSorted[j]
            
            for siteIndex in range(N):
                
                if PauliIn[siteIndex] == 'I':
                    PauliRes = PauliIn.copy()
                    
                    PauliRes[siteIndex] = 'Z'
                    PauliResCoeff = (GammaUD-GammaDU)
                        
                    # Calculate the weight of PauliRes:
                    weightPauliRes = PauliWeight(PauliRes)
                    
                    # If the resulting Pauli weight is greater than the 
                    # truncation cutoff we ignore it:
                    if weightPauliRes > weightCutoff:
                        continue
                    
                    # If it is smaller than the truncation cutoff, we store the
                    # coefficient in the appropriate entry of the MR:
                    else:
                        # Run over all Pauli strings that have the same weight
                        # as the resulting Pauli string:
                        for i in range(int(PauliListWeightIndex[weightPauliRes]),int(PauliListWeightIndex[weightPauliRes+1])):
                
                            PauliOut = PauliListSorted[i]
                
                            if PauliStringIsIdentical(PauliOut,PauliRes):
                                M[i,j] += PauliResCoeff                 
                                continue
                    
        
    
# We have now obtained the MR of the LME in the truncated Pauli basis.

fileName = "N" + str(N) + "_C" + str(weightCutoff) + "_TFIMcritical_GammaUD0_GammaDU0_GammaZ0.npy"
#np.save(fileName,M)