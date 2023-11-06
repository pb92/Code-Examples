### Preamble ###
# Philip Daniel Blocher
# Quantum Optics Group
# Department of Physics and Astronomy
# Aarhus University, Denmark
# Email: philipblocher@gmail.com
# Last updated: May 2019

### Description ###
# This function creates the Sx, Sy, Sz, Sp, and Sm operators for a given spin S and returns these as numpy arrays.

import numpy
import math

def CreateSpinOperators(S):
	# Create a list over m_I values:
	mList = [-x/2 for x in range(-int(S*2),int(S*2)+2,2)]
	N = len(mList)

	# Iz matrix:
	Sz = numpy.diag(mList)

	# Ip matrix:
	Sp = numpy.zeros((N,N))
	for i in range(0,N-1):
		Sp[i,i+1] = math.sqrt(S*(S+1) - mList[i+1]*(mList[i+1]+1))

	# Im matrix:
	Sm = numpy.zeros((N,N))
	for i in range(0,N-1):
		Sm[i+1,i] = math.sqrt(S*(S+1) - mList[i]*(mList[i]-1))

	# Ix & Iy matrices:
	Sx = (Sp + Sm) / 2
	Sy = (Sp - Sm) / (2j)

	return [Sx, Sy, Sz, Sp, Sm]