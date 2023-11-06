### Preamble ###
# Philip Daniel Blocher
# Quantum Optics Group
# Department of Physics and Astronomy
# Aarhus University, Denmark
# Email: philipblocher@gmail.com
# Last updated: May 2019

### Description ###
# This class creates the Floquet time evolution operator F for use in systems with periodic drive.

import math
import numpy as np
from scipy import linalg

class FloquetOperator:

	def __init__(self, Hamiltonian, tStart, f, NSegments, args):
		self.f = f
		self.tStart = tStart

		tPeriod = 1/self.f
		H = Hamiltonian(self.tStart,args)

		self.F = np.identity(H.shape[0])

		for i in range(NSegments):
			t = self.tStart+(i+1)*tPeriod/NSegments
			H = Hamiltonian(t,args)
			self.F = linalg.expm(-1j*H*(tPeriod/NSegments))@self.F


	def GetFloquetOperator(self):
		return self.F

	def GetTimeEvolutionOperator(self, t, tPrior = None, UPrior = None):
		if tPrior == None:
			tPrior = self.tStart
			UPrior = np.identity(self.F.shape[0])
			print('CreateFloquetOperator.py: No tPrior specified. Using default setting tPrior = tStart.')

		NPeriods = round((t-tPrior)*self.f)

		if NPeriods == 0:
			return UPrior
		else:
			return self.__MultiplyFloquetOperators(NPeriods)@UPrior

	def __MultiplyFloquetOperators(self,N):
		n = math.log2(N)
		n = int(n)
		remainder = N - 2**n

		F = np.copy(self.F)

		for i in range(n):
			F = F@F


		Fprime = np.identity(F.shape[0])
		if remainder > 1:
			Fprime = self.__MultiplyFloquetOperators(remainder)
		elif remainder == 1:
			Fprime = np.copy(self.F)


		return F@Fprime