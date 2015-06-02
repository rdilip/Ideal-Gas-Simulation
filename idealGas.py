from molecule import molecule
import random
from math import sqrt, pow, exp
import numpy as np
class idealGas:
	# CONSTANTS
	def __init__(self, n, width, length, height):
		gas = []
		for i in xrange(0, n):
			gas.append(molecule(random.random() * width, \
				random.random() * length, random.random() \
				* height))
		"""
		for i in xrange(0, n):
			gas.append(molecule(width / 2 + random.random() * 20,\
			length/2 + random.random() * 20,\
			height/2 + random.random() * 20))
		"""
		#self.gas = gas
		self.gas = np.array(gas)
		self.width = width
		self.length = length
		self.height = height
		self.n = n

	def __str__(self):
		retStr = ""
		for i in xrange(0, len(self.gas)):
			retStr += ("Molecule %d at position (%.3f, %.3f, %.3f)\n" % (i, \
				self.gas[i].x, self.gas[i].y, self.gas[i].z))	
		return (retStr)

	def distance(self, i, j):
		x_i, y_i, z_i = self.gas[i].pos
		x_j, y_j, z_j = self.gas[j].pos
		dist = (x_i - x_j)**2 + (y_i - y_j)**2 + (z_i - z_j)**2
		return(dist**(0.5))

	def potEnergy(self, k):
		""" Returns the potential energy of the kth molecule due to all 
		other molecules in the system
		"""

		k_b = 9.36E-18
		A, B, E = 1, 1, 1
		
		gas = self.gas
		potEnergy = 0
		for i in xrange(0, len(gas)):
			if i != k:
				r_ki = self.distance(k, i)
				potEnergy += (1.0 / float(r_ki))*4*E*((float(A) / float(r_ki))**6 - (float(B) / float(r_ki)**3))
		return potEnergy
		
	def totEnergy(self):
		""" Returns the total potential energy of all molecules """

		gas = self.gas
		totEnergy = 0
		for i in xrange(0, len(gas)):
			totEnergy += self.potEnergy(i)
		return totEnergy / 2

	def posChange(self, temperature, i, dr, result):
		""" Changes the position of the kth molecule, subject to an
		acceptance probability based on the potential energy change
		Returns the indices and new positions of all changed molecules
		Also prints stuff
		"""
		k_b = 9.36E-18
		A, B, E = 1, 1, 1
		
		gas = self.gas


		U_i = self.totEnergy()

		dx, dy, dz = (-1.0 + 2 * random.random()) * dr,\
			(-1.0 + 2 * random.random()) * dr,\
			(-1.0 + 2 * random.random()) * dr 
		x_i, y_i, z_i = gas[i].x, gas[i].y, gas[i].z
		gas[i].x = (x_i + dx) 
		gas[i].y = (y_i + dy)
		gas[i].z = (z_i + dz) 
		gas[i].pos = (gas[i].x, gas[i].y, gas[i].z)

		U_f = self.totEnergy()
		if (U_f - U_i <= float(0.0)): 
			# accept it
			result.append((i, gas[i].x, gas[i].y, gas[i].z))
		else:
			prob = exp((-1.0 * (U_f - U_i)) / (k_b * temperature))
			if (random.random() < prob):
				# accept it
				result.append((i, gas[i].x, gas[i].y, gas[i].z))
			else:
				gas[i].x = x_i
				gas[i].y = y_i
				gas[i].z = z_i
		self.gas[i].pos = (gas[i].x, gas[i].y, gas[i].z)
		self.gas = gas

	def posUpdate(self, temperature, dr):
		""" Uses posChange() function to update all molecules in the system """
		gas = self.gas
		result = []
		#i = int(random.random() * len(gas))	
		for i in xrange(0, self.n):
			self.posChange(temperature, i, dr, result)
		return (result)

	def simulate(self, temperature, dr, trials):
		count = 0
		for i in xrange(0, trials):
			count += len(self.posUpdate(temperature, dr))
		return float(count) / float(self.n * trials)

