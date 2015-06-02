class molecule:
	def __init__(self, xpos, ypos, zpos):
		self.pos = (xpos, ypos, zpos)
		self.x = xpos
		self.y = ypos
		self.z = zpos

	def __str__(self):
		return ("Position: (%.3f, %.3f, %.3f)" % (self.x, self.y, self.z))

