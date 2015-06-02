import Tkinter as tk
import random
import time
from idealGas import idealGas
class App(tk.Tk):
	def __init__(self, num_molecules, size, temperature, dr):
		# scale and scale1 differ for visual purposes. scale1 refers to 
		# the visual size of the molecules.
		scale = 500/20
		scale1 = 10

		tk.Tk.__init__(self)

		self.canvas = tk.Canvas(self, height=size * scale, width=size * scale,\
			borderwidth=0, highlightthickness=0)

		self.canvas.pack(side="top", fill="both", expand="true")
		
		self.cellwidth= 1 * scale
		self.cellheight= 1 * scale

		self.oval = {}

		xenon = idealGas(num_molecules, size, size, size) 

		for i in xrange(0, xenon.n):
			x = xenon.gas[i].x * scale
			y = xenon.gas[i].y * scale

			self.oval[x, y] = self.canvas.create_oval(x+2*scale1\
				, y+2*scale1, x-2*scale1, y-2*scale1,fill="blue",\
				tags = "oval" + str(i))
		print xenon.simulate(temperature, dr, 1000)

		self.sim(temperature, dr, xenon, num_molecules, 0, scale, scale1)
	
	def sim(self, temperature, dr, gas, num_molecules, count, scale,\
	 	scale1):
		# block refers to the size of each step. Each posUpdate 
		# runs through all n molecules. The size of block refers
		# to the number of posUpdate functions called per visual update
		block = 1

		# delay refers to the time delay between cycles, measured 
		# in milliseconds
		delay = 100


		movements = gas.posUpdate(temperature, dr)	

		if (count % block == 0):
			for k in xrange(0, len(movements)):
				i, x, y, z = movements[k] 
				x, y, z = x * scale, y * scale, z * scale

				self.canvas.coords("oval" + str(i), x+0.7*scale1\
				, y+0.7*scale1, x-0.7*scale1, y-0.7*sscale1)
			#	self.canvas.move("oval" + str(i), dx * scale, dy * scale)
	#	print gas
		

		self.after(100, lambda: self.sim(temperature,\
			dr, gas, num_molecules, count + 1, scale, scale1))

		
if __name__ == '__main__':
	app = App(600, 100, 300, 2)
	app.mainloop()
