import Tkinter as tk
import csv
import random
import time
class App(tk.Tk):
	def __init__(self):
		# scale and scale1 differ for visual purposes. scale1 refers to 
		# the visual size of the molecules.
		scale = 500/125
		scale1 = 2

		tk.Tk.__init__(self)

		self.canvas = tk.Canvas(self, height=125 * scale, width=125 * scale,\
			borderwidth=0, highlightthickness=0)

		self.canvas.pack(side="top", fill="both", expand="true")
		
		self.cellwidth= 1 * scale
		self.cellheight= 1 * scale

		self.oval = {}
		positions = list(csv.reader(open('final_config.txt', 'rb'), delimiter='\t'))
		for i in positions:
			for j in xrange(3):
				i[j] = float(i[j])

		
		for i in xrange(0, 50):
			x = positions[i][0]
			y = positions[i][1]

			self.oval[x, y] = self.canvas.create_oval(x+2*scale1\
				, y+2*scale1, x-2*scale1, y-2*scale1,fill="blue",\
				tags = "oval" + str(i))


		self.sim()
	
	def sim(self):
		scale = 500/125
		scale1 = 2


		# block refers to the size of each step. Each posUpdate 
		# runs through all n molecules. The size of block refers
		# to the number of posUpdate functions called per visual update

		# delay refers to the time delay between cycles, measured 
		# in milliseconds
		delay = 100

		positions = list(csv.reader(open('final_config.txt', 'rb'), delimiter='\t'))
		for i in positions:
			for j in xrange(3):
				i[j] = float(i[j])

		
		for i in xrange(0, 50):
			x = positions[i][0]
			y = positions[i][1]

			self.oval[x, y] = self.canvas.create_oval(x+2*scale1\
				, y+2*scale1, x-2*scale1, y-2*scale1,fill="blue",\
				tags = "oval" + str(i))


		self.after(delay, lambda: self.sim())
		
if __name__ == '__main__':
	app = App()
	app.mainloop()
