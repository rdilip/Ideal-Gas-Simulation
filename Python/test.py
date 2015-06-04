import Tkinter as tk
import random
import time
class App(tk.Tk):
	def __init__(self):
		tk.Tk.__init__(self)
		size = 500
		self.canvas = tk.Canvas(self, height=size, width= size,\
			borderwidth=0, highlightthickness=0)
		self.canvas.pack(side="top", fill="both", expand="true")
		
		self.cellwidth= 0.01
		self.cellheight= 0.01
		x = 250
		y = 300
		self.oval = {}
		self.oval[x, y] = self.canvas.create_oval(x+3, y+3, x-3, y-3,\
				fill="blue", tags = "oval")

		"""
		for i in xrange(1, 50):
			x = i * self.cellwidth
			y = i * self.cellheight
			self.oval[x, y] = self.canvas.create_oval(x+3, y+3, x-3, y-3,\
				fill="blue", tags = "oval")
"""
	
	def sim(self):
		dx = 1
		dy = 1
		self.canvas.move("oval", dx, dy)
		
		self.after(1, lambda: self.sim())

if __name__ == '__main__':
	app = App()
	app.mainloop()
