# Ideal Gas Simulation of Helium

This respository contains an ongoing attempt to simulate the Helium-Helium gas interactions using a Monte-Carlo implementation of the Metropolis algorithm. This algorithm works by moving each molecule in a random direction dr, then probabilitically accepting or rejecting the movement based on the Maxwell-Boltzmann distribution factor (exp(-dU / kB x T)). 

There are two implementations in here. One is a python implementation, which contains three files -- molecule.py, idealGas.py, and visual.py. Molecule is a class that simply contains position information about a molecule. IdealGas imports molecule and creates an array of molecules -- this is where most of the computation is done. Visual provides a visual simulation of the molecules using Tkinter. These files may contain some errors or computational goofs.

For all practical purposes, I am updating idealGas.c as a main file, which creates a struct molecule and a struct idealGas, then proceeds from there. This is much faster and can handle many more molecules, which is highly beneficial. The algorithms are more or less the same. For this simulation, I implement periodic boundary conditions, and create a radius of inclusion given by half the size of the enclosing box. 

Update 2015/6/2: The C code now accepts more values directly from the command line, several values are constants within the code itself. There are stilla few bugs which I am currently working on fixing. 

If you have any questions, don't email me, just figure it out yourself. 
