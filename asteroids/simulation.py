# simulation.py

# This module defines the functions required to run the animation and then 
# using these functions the simulation is run.

print("loading...")

from SolarSystemRotatingFrame import *
from visual import *
import threading, thread
import os, sys

################################################################################
###	Function definitions
################################################################################
	
def intialisebodies(system, Trails):
	DictOfBodies = {}
	DictOfTrails = {}
	
	global data

	i=0
	for body in system.bodies:
		name = body.name
		
		if name == 'Sun':
			DictOfBodies[name]=sphere(pos=body.position,\
					radius=0.5, color=color.yellow)
		
		elif name == 'Jupiter':
			DictOfBodies[name]=sphere(pos=body.position,\
					radius=0.2, color=color.red)
			
		else:
			if i < system.N_asteroids/2:
				DictOfBodies[name]=sphere(pos=body.position,\
					radius=0.05, color=color.white)
				# asteroids at L4
				i+=1
			else:
				DictOfBodies[name]=sphere(pos=body.position,\
					radius=0.05, color=color.blue)
				# asteroids at L5
				i+=1
				
		if Trails:
			DictOfTrails[name]=curve(color=DictOfBodies[name].color)
		
	return DictOfBodies, DictOfTrails
		
def visualisation(system, DictOfBodies, DictOfTrails, Trails, RFOR, \
	Trail_limit, Rate):
	
	scene.visible = True

	global data
	global simulation
	
	rotation_angle = 2*np.pi*system.dt/system.period
	# The rotation angle per time step. This variable is only used if
	# RFOR = False
	
	L4 = system.L4
	L5 = system.L5
	
	LPs = points(pos=[L4, L5], size=2, color=color.green)
	# Create points showing where the lagrange points are

	i=0
	while True:
		
		rate(Rate)

		j = -2	# changed from 0 to -2 to accomodate the sun and planet
		for body in system.bodies:
			
			name = body.name
			
			loop = True
			while loop==True:
				
				current_positions = []
				
				try:
					if RFOR:
						
						if j==-2: # the sun
							pass
						elif j==-1: # the planet
							pass
						else:
							DictOfBodies[name].pos = list(data[i,j,0,:])	

					else:
						rotate = system.rotationmatrix(rotation_angle*i)
						if j==-2: # the sun
							DictOfBodies[name].pos = np.dot(rotate, system.sunposition)
							
						elif j==-1: # the planet
							DictOfBodies[name].pos = np.dot(rotate, system.planetposition)
							
						else:
							DictOfBodies[name].pos = np.dot(rotate, list(data[i,j,0,:]))
					loop = False
					
				except IndexError:
					print("loading...")
					print("Tip: try with fewer asteroids \
and with trails off to make the simulation run faster.")
					time.sleep(1)
				
			if Trails:
				try:
					DictOfTrails[name].append(pos = DictOfBodies[name].pos)
				except KeyError:
					DictOfTrails[name] = curve(color=DictOfBodies[name].color)
					DictOfTrails[name].pos = DictOfBodies[name].pos
				
				if Trail_limit != None:
				#This places a limit on the trail length
					if len(list(DictOfTrails[name].pos)) > Trail_limit:
						l = list(DictOfTrails[name].pos)
						l.pop(0)
						DictOfTrails[name].pos = l

				
			if j>=0:
				data[i,j,:,:] = np.NaN 	
			# remove data that has already been used.
			j+=1
		
		if RFOR==False: 
			# Rotate positions of lagrange points if necessary
			rotate = system.rotationmatrix(rotation_angle*i)
			L4 = np.dot(rotate, system.L4)
			L5 = np.dot(rotate, system.L5)
			
			LPs.pos = [L4, L5]	
			# adjust positions of green points showing the LPs
		
		i+=1
				
		if scene.kb.keys: # event waiting to be processed?
			s = scene.kb.getkey()
			print(s)
		
		try:
			if s == 's': 
				# to centre the camera on the sun, press 's'
				scene.center = DictOfBodies['Sun'].pos	
			elif s == 'j': 
				# to centre the camera on jupiter,press 'j'
				scene.center = DictOfBodies['Jupiter'].pos
			elif s == '4': 
				# to centre on L4, press '4'
				scene.center = L4
			elif s == '5': 
				# to centre on L5, press '5'
				scene.center = L5
			
			if s == 'o': # turn trails off by pressing 'o'
				Trails = False
			elif s == 'i': # turn trails on by pressing 'i'
				Trails = True
			
			if s == 'd': # delete trails
				for body in system.bodies:
					name = body.name
					DictOfTrails[name].pos = []
					
			if s == 'ctrl+c': # make a way of exiting the simulation
				print("The program is exiting.")
				simulation = False
				os._exit(0)
				
		except NameError:
			pass
		

def computation(system):
	global data
	global simulation
	
	while simulation==True:
		t, results = system.results()
		
		data = np.concatenate((data, results))
		

################################################################################
###	System Parameters
################################################################################

Jmass = 0.001		# The mass of Jupiter in units of the solar mass
R = 5.2			# The mean distance between the Sun and Jupiter in AU
dt = 0.01		# Time step in Earth years
simulationtime = 10	# simulation time step in Earth years
N_asteroids = 2		# Total number of asteroids, half at each L.P.
offset = [0,0,0]	# The change in the initial position of the asteroids 
			# from the L.P., in cylindrical polars, in AU.
randomamplitude = 0.05

RFOR = True		# RFOF = Rotating frame of reference
# If True, then the simulation will be displayed in the frame of reference that 
# rotates with the sun and planet. If false, then the simulation is viewed in 
# the solar system inertial frame.

Trail_limit = None	# The time limit of points included in the trails

Rate = 1000		# The simulation rate. If the simulation is jittery, try
			# reducing the simulation rate

################################################################################
### 	The main simulation starts here
################################################################################

# first, create an instance of SolarSystem with the desired parameter values.
system = SolarSystem(Jmass, R, dt, simulationtime, N_asteroids, offset, \
	randomamplitude)	 
system.setup()

scene.visible = False
Trails = False 			# Trails off by default
	    
DictOfBodies, DictOfTrails = intialisebodies(system, Trails)

global data
global simulation

t, data = system.results()	# pre-load some data

simulation = True

print("The simulation is starting.")

# Create two threads. thread 1 to display the simulation and thread 2 to compute
# the future positions of the asteroids.
thread1 = threading.Thread(target=visualisation, args=(system, DictOfBodies, \
	DictOfTrails, Trails, RFOR, Trail_limit, Rate))
thread2 = threading.Thread(target=computation, args=(system,))

thread1.start()
thread2.start()

