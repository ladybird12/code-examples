# SolarSystemRotatingFrame.py

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random
import time

class SolarSystem:
	'''
	The SolarSystem class is the superclass that represents the whole 
	system. This class sets up the system parameters and the subclass 
	Body adds in the sun, planet and asteroids. In this program, the frame 
	of reference is a rotating one in which both the sun and the planet are 
	at fixed positions (remember that we are making the approximation that 
	the sun and planet undergo circular motion).
	'''
	asteroids = []
	# the list of asteroids, set to be empty for now.
	bodies = []
	# the list of all the bodies in the system, sun and planet included.
	
	def __init__(self, planetmass, R, dt, simulationtime, N_asteroids, \
		offset, randomamplitude):

		if N_asteroids%2 == 1: 	
		# ensure that N_asteroids is even by increasing N_asteroids by 
		# one if it is odd.
			N_asteroids += 1

		self.planetmass = planetmass
		self.R = R
		self.dt = dt
		self.simulationtime = simulationtime
		self.N_steps = simulationtime/dt
		self.N_asteroids = N_asteroids
		self.offset = offset
		# the initial displacement of the asteroids from the LPs in 
		# cylindrical polar form.
		
		self.randomamplitude = randomamplitude	
		# add in a random factor to the asteroid position
		
		self.sunmass = 1	# in units of the solar mass!
		self.period = (R**3/(self.sunmass+self.planetmass))**0.5
		# The period is derived from the condition that the net force on
		# an asteroid at either L4 or L5 must be zero.
		self.sunposition = np.array([-R*self.planetmass/(self.sunmass +\
			self.planetmass), 0, 0])	
		self.planetposition = np.array([R*self.sunmass/(self.sunmass + \
			self.planetmass), 0, 0])
		self.omega = np.array([0, 0, 2*np.pi/self.period])
		self.L4 = None
		self.L5 = None
	
	def setup(self):
		self.lagrangepoints()	
		# specify the positions of the lagrange points
		
		Body('Sun', self.sunposition, np.zeros(3), mass=self.sunmass)
		Body('Jupiter', self.planetposition, np.zeros(3), \
			mass=self.planetmass)
		
		L4 = self.L4	# L4 in cartesian coordinates
		L5 = self.L5	# L5 in cartesian coordinates
		
		L4cp = np.array([(L4[0]**2 + L4[1]**2)**0.5, \
		np.arctan(L4[1]/L4[0]), L4[2]])	# L4 in cylindrical polars (cp)
		L5cp = np.array([(L5[0]**2 + L5[1]**2)**0.5, \
		np.arctan(L5[1]/L5[0]), L5[2]])	# L5 in cylindrical polars		
		
		for j in range(int(self.N_asteroids/2)): # set up the Greeks
			offset = self.offset
			adjustedoffset = np.array([offset[0], \
			offset[1]/L4cp[0], offset[2]])
			poscp = L4cp + adjustedoffset # asteroid position in cp
			
			position = np.array([poscp[0]*np.cos(poscp[1]),poscp[0]\
			*np.sin(poscp[1]), poscp[2]]) + self.randomvector()
			# convert back to cartesian coordinates and add the
			# random vector
			
			Body('Greek{}'.format(j+1), position, np.zeros(3))
			
		for j in range(int(self.N_asteroids/2)): # set up the Trojans
			offset = self.offset
			adjustedoffset = np.array([offset[0], \
			-offset[1]/L5cp[0], offset[2]])
			# include minus sign to preserve symmetry in system
			poscp = L5cp + adjustedoffset # asteroid position in cp
			
			position = np.array([poscp[0]*np.cos(poscp[1]),poscp[0]\
			*np.sin(poscp[1]), poscp[2]]) + self.randomvector()
			# convert back to cartesian coordinates and add the
			# random vector
			
			Body('Trojan{}'.format(j+1), position, np.zeros(3))
				
		
	def rotationmatrix(self, angle):
		# a matrix which rotates a vector in the positive phi direction
		# (anticlockwise in the x-y plane when viewed from above) by
		# an angle 'angle'.
		return np.array([[np.cos(angle), -np.sin(angle), 0],\
		[np.sin(angle), np.cos(angle), 0],[0,0,1]])
		
	def lagrangepoints(self):
		angle = [np.pi/3, -np.pi/3]
		# the L4 and L5 lagrange points are located pi/3 ahead and pi/3
		# behind the planet in its orbit respectively. The sun, the 
		# planet and the lagrange points form equilateral triangles.
		
		suntoplanet = self.planetposition - self.sunposition	
		# the displacement vector from the sun to the planet
		
		self.L4 = np.dot(self.rotationmatrix(angle[0]), suntoplanet) + \
			self.sunposition
		self.L5 = np.dot(self.rotationmatrix(angle[1]), suntoplanet) + \
			self.sunposition
			
	def randomvector(self): 
	# generates a random vector in 3D with a length given by randomamplitude.
	# use spherical polar coordinates
		r = self.randomamplitude
		theta = random.uniform(0,np.pi)
		phi = random.uniform(0,2*np.pi)
		
		return np.array([r*np.sin(theta)*np.cos(phi), \
				r*np.sin(theta)*np.sin(phi), r*np.cos(theta)])
	
	def displacement(self, position1, position2):
		# displacement vector from position1 to position2.
		return position2 - position1
	
	def distance(self, position1, position2):
		x = self.displacement(position1, position2)
		return np.sqrt(x.dot(x))
		
	def derivatives(self, y, t):
		# This method computes the velocity and accelerations of all
		# the asteroids in the system for the odeint function in the
		# results method.

		G = 4*np.pi**2		
		# The value of big G in our current units.
		
		omega = self.omega
		
		D = np.array([]) 
		# D is the array of derivatives
		
		i=0
		for ast in self.asteroids:
			ast.position = y[6*i:6*i+3]
			ast.velocity = y[6*i+3:6*i+6]
			
			gtosun = G*self.sunmass*\
			self.displacement(ast.position,self.sunposition)/\
			self.distance(ast.position,self.sunposition)**3
			# gtosun is the gravitational acceleration of the 
			# asteroid in question towards the sun.
			
			gtoplanet = G*self.planetmass*\
			self.displacement(ast.position,self.planetposition)/\
			self.distance(ast.position,self.planetposition)**3
			# gtoplanet is the gravitational acceleration of the
			# asteroid towards the planet.
		
			centrifugal=-np.cross(omega,np.cross(omega,ast.position))
			# centrifugal acceleration
			
			coriolis = -2*np.cross(omega, ast.velocity)
			# coriolis acceleration
			
			# Since we are in a rotating frame the ficticious forces
			# must be taken into account. These are the centrifugal
			# force and the coriolis force. The Euler force is zero
			# because we make the assumption that omega is constant.
			
			velocity = ast.velocity
			acceleration = gtosun + gtoplanet + centrifugal+coriolis
			
			D = np.concatenate((D, velocity, acceleration))
			# add on the velocity and acceleration of the asteroid
			# to the derivatives array
			i+=1
			
		return D
	
	def results(self):
		# This method computes the results when the system is evolved
		# for the simulation time. 
		
		y0 = np.array([])
		# Initiate an array to store the initial conditions.
		for ast in self.asteroids:
			y0 = np.concatenate((y0, ast.position, ast.velocity))
	
		t = np.linspace(0, self.simulationtime, num=self.N_steps+1)
		
		z = odeint(self.derivatives, y0, t)
		
		result = np.reshape(z, (self.N_steps+1, self.N_asteroids, 2, 3))
		# reshape the array so that the results are easier to interpret.
		# [time, asteroid, pos/vel, x/y/z]
		
		j=0
		for ast in self.asteroids:	
			# update the positions and velocities of the asteroids.
			ast.position = result[-1,j,0,:]
			ast.velocity = result[-1,j,1,:]
			j+=1
		
		return t, result
	
	def wander(self, data):
		# This method generates an array, wandering[i,j], that gives the
		# distance of each asteroid (j) at each point in time (i).
		
		L4 = self.L4		# The coordinates of L4
		L5 = self.L5		# The coordinates of L5
	
		data_shape = list(data.shape)

		newshape = [data_shape[0],data_shape[1],data_shape[3]]
		# At this point we are only interested in the positions, not the
		# velocities, of the asteroids. So, we define a new array shape.
	
		lp_array = np.empty(newshape)	
		# Initiate an array which will store the lagrange point 
		# associated with each asteroid at each point in time.
	
		j = newshape[1]/2
		
		for k in range(3):	
		# Fill the arrays with the relevant lagrange points.
			lp_array[:,:j,k].fill(L4[k])	
			lp_array[:,j:,k].fill(L5[k])
		# half of the asteroids are initially close to L4 and the other 
		# half are close to L5.
	
		displacement = data[:,:,0,:] - lp_array
		# The displacement array gives the displacement of each asteroid
		# from its initial position/the lagrange point at each point in 
		# time. 
		
		wandering = np.sqrt(np.einsum('ijk,ijk->ij', displacement,\
			displacement))
		# use the numpy einsum function to quickly compute the distance 
		#of each asteroid from the lagrange point at each point in time.
	
		# Using the einsum function is of the order of hundreds of times
		# faster than using multiple for loops.

		return wandering
		
class Body(SolarSystem):
	
	def __init__(self, name, position, velocity, mass=0):
		self.name = name
		self.initialposition = None
		self.position = position
		self.velocity = velocity
		self.mass = mass
	
		if mass==0:
			SolarSystem.asteroids.append(self)
		
		SolarSystem.bodies.append(self)
		

