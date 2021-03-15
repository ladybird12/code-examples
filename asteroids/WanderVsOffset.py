# WanderVsOffset.py

from SolarSystemRotatingFrame import *


figname = 'AsteroidWanderVsOffset7'# specify a filename to save the graph to

Jmass = 0.001		# The mass of Jupiter in units of the solar mass
R = 5.2			# The mean distance between the Sun and Jupiter in AU
dt = 0.01		# Time step in Earth years
simulationtime = 1200	# Total simulation time in Earth years
N_asteroids = 2		# Total number of asteroids, half at each L.P.
randomamplitude = 0

offsetlist = np.linspace(-0.01,0.01, num=100)
# Create an array for the range of offsets to test. The offset is the magnitude
# of the initial displacement of the asteroids from the L.P.s, in AU.

asteroidoffsets = np.array([[r, 0, 0] for r in offsetlist])

systems = [SolarSystem(Jmass, R, dt, simulationtime, N_asteroids, \
	offset, randomamplitude) for offset in asteroidoffsets]
	# Set up the systems with the different planet masses

meanwander = []
wandererror = []

number = 0
for system in systems:
	
	system.setup()
	t, results = system.results()		
	wandering = system.wander(results)
	
	asteroidwander =[max(wandering[:,j]) for j in range(system.N_asteroids)]
	meanwander.append(np.mean(asteroidwander))
	wandererror.append(np.std(asteroidwander))

	print("Number ",number+1," completed.")
	
	SolarSystem.asteroids = []	# reset asteroids
	SolarSystem.bodies = []		# reset bodies
	number+=1

rawdata, = plt.plot(offsetlist, meanwander, 'bo')
	
plt.errorbar(offsetlist, meanwander, yerr=wandererror, fmt=None)

#plt.xscale('log')
#plt.yscale('log')
plt.xlabel(r'Time / years')
plt.ylabel('Distance from Lagrange point / AU')
plt.title('The relationship between the computed \n wander and the initial \
offset', fontsize=18)

plt.savefig(figname, format='png', dpi=1000)
plt.show()
