# WanderVsMass.py

from SolarSystemRotatingFrame import *


figname = 'AsteroidWanderVsPlanetMass8'# specify a filname to save the graph to

R = 5.2			# The mean distance between the Sun and Jupiter in AU
dt = 0.01		# Time step in Earth years
simulationtime = 1200	# Total simulation time in Earth years
N_asteroids = 2		# Total number of asteroids, half at each L.P.
offset = [[0.01, 0, 0], [0, 0.01, 0], [0, 0, 0.01]]	
			# The change in the initial position of the asteroids 
			# from the LP, in the rho, phi, z directions, in AU.
randomamplitude = 0	# The magnitude of the random vector also added to the
			# asteroid initial position, in AU.
colours = ['bo','g.','r+'] 
# Give the plots for each of the offsets different colours.
			
planetmasses = np.linspace(0.0001, 0.04, num=100)
# Create an array for the range of planet masses to test. The units are solar 
# masses. These planets are all going to be at Jupiter's radius from the sun.

datalist = []

for i in range(3):

	systems = [SolarSystem(mass, R, dt, simulationtime, N_asteroids, offset[i], \
randomamplitude) for mass in planetmasses]
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
							       
		print(system.planetmass, "Number ",number+1," completed.")
	
		SolarSystem.asteroids = []	# reset asteroids
		SolarSystem.bodies = []		# reset bodies
		number+=1
	
	rawdata, = plt.plot(planetmasses, meanwander/min(meanwander), colours[i])

	plt.errorbar(planetmasses, meanwander/min(meanwander), yerr=wandererror, fmt=None)

	datalist.append(rawdata)
	
plt.legend(datalist, [r'$\rho$ displacement', r'$\phi$ displacement', \
r"$\mathit{z}$ displacement"], loc=2)
plt.xlim([0,0.045]), plt.ylim([0,10])
plt.xlabel('Planet Mass / Solar masses')
plt.ylabel('Wander / Minimum Wander')
plt.title('The relationship between the computed \n wander and the planet mass',
	fontsize=18)

plt.savefig(figname, format='png', dpi=1000)
plt.show()
