import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

N = 1000
n=20
pos11 = np.zeros((N-1,3))
pos12 = np.zeros((N-1,3))

pos21 = np.zeros((N-1,3))
pos22 = np.zeros((N-1,3))

pos31 = np.zeros((n-1,3))
pos32 = np.zeros((n-1,3))

myfile1 = open('build-SolarSystem-Desktop-Debug/positions_nasa_verlet_dtE-3.txt')
myfile2 = open('build-SolarSystem-Desktop-Debug/positions_nasa_verlet_dtE-2.txt')
myfile3 = open('build-SolarSystem-Desktop-Debug/positions_nasa_verlet_dtE-1.txt')


numbers1 = []
numbers2 = []
numbers3 = []

for line in myfile1:
	if line.startswith("C") == False:   
		line = line.split()
		if line:
			line = [float(i) for i in line]
			numbers1.append(line)

for line in myfile2:
	if line.startswith("C") == False:   
		line = line.split()
		if line:
			line = [float(i) for i in line]
			numbers2.append(line)

for line in myfile3:
	if line.startswith("C") == False:   
		line = line.split()
		if line:
			line = [float(i) for i in line]
			numbers3.append(line)

#j = 2k, 2k+1
for k in range(N-1):
	row11 = numbers1[2*k]
	row12 = numbers1[2*k+1]

	row21 = numbers2[2*k]
	row22 = numbers2[2*k+1]

	pos11[k][0] = row11[1]
	pos11[k][1] = row11[2]
	pos11[k][2] = row11[3]

	pos12[k][0] = row12[1]
	pos12[k][1] = row12[2]
	pos12[k][2] = row12[3]

	pos21[k][0] = row21[1]
	pos21[k][1] = row21[2]
	pos21[k][2] = row21[3]

	pos22[k][0] = row22[1]
	pos22[k][1] = row22[2]
	pos22[k][2] = row22[3]

for k in range(n-1):
	row31 = numbers3[2*k]
	row32 = numbers3[2*k+1]

	pos31[k][0] = row31[1]
	pos31[k][1] = row31[2]
	pos31[k][2] = row31[3]

	pos32[k][0] = row32[1]
	pos32[k][1] = row32[2]
	pos32[k][2] = row32[3]

plt.plot(pos12[:,0],pos12[:,1], pos22[:,0], pos22[:,1], pos32[:,0], pos32[:,1])
plt.title('Earth position (Verlet)',size='x-large')
plt.legend(['N=1000, dt=0.001', 'N=1000, dt=0.01','N=20, dt=0.1'],loc='lower right')
plt.xlabel('x')
plt.ylabel('y')
plt.savefig('Verlet_stability_dt.png')
plt.show()






















