import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

#for setting background
from scipy.misc import imread
import matplotlib.cbook as cbook


N = 1500
pos1 = np.zeros((N-1,3))
pos2 = np.zeros((N-1,3))

myfile = open('build-SolarSystem-Desktop-Debug/positions.txt')

numbers = []
for line in myfile:
	if line.startswith("C") == False:   
		line = line.split()
		if line:
			line = [float(i) for i in line]
			numbers.append(line)

#j = 2k, 2k+1
for k in range(N-1):
	row1 = numbers[2*k]
	row2 = numbers[2*k+1]

	pos1[k][0] = row1[1]
	pos1[k][1] = row1[2]
	pos1[k][2] = row1[3]

	pos2[k][0] = row2[1]
	pos2[k][1] = row2[2]
	pos2[k][2] = row2[3]


#plot in 3D

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

plt.plot(pos2[:,0], pos2[:,1], pos2[:,2], 'yo', linewidth=50.0)
plt.plot(pos1[:,0], pos1[:,1], pos1[:,2])
plt.xlim([-1.5,1.5])
plt.ylim([-1.5,1.5])
plt.title("Earth and sun's position (Verlet)", size='x-large')
plt.legend(['Earth','Sun'],loc='lower left')

#img = imread('uni.jpg')
#plt.imshow(img, extent=[-0.1,0.1,-0.1,0.1], alpha=0.5)

ax.set_xlabel('x-position')
ax.set_ylabel('y-position')
ax.set_zlabel('z-position')

plt.savefig('fig1_verlet.png')

plt.show()
















