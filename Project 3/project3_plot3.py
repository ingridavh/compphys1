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
pos3 = np.zeros((N-1,3))

myfile = open('build-SolarSystem-Desktop-Debug/positions_jup_1000.txt')

numbers = []
for line in myfile:
	if line.startswith("C") == False:   
		line = line.split()
		if line:
			line = [float(i) for i in line]
			numbers.append(line)

#j = 2k, 2k+1
for k in range(N-1):
	row1 = numbers[3*k]
	row2 = numbers[3*k+1]
	row3 = numbers[3*k+2]

	pos1[k][0] = row1[1]
	pos1[k][1] = row1[2]
	pos1[k][2] = row1[3]

	pos2[k][0] = row2[1]
	pos2[k][1] = row2[2]
	pos2[k][2] = row2[3]

	pos3[k][0] = row3[1]
	pos3[k][1] = row3[2]
	pos3[k][2] = row3[3]


#plot in 3D

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

plt.plot(pos2[:,0], pos2[:,1], pos2[:,2],'w')
plt.plot(pos1[:,0], pos1[:,1], pos1[:,2], 'yo', linewidth=50.0)
plt.plot(pos3[:,0], pos3[:,1], pos3[:,2], '-b')
plt.xlim([-10,10])
plt.ylim([-10,10])
plt.title("Earth, Jupiter and sun's position (Verlet)", size='x-large', color='white')
plt.legend(['Earth','Sun', 'Jupiter'],loc='lower left')

img = imread('uni2.jpg')
plt.imshow(img, extent=[-0.1,0.1,-0.1,0.1], alpha=0.7)

ax.set_xlabel('x-position', color='white')
ax.set_ylabel('y-position',color='white')
ax.set_zlabel('z-position',color='white')

plt.savefig('jup_N1500_1000.png')

plt.show()
















