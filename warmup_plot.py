import numpy as np
import matplotlib.pyplot as plt

myfile = open('result.txt')

N = 10

numbers = []

h = np.zeros(N+1)
f_2 = np.zeros(N+1)
f_3 = np.zeros(N+1)
err_2 = np.zeros(N+1)
err_3 = np.zeros(N+1)

for line in myfile:
	line = line.split()
	if line:
		line = [float(i) for i in line]
		numbers.append(line)

for j in range(N+1):
	row = numbers[j]
	h[j] = row[0]
	f_2[j] = row[1]
	f_3[j] = row[2]
	err_2[j] = row[3]
	err_3[j] = row[4]

plt.plot(np.log10(h), np.log10(err_2), 'r', np.log10(h), np.log10(err_3), 'g')
plt.show()

