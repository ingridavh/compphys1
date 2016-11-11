import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

myfile = open('build-FYS4150_project4-Desktop-Debug/exp_mcs.txt')
numbers = []

for line in myfile:
	line = line.split()
	if line:
		for i in line:
			if i != ",":
				word = float(i)
				numbers.append(word)
N = len(numbers)

mcs = np.zeros(N/3)
E = np.zeros(N/3)
M = np.zeros(N/3)

for k in range(N/3):
	mcs[k] = numbers[3*k]
	E[k] = numbers[3*k+1]
	M[k] = numbers[3*k+2]

mcs = [np.log10(mcs[i]) for i in range(N/3)] 

plt.plot(mcs, E, mcs, M)
plt.xlabel('Number of Monte Carlo cycles [log(x)]')
plt.ylabel('Expectation value [log(y)]')
plt.title('Expectation values for Monte Carlo')
plt.legend(['Energy', 'Magnetization'])
plt.show()
