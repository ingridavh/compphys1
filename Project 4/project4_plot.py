import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#myfile = open('build-FYS4150_project4-Desktop-Debug/prettywoman.txt')
#myfile2 = open('build-FYS4150_project4-Desktop-Debug/prettywoman_.txt')
numbers = []
numbers2 = []

for line in myfile:
	line = line.split()
	if line:
		for i in line:
			if i != ",":
				word = float(i)
				numbers.append(word)
N = len(numbers)
n = N/4

mcs = np.zeros(n)
E = np.zeros(n)
M = np.zeros(n)
no_accept = np.zeros(n)

for k in range(n):
	mcs[k] = numbers[4*k]
	E[k] = numbers[4*k+1]
	M[k] = numbers[4*k+2]
	no_accept[k] = numbers[4*k+3]

for line in myfile2:
	line = line.split()
	if line:
		for i in line:
			if i != ",":
				word = float(i)
				numbers2.append(word)

mcs2 = np.zeros(n)
E2 = np.zeros(n)
M2 = np.zeros(n)
no_accept2 = np.zeros(n)

for k in range(n):
	mcs2[k] = numbers2[4*k]
	E2[k] = numbers2[4*k+1]
	M2[k] = numbers2[4*k+2]
	no_accept2[k] = numbers2[4*k+3]




"""
plt.subplot(2,1,1)
plt.plot(mcs[0:3], E[0:3])
plt.ylabel('Expectation value')
plt.legend(['Energy'], loc='lower right')
plt.title('Expectation values')

plt.subplot(2,1,2)
plt.plot(mcs[4:n], E[4:n])
plt.xlabel('Number of Monte Carlo cycles')
plt.ylabel('Expectation value')
plt.legend(['Energy'], loc='lower right')
plt.savefig('exp_E_T24_ordered_.pdf')
plt.show()
"""
"""
#plt.semilogx(mcs, no_accept)
plt.plot(mcs, no_accept)
plt.xlabel('Temperature')
plt.ylabel('Number of accepted configurations')
plt.title('Accepted configurations as a function of cycles')
plt.savefig('accept_cycl_t1.pdf')
plt.show()
"""

plt.subplot(2,1,1)
plt.title('P(E)')
plt.hist(E, bins=100)
plt.legend(['T=1.0'])

plt.subplot(2,1,2)
plt.hist(E2, bins=100, facecolor='red')
plt.legend(['T=2.4'])
plt.xlabel('Energy')
plt.ylabel('Frequency')
plt.savefig('P(E).pdf')
plt.show()










