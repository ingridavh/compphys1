import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#define rho vector
N = 100
rho = np.linspace(0.04, 4, N)

lambda1 = np.zeros(N)
lambda2 = np.zeros(N)
lambda3 = np.zeros(N)



myfile1 = open('build-FYS4150_project2-Desktop-Debug/eigenvec_omega=0.01.txt')
myfile2 = open('build-FYS4150_project2-Desktop-Debug/eigenvec_omega=0.5.txt')
myfile3 = open('build-FYS4150_project2-Desktop-Debug/eigenvec_omega=1.txt')
myfile4 = open('build-FYS4150_project2-Desktop-Debug/eigenvec_omega=5.txt')


#Choose myfile to we wanted file

myfile = myfile1


numbers = []
for line in myfile:
    line = line.split()
    if line:
        line = [float(i) for i in line]
        numbers.append(line)

for j in range(N):
    j = int(j)

    row = numbers[j]
    lambda1[j] = row[0]
    lambda2[j] = row[1]
    lambda3[j] = row[2]


plt.plot(rho, lambda1, rho, lambda2, rho, lambda3)
plt.show()
